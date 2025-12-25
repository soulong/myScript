#! /home/haohe/miniforge3/bin/conda run -n cellprofiler --no-capture-output python
#!/usr/bin/env python
"""
run_cellprofiler.py
-------------------
Run CellProfiler in parallel over all images in every ``*_Measurement X`` folder.

* Discovers measurements automatically
* Splits dataloader (CSV) into batches → one CellProfiler process per batch
* Each batch logs to ``<root>/logs/cellprofiler/<measurement>/batch_XXX.log``
* Merges per‑batch SQLite DBs into a single ``cp_result.db``
* Fixes ``Experiment_Properties`` and copies ``.properties`` files
* Full CLI with defaults (``--help``)
"""

from __future__ import annotations

import argparse
import math
import sys
import time
import logging
from datetime import datetime
from pathlib import Path
from typing import List, Optional
from natsort import natsorted
import pandas as pd
from sqlalchemy import create_engine
from subprocess import run
from concurrent.futures import ProcessPoolExecutor, as_completed

from image_helper import setup_logger, find_measurement_dirs


# ----------------------------------------------------------------------- #
# Run one batch of CellProfiler
# ----------------------------------------------------------------------- #
def _run_cellprofiler_batch(
    dataloader_path: Path,
    cp_project_path: Path,
    cp_exec_path: Path,
    output_dir: Path,
    log_path: Path,
) -> int:
    """
    Execute CellProfiler on a **single** dataloader CSV.

    Parameters
    ----------
    dataloader_path
        Path to batch CSV.
    cp_project_path
        Path to ``.cppipe`` or ``.cpproj``.
    cp_exec_path
        Path to ``cellprofiler`` executable.
    output_dir
        Where to write DB and log.
    log_path
        Full path to log file.

    Returns
    -------
    int
        Process exit code.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Activate conda env + run CellProfiler
    cmd = (
        f'. ~/miniforge3/etc/profile.d/conda.sh && '
        f'conda activate cellprofiler && '
        f'conda info && '
        f'"{cp_exec_path}" -c -r '
        f'-p "{cp_project_path}" '
        f'-o "{output_dir}" '
        f'--data-file "{dataloader_path}"'
    )

    start = datetime.now()
    print(f"[{start}] START batch → {dataloader_path.name} | log: {log_path.name}")

    with log_path.open("w", encoding="utf-8", buffering=1) as log_file:
        log_file.write(f"=== CellProfiler START: {start} ===\n")
        log_file.write(f"CMD: {cmd}\n\n")
        log_file.flush()

        try:
            result = run(
                cmd,
                shell=True,
                text=True,
                stdout=log_file,
                stderr=log_file,
            )
            end = datetime.now()
            duration = end - start
            status = "SUCCESS" if result.returncode == 0 else f"FAILED (code {result.returncode})"
            log_file.write(f"\n=== CellProfiler END: {end} | Duration: {duration} | {status} ===\n")
            print(f"[{end}] END batch → {status}")
            return result.returncode
        except Exception as e:
            error = f"EXCEPTION: {e}\n"
            log_file.write(error)
            print(f"[ERROR] {error}")
            return -1



# ----------------------------------------------------------------------- #
# DB merge helper
# ----------------------------------------------------------------------- #
def _merge_cellprofiler_dbs(
    job_dir: Path,
    merged_db_path: Path,
    measurement_dir: Path,
    logger: logging.Logger,
) -> None:
    """Merge per‑batch SQLite DBs into a single DB with correct ImageNumber offsets."""
    start = time.time()

    # Create merged DB
    merged_engine = create_engine(f"sqlite:///{merged_db_path}", connect_args={"uri": False})
    merged_conn = merged_engine.connect()
    # Find all batch DBs
    db_paths = natsorted(job_dir.rglob("*.db"))
    if not db_paths:
        raise FileNotFoundError("No batch DBs found")

    image_offset = 0
    object_tables: Optional[List[str]] = None

    for idx, db_path in enumerate(db_paths):
        engine = create_engine(f"sqlite:///{db_path}", connect_args={"uri": False})
        conn = engine.connect()

        # First DB: copy metadata tables
        if idx == 0:
            for table in ["Experiment", "Experiment_Properties", "Per_Experiment"]:
                df = pd.read_sql(table, con=conn)
                df.to_sql(table, merged_engine, if_exists="replace", index=False)

            # Discover object tables
            import sqlite3
            tmp = sqlite3.connect(str(db_path))
            cursor = tmp.cursor()
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            all_tables = [row[0] for row in cursor.fetchall()]
            object_tables = [
                t for t in all_tables
                if t not in {"Per_Image", "Experiment", "sqlite_sequence", "Experiment_Properties", "Per_Experiment"}
            ]
            tmp.close()

        # Append Per_Image with offset
        per_image = pd.read_sql("Per_Image", con=conn)
        per_image["ImageNumber"] += image_offset
        per_image.to_sql("Per_Image", merged_engine, if_exists="append", index=False)

        # Append object tables
        for obj in object_tables:
            df = pd.read_sql(obj, con=conn)
            df["ImageNumber"] += image_offset
            df.to_sql(obj, merged_engine, if_exists="append", index=False)

        image_offset += len(per_image)
        conn.close()

    # Fix Experiment_Properties → db path
    exp_prop = pd.read_sql("Experiment_Properties", con=merged_engine)
    # print(exp_prop)
    exp_prop.loc[exp_prop["field"] == "db_sqlite_file", "value"] = str(merged_db_path)
    exp_prop.to_sql("Experiment_Properties", merged_engine, if_exists="replace", index=False)

    # Copy and fix .properties files
    first_job_dir = job_dir / "job_000"
    for obj in object_tables:
        src = first_job_dir / f"cp_result_{obj.replace('Per_', '')}.properties"
        dst = measurement_dir / f"cp_result_{obj.replace('Per_', '')}.properties"
        if src.exists():
            dst.write_text(src.read_text().replace("jobs/job_000/", ""), encoding="utf-8")

    logger.info(f"DB merge took {time.time() - start:.1f}s")




# ----------------------------------------------------------------------- #
# Process one measurement folder
# ----------------------------------------------------------------------- #
def cellprofiler_process_measurement(
    measurement_dir: Path,
    cp_project_path: Path,
    cp_exec_path: Path,
    max_workers: int,
    overwrite: bool,
    test: bool,
    logger: logging.Logger,
) -> bool:
    """
    Run CellProfiler on **one** measurement.

    Returns
    -------
    bool
        ``True`` if all batches succeeded and DB merged.
    """
    dataloader_path = measurement_dir / "cp_dataloader.csv"
    if not dataloader_path.exists():
        logger.warning(f"cp_dataloader.csv not found in {measurement_dir.name} → skip")
        return True

    merged_db_path = measurement_dir / "cp_result.db"
    if merged_db_path.exists() and not overwrite:
        logger.info(f"cp_result.db already exists → skip {measurement_dir.name}")
        return True

    # ------------------------------------------------------------------- #
    # Load and optionally filter dataloader
    # ------------------------------------------------------------------- #
    df = pd.read_csv(dataloader_path)
    if df.empty:
        logger.info(f"Empty dataloader → skip {measurement_dir.name}")
        return True

    if test:
        df = df.sample(n=min(100, len(df)), random_state=42)
        logger.info(f"TEST MODE: using {len(df)} rows")

    # ------------------------------------------------------------------- #
    # Split into batches
    # ------------------------------------------------------------------- #
    batch_size = math.ceil(len(df) / max_workers)
    batches: List[pd.DataFrame] = [
        df.iloc[i:i + batch_size] for i in range(0, len(df), batch_size)
    ]
    logger.info(f"Splitting {len(df)} rows into {len(batches)} batches (size ≈ {batch_size})")

    # Save batch CSVs
    job_dir = measurement_dir / "jobs"
    batch_paths: List[Path] = []
    for idx, batch_df in enumerate(batches):
        batch_dir = job_dir / f"job_{idx:03d}"
        batch_csv = batch_dir / f"job_{idx:03d}.csv"
        batch_dir.mkdir(parents=True, exist_ok=True)
        batch_df.to_csv(batch_csv, index=False)
        batch_paths.append(batch_csv)

    # ------------------------------------------------------------------- #
    # Run in parallel
    # ------------------------------------------------------------------- #
    log_dir = measurement_dir / "logs" / "cellprofiler"
    failed = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_idx = {
            executor.submit(
                _run_cellprofiler_batch,
                dataloader_path=batch_csv,
                cp_project_path=cp_project_path,
                cp_exec_path=cp_exec_path,
                output_dir=batch_csv.parent,
                log_path=log_dir / f"batch_{idx:03d}.log",
            ): idx
            for idx, batch_csv in enumerate(batch_paths)
        }
        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            code = future.result()
            if code != 0:
                failed.append(idx)
    if failed:
        logger.error(f"Failed batches in {measurement_dir.name}: {failed}")
        return False

    # ------------------------------------------------------------------- #
    # Merge SQLite DBs
    # ------------------------------------------------------------------- #
    logger.info(f"Merging {len(batches)} DBs → {merged_db_path.name}")
    try:
        _merge_cellprofiler_dbs(
            job_dir=job_dir,
            merged_db_path=merged_db_path,
            measurement_dir=measurement_dir,
            logger=logger,
        )
        logger.info("DB merge complete")
    except Exception as e:
        logger.error(f"DB merge failed: {e}, remove existed DB file", exc_info=True)
        # delete already generated DB file if failed
        try:
            logger.info(f"remove DB file", exc_info=True)
            merged_db_path.unlink()
            print(f"File '{merged_db_path}' removed successfully.")
        except Exception as e:
            print(f"An error occurred: {e}")
        return False

    return True



# ----------------------------------------------------------------------- #
# CLI
# ----------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(description="Run CellProfiler in parallel over multiple Measurement folders", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("root_dir", type=str, help="Root folder containing *_Measurement X folders")
    parser.add_argument("--cp-project", type=str, required=True, help="Path to CellProfiler pipeline (.cppipe or .cpproj)")
    parser.add_argument("--cp-exec", type=str, default="/home/haohe/miniforge3/envs/cellprofiler/bin/cellprofiler", help="Path to CellProfiler executable")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel CellProfiler processes")
    parser.add_argument("--overwrite", action="store_true", help="Reprocess even if cp_result.db exists")
    parser.add_argument("--test", action="store_true", help="Process only ~100 random rows per measurement")
    args = parser.parse_args()

    # ------------------------------------------------------------------- #
    # Setup
    # ------------------------------------------------------------------- #
    root_dir = Path(args.root_dir)
    if not root_dir.is_dir():
        print(f"Root directory not found: {root_dir}")
        sys.exit(1)

    cp_project = Path(args.cp_project)
    if not cp_project.exists():
        print(f"CellProfiler project not found: {cp_project}")
        sys.exit(1)

    cp_exec = Path(args.cp_exec)
    if not cp_exec.exists():
        print(f"CellProfiler executable not found: {cp_exec}")
        sys.exit(1)

    # Logging
    log_dir = root_dir / "logs"
    logger = setup_logger(log_dir, name="cellprofiler")
    logger.info(f"CellProfiler pipeline started – root: {root_dir}")
    logger.info(f"Project: {cp_project} | Exec: {cp_exec} | Workers: {args.workers}")

    # Discover measurements
    measurements = find_measurement_dirs(root_dir)
    if not measurements:
        logger.error("No *_Measurement X folders found")
        sys.exit(1)
    logger.info(f"Found {len(measurements)} measurement directories")

    # ------------------------------------------------------------------- #
    # Run per measurement
    # ------------------------------------------------------------------- #
    overall_success = True
    global_start = time.time()

    for mdir in measurements:
        logger.info(f"--- Processing {mdir.name} ---")
        try:
            ok = cellprofiler_process_measurement(
                measurement_dir=mdir,
                cp_project_path=cp_project,
                cp_exec_path=cp_exec,
                max_workers=args.workers,
                overwrite=args.overwrite,
                test=args.test,
                logger=logger,
            )
            if not ok:
                overall_success = False
        except Exception as e:
            logger.error(f"Unexpected error in {mdir.name}: {e}", exc_info=True)
            overall_success = False

    total_time = time.time() - global_start
    logger.info(f"CellProfiler pipeline finished in {total_time:.1f}s")
    sys.exit(0 if overall_success else 1)


if __name__ == "__main__":
    main()