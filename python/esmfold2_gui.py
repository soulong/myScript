#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
esmfold2_gui.py — ESMFold2 web GUI for all-atom complex structure prediction.

conda env : bio  (win11, python 3.12, RTX 5070 Ti 16GB, torch cu128+)
required  : esm (Biohub), gradio, pandas, matplotlib, rdkit

Run:
    conda activate bio
    python esmfold2_gui.py          # opens the web page in the browser

Features
    * Fold single or multiple complexes of protein / DNA / RNA / small molecule
      (SMILES) / modified residues (CCD codes) with Biohub ESMFold2.
    * Two checkpoints:
        - ESMFold2-Fast (default) : inference-optimized single-sequence model,
          several times faster than ESMFold2.
        - ESMFold2 (optional)     : MSA-capable model; upload .a3m MSAs to
          condition protein/RNA chains.
      Both share the ESMC-6B backbone weights.
    * Inputs are dynamic molecule rows (AlphaFold-3-server style):
        - Single complex : rows are added with "+ Add molecule"; each row has
          its own type selector (`protein`/`dna`/`rna`/`ligand`) and one text
          field that holds AA/bases, or the SMILES for ligands. All rows are
          folded as ONE complex.
        - N x N combination : two boxes (List A / List B), each with its own
          molecule type; every component of A is crossed with every component
          of B.
    * Metrics: native pLDDT / pTM / ipTM / PAE from the model, PAE heatmap,
      per-chain mean pLDDT bars, per-residue pLDDT CSV.
      (reference-free pDockQ via the DockQ package is temporarily skipped)
    * Results shown interactively (3Dmol.js viewer, colored by pLDDT) and
      downloadable per job (zip).

Memory notes
    * The 6B model needs ~13 GB RAM + ~13 GB VRAM. Close other heavy apps
      before folding, otherwise loading swaps and looks like a hang.
    * The ~400 MB CCD pickle (needed only for ligands / modified residues) is
      loaded before the model and SKIPPED entirely for ligand-free jobs.
    * Stay on the officially supported dtype mix: ESMC backbone bf16 (default
      `esmc_precision="bf16"`), folding head fp32. Do NOT cast the whole model
      to bf16 via `dtype=` — the processor emits fp32 inputs and a bf16 head
      raises "mat1 and mat2 must have the same dtype".
"""

import os
import re
import json
import zipfile
import datetime
import urllib.request
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gradio as gr

# === CONFIG ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
STATIC_DIR = os.path.join(BASE_DIR, "3dmol_static")
MAX_TOKENS = 2048          # hard length cap (tokens), safe on 16 GB VRAM
VALID_TYPES = ["protein", "dna", "rna", "ligand"]
MAX_SINGLE = 6             # max molecules in the single-complex list
MAX_COMBO = 4              # max components per box in N x N mode
MAIN_REPO = "biohub/ESMFold2"       # MSA-capable checkpoint (owns ccd.pkl)
DEFAULT_MODEL = "ESMFold2-Fast (default, single-sequence)"
CHECKPOINTS = {
    "ESMFold2-Fast (default, single-sequence)": "biohub/ESMFold2-Fast",
    "ESMFold2 (MSA-capable)": MAIN_REPO,
}
# The huggingface.co endpoints can be unreachable (CN networks); weights may
# instead be staged locally from a mirror into the HF cache layout (see
# esmfold2_runs/weights_note.txt). When a staged snapshot exists it is used
# directly, otherwise from_pretrained falls back to the network.
HF_HUB_CACHE = os.path.join(os.environ.get("USERPROFILE") or os.path.expanduser("~"),
                            ".cache", "huggingface", "hub")
SNAP_REV = "modelscope-mirror"

# === 0) shared state ===
STOP_FLAG = {"stop": False}
MODEL_CACHE = {}           # checkpoint id -> loaded model
BUILDER_CACHE = {}         # "ccd"/"noccd" -> ESMFold2InputBuilder
RESULTS = []               # list of job dicts, refreshed from disk on demand
_3DMOL_JS = None           # cached JS path


# === 1) input parsing / validation ===
def _make_row(group, typ, txt, ccd, name, copies, default_chain):
    """One molecule row -> plain dict; None when the row is empty."""
    typ = (typ or "protein").strip().lower()
    txt = (txt or "").strip().replace(" ", "").replace("\n", "")
    ccd = (ccd or "").strip()
    name = (name or "").strip()
    try:
        copies = max(1, int(copies or 1))
    except (ValueError, TypeError):
        copies = 1
    if not (txt or ccd):
        return None
    chain = name or default_chain
    is_lig = typ == "ligand"
    return {"group": group, "component": chain, "chain": chain, "type": typ,
            "sequence": "" if is_lig else txt, "copies": copies,
            "ccd": ccd, "smiles": txt if is_lig else ""}


def _slots(flat, per, n):
    """Split a flat component-value list into per-slot tuples."""
    return [tuple(flat[k * per:(k + 1) * per]) for k in range(n)]


def parse_single_slots(slots, flags):
    """Single-complex mode: (type, text, ccd, name, copies) rows -> row dicts."""
    rows = []
    flags = flags or []
    for i, tup in enumerate(slots):
        if i >= len(flags) or not flags[i] or tup is None:
            continue
        if all(v is None or str(v).strip() == "" for v in tup):
            continue
        row = _make_row("A", tup[0], tup[1], tup[2], tup[3], tup[4],
                        f"chain{i + 1}")
        if row:
            rows.append(row)
    return rows


def parse_combo_slots(group, box_type, slots, flags):
    """One N x N box: (text, ccd, name, copies) rows -> row dicts."""
    rows = []
    flags = flags or []
    for i, tup in enumerate(slots):
        if i >= len(flags) or not flags[i] or tup is None:
            continue
        if all(v is None or str(v).strip() == "" for v in tup):
            continue
        row = _make_row(group, box_type, tup[0], tup[1], tup[2], tup[3],
                        f"{group}mol{i + 1}")
        if row:
            rows.append(row)
    return rows


def collect_rows(mode, single_slots, single_flags,
                 type_a, a_slots, a_flags, type_b, b_slots, b_flags):
    """Pick the row builder active for the current mode."""
    if mode == "Single complex":
        return parse_single_slots(single_slots, single_flags)
    return (parse_combo_slots("A", type_a, a_slots, a_flags) +
            parse_combo_slots("B", type_b, b_slots, b_flags))


def needs_ccd(rows):
    """True when any chain needs the CCD dictionary (ligands / modifications)."""
    return any(r["type"] == "ligand" or r["ccd"] for r in rows)


def parse_modifications(ccd_str):
    """Parse 'SEP@5,PHO@10' into a list of (ccd_code, 0-based position)."""
    mods = []
    for item in ccd_str.split(","):
        item = item.strip()
        if not item:
            continue
        if "@" in item:
            code, pos = item.rsplit("@", 1)
            mods.append((code.strip(), int(pos) - 1))    # 1-based in the GUI
        else:
            mods.append((item, -1))                       # position unknown
    return mods


def estimate_tokens(row):
    """Token estimate: residues/bases for polymers, heavy atoms for ligands."""
    if row["type"] == "ligand":
        n = 1
        if row["smiles"]:
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(row["smiles"])
                if mol is not None:
                    n = mol.GetNumHeavyAtoms()
            except Exception:
                pass
        return n * row["copies"]
    return len(row["sequence"]) * row["copies"]


def validate_rows(rows, mode):
    """Raise ValueError with a clear message on invalid input."""
    if not rows:
        raise ValueError("No molecules defined — add them with '+ Add molecule'.")
    names = [r["chain"] for r in rows]
    if len(names) != len(set(names)):
        dups = sorted({c for c in names if names.count(c) > 1})
        raise ValueError(f"Duplicate chain names: {dups} — give unique names.")
    bad_type = [r["chain"] for r in rows if r["type"] not in VALID_TYPES]
    if bad_type:
        raise ValueError(f"Unknown molecule type for {bad_type} (valid: {VALID_TYPES}).")
    for r in rows:
        if r["type"] in ("protein", "dna", "rna"):
            if not r["sequence"]:
                raise ValueError(f"Chain '{r['chain']}' needs a sequence.")
            if r["smiles"]:
                raise ValueError(f"Chain '{r['chain']}' is {r['type']}; remove the SMILES.")
        if r["type"] == "ligand":
            if not (r["smiles"] or r["ccd"]):
                raise ValueError(f"Ligand chain '{r['chain']}' needs a SMILES or CCD code.")
            if r["sequence"]:
                raise ValueError(f"Ligand chain '{r['chain']}' must not have a sequence.")
    total = sum(estimate_tokens(r) for r in rows)
    if total > MAX_TOKENS:
        raise ValueError(
            f"Complex too long: ~{total} tokens > {MAX_TOKENS} cap. "
            "Shorten sequences or reduce copies."
        )
    if mode == "N x N combination":
        groups = sorted({r["group"] for r in rows})
        if len(groups) != 2 or any(not any(x["group"] == g for x in rows) for g in groups):
            raise ValueError(
                f"Combination mode needs BOTH input boxes non-empty "
                f"(found groups {groups})."
            )


def build_jobs(rows, mode):
    """Expand rows into the list of jobs to fold.

    Returns a list of jobs; each job is a dict:
        name   : str, unique job name
        chains : list of row dicts included in this complex
    """
    jobs = []
    if mode == "N x N combination":
        groups = sorted({r["group"] for r in rows}, key=str.casefold)
        g0, g1 = groups
        comps0 = sorted({r["component"] for r in rows if r["group"] == g0},
                        key=str.casefold)
        comps1 = sorted({r["component"] for r in rows if r["group"] == g1},
                        key=str.casefold)
        for c0 in comps0:
            for c1 in comps1:
                sel = [r for r in rows
                       if (r["group"] == g0 and r["component"] == c0)
                       or (r["group"] == g1 and r["component"] == c1)]
                jobs.append({"name": _slug(c0 + "__" + c1), "chains": sel})
    else:
        jobs.append({"name": "complex", "chains": rows})
    # make job names unique
    seen = {}
    for job in jobs:
        base = job["name"]
        seen[base] = seen.get(base, 0) + 1
        job["name"] = base if seen[base] == 1 else f"{base}_{seen[base]}"
    return jobs


def _slug(s):
    s = re.sub(r"[^0-9A-Za-z_\-]+", "_", s).strip("_")
    return s or "job"


def job_tokens(job):
    return sum(estimate_tokens(r) for r in job["chains"])


def load_msa_files(files):
    """Map file stem -> local path for uploaded MSA files."""
    msa_map = {}
    for f in files or []:
        if hasattr(f, "name"):
            path = f.name
            stem = os.path.splitext(os.path.basename(path))[0]
            msa_map[stem] = path
    return msa_map


# === 2) sequence building (ESMFold2 input) ===
def build_structure_input(job, msa_map, use_msa):
    """Build an esm StructurePredictionInput for one job."""
    from esm.models.esmfold2 import (
        ProteinInput, RNAInput, DNAInput, LigandInput, Modification,
        StructurePredictionInput,
    )

    sequences = []
    for r in job["chains"]:
        ccd = r["ccd"]
        mods = None
        if r["type"] in ("protein", "dna", "rna") and ccd:
            parsed = parse_modifications(ccd)
            mods = [Modification(position=pos, ccd=code)
                    for code, pos in parsed if pos >= 0]
        ids = [f"{r['chain']}_{i}" for i in range(1, r["copies"] + 1)]

        if r["type"] == "protein":
            msa = None
            if use_msa and r["chain"] in msa_map:
                msa = build_msa(msa_map[r["chain"]])
            sequences.append(ProteinInput(id=ids, sequence=r["sequence"],
                                          modifications=mods, msa=msa))
        elif r["type"] == "rna":
            msa = None
            if use_msa and r["chain"] in msa_map:
                msa = build_msa(msa_map[r["chain"]])
            sequences.append(RNAInput(id=ids, sequence=r["sequence"],
                                      modifications=mods, msa=msa))
        elif r["type"] == "dna":
            sequences.append(DNAInput(id=ids, sequence=r["sequence"],
                                      modifications=mods))
        else:  # ligand
            if ccd:
                sequences.append(LigandInput(id=ids, ccd=[c.strip() for c in ccd.split(",")]))
            else:
                sequences.append(LigandInput(id=ids, smiles=r["smiles"]))

    return StructurePredictionInput(sequences=sequences)


def build_msa(a3m_path):
    """Build an esm MSA object from a .a3m file."""
    from esm.utils.msa import MSA
    return MSA.from_a3m(a3m_path)


def local_repo_dir(repo_id):
    """Path of a locally-staged snapshot for an HF repo id (mirror layout)."""
    return os.path.join(HF_HUB_CACHE, "models--" + repo_id.replace("/", "--"),
                        "snapshots", SNAP_REV)


def get_model(checkpoint):
    """Load (and cache) the ESMFold2 model for a checkpoint id.

    Uses the locally staged snapshot when present (huggingface.co may be
    unreachable), otherwise falls back to from_pretrained over the network.
    Keeps the officially supported dtype mix: ESMC backbone bf16 (default
    esmc_precision="bf16"), folding head fp32 (do not pass dtype=).
    """
    if checkpoint not in MODEL_CACHE:
        from esm.models.esmfold2 import EsmFold2Model, EsmFold2Config
        from esm.models.esmfold2.config import default_module_flags
        local = local_repo_dir(checkpoint)
        if os.path.isdir(local) and os.path.isfile(os.path.join(local, "config.json")):
            print(f"Loading {checkpoint} from local mirror snapshot ...", flush=True)
            cfg = EsmFold2Config.from_pretrained(local, **default_module_flags(local))
            esmc_local = local_repo_dir(cfg.esmc_id)
            if os.path.isdir(esmc_local):
                cfg.esmc_id = esmc_local      # ESMC backbone lives next to it
            model = EsmFold2Model.from_pretrained(local, config=cfg,
                                                  device="cuda").eval()
        else:
            print(f"Loading {checkpoint} ...", flush=True)
            model = EsmFold2Model.from_pretrained(checkpoint, device="cuda").eval()
        MODEL_CACHE[checkpoint] = model
    return MODEL_CACHE[checkpoint]


def get_builder(need_ccd):
    """Return a cached ESMFold2InputBuilder.

    Its constructor loads the ~400 MB CCD pickle, which is only required when
    the batch contains ligands / modified residues. CCS lives in the MSA model
    repo (biohub/ESMFold2); ligand-free batches use a builder created with an
    empty CCD dict (fast, no 400 MB load / no RAM spike); the two cached
    builders coexist safely.
    """
    key = "ccd" if need_ccd else "noccd"
    if key not in BUILDER_CACHE:
        from esm.models.esmfold2 import ESMFold2InputBuilder
        from esm.models.esmfold2 import conformers as _cf
        main_local = local_repo_dir(MAIN_REPO)
        if need_ccd and os.path.isfile(os.path.join(main_local, "ccd.pkl")):
            print("Loading CCD dictionary ...", flush=True)
            BUILDER_CACHE[key] = ESMFold2InputBuilder(ccd_cache=main_local)
        else:
            # empty CCD placeholder - only safe when nothing needs CCD
            saved = _cf._CCD_MOLECULES
            _cf._CCD_MOLECULES = {}
            try:
                BUILDER_CACHE[key] = ESMFold2InputBuilder()
            finally:
                _cf._CCD_MOLECULES = saved
    return BUILDER_CACHE[key]


def fold_job(job, checkpoint, num_loops, num_sampling_steps, seed=0,
             msa_map=None, use_msa=False):
    """Run a fold for one job; returns (fold result, chain id list)."""
    model = get_model(checkpoint)
    spi = build_structure_input(job, msa_map, use_msa)
    builder = get_builder(needs_ccd(job["chains"]))
    result = builder.fold(
        model, spi,
        num_loops=int(num_loops),
        num_sampling_steps=int(num_sampling_steps),
        num_diffusion_samples=1,
        seed=seed,
    )
    chain_ids = [f"{r['chain']}_{i}" for r in job["chains"]
                 for i in range(1, r["copies"] + 1)]
    return result, chain_ids


def chain_token_counts(job):
    """Per-chain token counts in the same order as the chain ids."""
    counts = []
    for r in job["chains"]:
        n = 1 if r["type"] == "ligand" else len(r["sequence"])
        counts.extend([n] * r["copies"])
    return counts


# === 3) metrics / artifact writing ===
def write_job_artifacts(job_dir, result, chain_ids, token_counts_):
    """Dump structure + metrics + plots for one job. Returns a summary dict."""
    os.makedirs(job_dir, exist_ok=True)

    # structure (mmCIF)
    cif_path = os.path.join(job_dir, "complex.cif")
    with open(cif_path, "w") as f:
        f.write(result.complex.to_mmcif())

    # per-token pLDDT (the model returns 0..1 floats; store everything 0..100)
    plddt = np.asarray(result.plddt, dtype=float).reshape(-1)
    if plddt.max() <= 1.0:
        plddt = plddt * 100.0
    n_total = plddt.shape[0]

    # mean pLDDT per chain (polymers: 1 token/residue; ligand: 1 token)
    chain_plddt = {}
    start = 0
    for cid, n in zip(chain_ids, token_counts_):
        seg = plddt[start:start + n]
        chain_plddt[cid] = float(seg.mean()) if seg.size else float("nan")
        start += n
    mean_plddt = float(plddt.mean()) if plddt.size else float("nan")

    # per-residue pLDDT CSV (token-level, includes the dict of chains order)
    pd.DataFrame([{"index": i, "plddt": float(v)} for i, v in enumerate(plddt)]
                 ).to_csv(os.path.join(job_dir, "plddt.csv"), index=False)

    pd.DataFrame(
        [{"chain": cid, "n_tokens": n, "mean_plddt": chain_plddt[cid]}
         for cid, n in zip(chain_ids, token_counts_)]
    ).to_csv(os.path.join(job_dir, "chain_plddt.csv"), index=False)

    # global metrics (pDockQ via DockQ package is temporarily skipped)
    ptm = getattr(result, "ptm", None)
    iptm = getattr(result, "iptm", None)
    ptm = float(ptm) if ptm is not None else float("nan")
    iptm = float(iptm) if iptm is not None else float("nan")

    pd.DataFrame(
        [{"metric": "pTM", "value": f"{ptm:.3f}"},
         {"metric": "ipTM", "value": f"{iptm:.3f}"},
         {"metric": "mean pLDDT", "value": f"{mean_plddt:.3f}"},
         {"metric": "n_tokens", "value": n_total}]
    ).to_csv(os.path.join(job_dir, "metrics.csv"), index=False)

    # PAE heatmap + matrix
    pae = getattr(result, "pae", None)
    pae_path = None
    if pae is not None:
        pae = np.asarray(pae, dtype=float)
        try:
            np.save(os.path.join(job_dir, "pae.npy"), pae)
        except Exception:
            pass
        pae_path = plot_pae(job_dir, pae, token_counts_)

    # per-chain pLDDT bar plot
    plot_chain_plddt(job_dir, chain_plddt)

    return {"ptm": ptm, "iptm": iptm, "mean_plddt": mean_plddt,
            "pae": pae_path, "cif": cif_path,
            "metrics_csv": os.path.join(job_dir, "metrics.csv")}


def plot_pae(job_dir, pae, token_counts_):
    """PAE heatmap PNG with chain boundary ticks."""
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(pae, cmap="Greens_r", vmin=0, vmax=30, origin="upper")
    boundaries = np.cumsum(token_counts_)[:-1]
    for b in boundaries:
        ax.axhline(b - 0.5, color="k", linewidth=0.6)
        ax.axvline(b - 0.5, color="k", linewidth=0.6)
    ax.set_title("ESMFold2 PAE")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(im, label="PAE (Å)", shrink=0.85)
    path = os.path.join(job_dir, "pae.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_chain_plddt(job_dir, chain_plddt):
    """Per-chain mean pLDDT bar plot."""
    items = sorted(chain_plddt.items())
    width = max(5.0, 0.35 * len(items))
    fig, ax = plt.subplots(figsize=(width, 4.2), dpi=120)
    ax.bar([c for c, _ in items], [v for _, v in items], color="#7fb3d5")
    ax.set_ylabel("mean pLDDT")
    ax.set_ylim(0, 100)
    for i, (_, v) in enumerate(items):
        ax.text(i, v + 1, f"{v:.0f}", ha="center", fontsize=8)
    ax.set_xticks(range(len(items)))
    ax.set_xticklabels([c for c, _ in items], rotation=45, ha="right")
    ax.set_title("mean pLDDT per chain")
    fig.tight_layout()
    path = os.path.join(job_dir, "chain_plddt.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def refresh_results(output_dir):
    """Recursively scan output_dir for job dirs (each holds a job.json).

    Jobs are stored at <output>/run_<ts>/<job>/job.json (two levels deep), so
    the scan walks the whole tree instead of listing one level only.
    """
    global RESULTS
    RESULTS = []
    if not os.path.isdir(output_dir):
        return []
    for root, _dirs, files in os.walk(output_dir):
        if "job.json" not in files:
            continue
        try:
            with open(os.path.join(root, "job.json"), "r", encoding="utf-8") as f:
                meta = json.load(f)
            meta["dir"] = root
            meta["name"] = os.path.basename(root)
            RESULTS.append(meta)
        except Exception:
            continue
    RESULTS.sort(key=lambda m: m.get("time", ""), reverse=True)
    return RESULTS


def zip_job(job_dir):
    """Create <job_dir>.zip next to the folder; returns the path."""
    zip_path = job_dir + ".zip"
    if os.path.isfile(zip_path):
        os.remove(zip_path)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(job_dir):
            for name in files:
                if name.endswith(".zip"):
                    continue
                full = os.path.join(root, name)
                z.write(full, os.path.relpath(full, os.path.dirname(job_dir)))
    return zip_path


# === 4) 3Dmol.js viewer ===
def ensure_3dmol_js():
    """Download 3Dmol-min.js locally on first run; returns the file path.

    The viewer page references /3dmol-static/3Dmol-min.js so no CDN is needed
    at display time (offline-friendly after the first download).
    """
    global _3DMOL_JS
    if _3DMOL_JS and os.path.isfile(_3DMOL_JS):
        return _3DMOL_JS
    os.makedirs(STATIC_DIR, exist_ok=True)
    path = os.path.join(STATIC_DIR, "3Dmol-min.js")
    if not (os.path.isfile(path) and os.path.getsize(path) > 100_000):
        urls = [
            "https://3dmol.org/build/3Dmol-min.js",
            "https://cdn.jsdelivr.net/npm/3dmol/build/3Dmol-min.js",
            "https://unpkg.com/3dmol/build/3Dmol-min.js",
        ]
        for url in urls:
            try:
                print("Downloading", url, "->", path, flush=True)
                urllib.request.urlretrieve(url, path)
                if os.path.getsize(path) > 100_000:
                    _3DMOL_JS = path
                    return path
            except Exception as exc:
                print(f"[warn] 3Dmol download failed ({url}): {exc}", flush=True)
        print(f"[warn] Place 3Dmol-min.js manually in: {path}", flush=True)
    else:
        _3DMOL_JS = path
    return path


def make_viewer_html(cif_path):
    """Interactive 3D viewer HTML (3Dmol.js), cartoon colored by pLDDT (B-factor)."""
    if not cif_path or not os.path.isfile(cif_path):
        return "<p>No structure file for this job.</p>"
    if not os.path.isfile(os.path.join(STATIC_DIR, "3Dmol-min.js")):
        return (f"<p>3Dmol viewer JS missing — download <em>3Dmol-min.js</em> into "
                f"<code>{STATIC_DIR}</code> (and restart) to enable the 3D view.</p>")
    with open(cif_path) as f:
        cif = f.read()
    js_cif = json.dumps(cif)   # a JSON string is a valid JS string literal
    return (
        '<div id="esm3d" style="width:100%;height:560px;position:relative;"></div>\n'
        '<script src="/3dmol-static/3Dmol-min.js"></script>\n'
        '<script>\n'
        "var cifData = " + js_cif + ";\n"
        "var viewer = null;\n"
        "function init_esm3d() {\n"
        "  if (viewer || !window.$3Dmol) { return; }\n"
        '  var el = document.getElementById("esm3d"); if (!el) { return; }\n'
        '  viewer = $3Dmol.createViewer("esm3d", {backgroundColor: "white"});\n'
        '  var m = viewer.addModel(cifData, "cif");\n'
        "  // ESMFold2 stores pLDDT*100 in the B-factor column; hue: 0(red,low)->240(blue,high)\n"
        "  viewer.setStyle({}, {cartoon: {colorfunc: function(atom) {\n"
        "    var p = (atom.b || 0) / 100.0;\n"
        "    return 'hsl(' + Math.round(Math.min(100, Math.max(0, p)) * 2.4) + ',70%,50%)';\n"
        "  }}});\n"
        "  viewer.addStyle({hetflag: true}, {stick: {radius: 0.3}});\n"
        "  viewer.zoomTo();\n"
        "  viewer.render();\n"
        "}\n"
        "(function poll(){ if (window.$3Dmol) { init_esm3d(); } else { setTimeout(poll, 150); } })();\n"
        "</script>\n"
    )


# === 5) job runner (Gradio generator) ===
def _unpack_run_args(args):
    """Split the flat Gradio input list back into its logical parts."""
    single_n = 5 * MAX_SINGLE
    combo_n = 4 * MAX_COMBO
    p = list(args)
    head = p[:7]                        # mode, model, use_msa, msa, loops, steps, outdir
    i = 7
    single_flat = p[i:i + single_n]; i += single_n
    single_flags = p[i]; i += 1
    type_a = p[i]; i += 1
    a_flat = p[i:i + combo_n]; i += combo_n
    a_flags = p[i]; i += 1
    type_b = p[i]; i += 1
    b_flat = p[i:i + combo_n]; i += combo_n
    b_flags = p[i]
    return (head[0], head[1], head[2], head[3], head[4], head[5], head[6],
            _slots(single_flat, 5, MAX_SINGLE), single_flags,
            type_a, _slots(a_flat, 4, MAX_COMBO), a_flags,
            type_b, _slots(b_flat, 4, MAX_COMBO), b_flags)


def _oom_hint(exc):
    """Human-readable hint when a run fails on the GPU / swap."""
    msg = str(exc)
    if "out of memory" in msg.lower() or "cuda" in msg.lower():
        return (" (VRAM full — close other GPU apps, lower `steps`, or restart "
                "the program to free memory)")
    return ""


def _emit(msg):
    """Return msg while mirroring it to the terminal (UI + console both get it)."""
    print(msg, flush=True)
    return msg


def on_run_click(*args, progress=gr.Progress()):
    """Submit all planned jobs; yield status strings (also printed to stdout)."""
    global STOP_FLAG
    (mode, model_label, use_msa_on, msa_files, num_loops, num_sampling_steps,
     output_dir, single_slots, single_flags, type_a, a_slots, a_flags,
     type_b, b_slots, b_flags) = _unpack_run_args(args)

    STOP_FLAG["stop"] = False
    checkpoint = CHECKPOINTS.get(model_label, CHECKPOINTS[DEFAULT_MODEL])
    is_fast = checkpoint.endswith("Fast")

    if not os.path.isdir(output_dir):
        yield _emit(f"ERROR: output dir does not exist: {output_dir}")
        return
    try:
        rows = collect_rows(mode, single_slots, single_flags,
                            type_a, a_slots, a_flags, type_b, b_slots, b_flags)
        validate_rows(rows, mode)
    except ValueError as exc:
        yield _emit(f"ERROR: {exc}")
        return

    use_msa = bool(use_msa_on) and not is_fast
    if bool(use_msa_on) and is_fast:
        yield _emit("Note: ESMFold2-Fast is single-sequence only — uploaded MSA files "
                    "are ignored. Switch to 'ESMFold2 (MSA-capable)' to use them.")
    msa_map = load_msa_files(msa_files)
    if use_msa and msa_map:
        known = {r["chain"] for r in rows}
        unknown = sorted(msa_map.keys() - known)
        if unknown:
            yield _emit(f"Note: MSA files without a matching chain are ignored: {unknown}")

    jobs = build_jobs(rows, mode)
    total = len(jobs)
    ccd_needed = needs_ccd(rows)

    # Load CCD first (only when ligands/mods present) while RAM is roomy,
    # then the model; both cached process-wide. Visible progress at each step.
    if ccd_needed:
        yield _emit(f"Planned {total} job(s) with {model_label}. Building CCD dictionary (only once) ...")
        get_builder(need_ccd=True)
        yield _emit("CCD dictionary ready.")
    else:
        yield _emit(f"Planned {total} job(s) with {model_label}. No ligands — skipping the "
                    "~400 MB CCD dictionary entirely.")
    yield _emit("Loading model into memory (~13 GB RAM + VRAM). First time takes "
                "1-2 min; meanwhile close other heavy apps if the disk churns ...")
    get_model(checkpoint)
    yield _emit("Model loaded. Folding ...")

    run_id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(output_dir, "run_" + run_id)
    os.makedirs(run_dir, exist_ok=True)

    progress(0.0, desc=f"planned {total} job(s)")
    for i, job in enumerate(jobs, 1):
        if STOP_FLAG["stop"]:
            yield _emit(f"STOPPED after {i - 1}/{total} job(s) — partial results are in {run_dir}.")
            return
        desc = f"folding {i}/{total}: {job['name']} (~{job_tokens(job)} tokens)"
        progress((i - 1) / total, desc=desc)
        _emit(desc)
        job_dir = os.path.join(run_dir, job["name"])
        try:
            result, chain_ids = fold_job(job, checkpoint, num_loops,
                                         num_sampling_steps,
                                         msa_map=msa_map, use_msa=use_msa)
            summary = write_job_artifacts(job_dir, result, chain_ids,
                                          chain_token_counts(job))
            summary.update({"time": datetime.datetime.now().isoformat(),
                            "status": "done", "tokens": sum(chain_token_counts(job))})
            with open(os.path.join(job_dir, "job.json"), "w", encoding="utf-8") as f:
                json.dump(summary, f, ensure_ascii=False)
            yield _emit(f"DONE {i}/{total}: {job['name']} — pTM={summary['ptm']:.3f}, "
                        f"ipTM={summary['iptm']:.3f}, mean pLDDT={summary['mean_plddt']:.1f}")
        except Exception as exc:
            os.makedirs(job_dir, exist_ok=True)
            with open(os.path.join(job_dir, "job.json"), "w", encoding="utf-8") as f:
                json.dump({"status": "failed", "time": datetime.datetime.now().isoformat(),
                           "error": str(exc), "dir": job_dir}, f, ensure_ascii=False)
            yield _emit(f"FAILED {i}/{total}: {job['name']} — {exc}{_oom_hint(exc)}. Skipping.")
        progress(i / total, desc=f"{i}/{total} done")
    yield _emit(f"All {total} job(s) finished. See the Results tab (run dir: {run_dir}).")


def stop_jobs():
    STOP_FLAG["stop"] = True
    return "Stop requested — the current job finishes first, then the batch halts."


def on_preview_click(*args):
    """Batch preview table for the Setup tab."""
    (mode, _, _, _, _, _, _, single_slots, single_flags,
     type_a, a_slots, a_flags, type_b, b_slots, b_flags) = _unpack_run_args(args)
    try:
        rows = collect_rows(mode, single_slots, single_flags,
                            type_a, a_slots, a_flags, type_b, b_slots, b_flags)
        validate_rows(rows, mode)
        jobs = build_jobs(rows, mode)
        return pd.DataFrame(
            [{"job": j["name"],
              "molecules": ", ".join(sorted({r["component"] for r in j["chains"]})),
              "chains": len(j["chains"]), "est_tokens": job_tokens(j)}
             for j in jobs])
    except ValueError as exc:
        return pd.DataFrame({"error": [str(exc)]})
    except Exception as exc:
        return pd.DataFrame({"error": [str(exc)]})


def _toggle_mode(mode):
    """Show only the input block matching the selected mode."""
    return (gr.update(visible=mode == "Single complex"),
            gr.update(visible=mode == "N x N combination"))


def on_add_slot(flags, n_slots):
    """'+ Add molecule': reveal the first hidden row.

    Returns (flags, *row_updates, message) matching the event outputs
    [state, *rows, status].
    """
    flags = list(flags or [False] * n_slots)
    msg = ""
    if all(flags):
        msg = f"max {n_slots} rows reached"
    else:
        i = flags.index(False)
        flags[i] = True
    updates = [gr.update(visible=bool(flags[k])) for k in range(n_slots)]
    return tuple([flags, *updates, msg])


def on_remove_slot(i, flags, n_slots):
    """'✕': hide one row slot."""
    flags = list(flags or [False] * n_slots)
    flags[i] = False
    updates = [gr.update(visible=bool(flags[k])) for k in range(n_slots)]
    return tuple([flags, *updates, ""])


# === 6) Gradio UI ===
def on_refresh(outdir):
    """Rescan the output dir and repopulate the results tab."""
    refresh_results(outdir)
    if not RESULTS:
        return gr.update(choices=[], value=None), RESULTS, ""
    return gr.update(choices=[r["name"] for r in RESULTS],
                     value=RESULTS[0]["name"]), RESULTS, f"{len(RESULTS)} job(s) found."


def on_select(name, results_state):
    """Fill the results tab with the selected job's artifacts."""
    if not results_state:
        return (None, None, None, None, None, "")
    meta = next((r for r in results_state if r["name"] == name), None)
    if meta is None:
        return (None, None, None, None, None, "")
    job_dir = meta.get("dir")
    cif = os.path.join(job_dir, "complex.cif")
    cif = cif if os.path.isfile(cif) else None
    pae_png = os.path.join(job_dir, "pae.png")
    pae_png = pae_png if os.path.isfile(pae_png) else None
    bar_png = os.path.join(job_dir, "chain_plddt.png")
    bar_png = bar_png if os.path.isfile(bar_png) else None
    metrics = None
    mpath = os.path.join(job_dir, "metrics.csv")
    if os.path.isfile(mpath):
        metrics = pd.read_csv(mpath)
    zip_path = zip_job(job_dir)
    status = f"status: {meta.get('status', '?')}"
    if meta.get("status") == "failed":
        status += f" — error: {meta.get('error', 'unknown')}"
    return make_viewer_html(cif), pae_png, bar_png, metrics, zip_path, status


def make_ui():
    with gr.Blocks(title="ESMFold2 Predictor") as demo:
        gr.Markdown(
            "# ESMFold2 structure prediction (protein / DNA / RNA / ligand)\n"
            "**Single complex**: `+ Add molecule` → choose **type** → enter **AA sequence** "
            "(or SMILES for `ligand`). **N×N combination**: List A × List B, each box with its "
            "own type. `copies` ≥ 2 builds homomultimers. "
            "Optional MSA (ESMFold2 model only): `.a3m` files whose **stem equals the chain name**."
        )

        with gr.Tab("Setup"):
            with gr.Row():
                mode = gr.Radio(["Single complex", "N x N combination"],
                                value="Single complex", label="Mode", scale=3)
                model = gr.Radio(list(CHECKPOINTS), value=DEFAULT_MODEL,
                                 label="Model", scale=3)

            # ---------- Single complex: dynamic molecule rows ----------
            with gr.Column(visible=True) as single_block:
                add_s = gr.Button("+ Add molecule", variant="secondary")
                single_rows = []
                single_ctrls = []     # per slot: [type, text, ccd, name, copies]
                single_rm = []
                for k in range(MAX_SINGLE):
                    with gr.Row(visible=(k == 0)) as row:
                        td = gr.Dropdown(VALID_TYPES, value="protein",
                                         label=f"m{k + 1} type", min_width=110, scale=1)
                        tx = gr.Textbox(label="sequence / SMILES",
                                        placeholder="AA or bases (SMILES for ligand)",
                                        scale=5)
                        cc = gr.Textbox(label="ccd / mod", placeholder="SEP@5", scale=1)
                        nm = gr.Textbox(label="name (opt.)", scale=2)
                        cp = gr.Number(value=1, label="copies", precision=0, minimum=1,
                                       scale=1)
                        rm = gr.Button("\u2715", min_width=40, scale=0)
                    single_rows.append(row)
                    single_ctrls.append([td, tx, cc, nm, cp])
                    single_rm.append(rm)
                gr.Markdown(
                    "`ccd / mod`: ligand CCD code (e.g. `SAH`) or polymer modification "
                    "`CCD@POSITION` (1-based, e.g. `SEP@5`; comma-separated). "
                    "`copies` ≥ 2 makes a homomultimer (`<name>_1.._n`)."
                )

            # ---------- N x N combination: two typed boxes ----------
            with gr.Column(visible=False) as combo_block:
                with gr.Row():
                    with gr.Column(scale=1):
                        type_a = gr.Dropdown(VALID_TYPES, value="protein", label="List A type")
                        add_a = gr.Button("+ Add to List A", variant="secondary")
                        rows_a = []
                        ctrls_a = []
                        rm_a = []
                        for k in range(MAX_COMBO):
                            with gr.Row(visible=False) as row:
                                tx = gr.Textbox(label="sequence / SMILES",
                                                placeholder="AA, bases or SMILES",
                                                scale=5)
                                cc = gr.Textbox(label="ccd / mod", scale=1)
                                nm = gr.Textbox(label="name (opt.)", scale=2)
                                cp = gr.Number(value=1, label="copies", precision=0,
                                               minimum=1, scale=1)
                                rm = gr.Button("\u2715", min_width=40, scale=0)
                            rows_a.append(row)
                            ctrls_a.append([tx, cc, nm, cp])
                            rm_a.append(rm)
                    with gr.Column(scale=1):
                        type_b = gr.Dropdown(VALID_TYPES, value="ligand", label="List B type")
                        add_b = gr.Button("+ Add to List B", variant="secondary")
                        rows_b = []
                        ctrls_b = []
                        rm_b = []
                        for k in range(MAX_COMBO):
                            with gr.Row(visible=False) as row:
                                tx = gr.Textbox(label="sequence / SMILES",
                                                placeholder="AA, bases or SMILES",
                                                scale=5)
                                cc = gr.Textbox(label="ccd / mod", scale=1)
                                nm = gr.Textbox(label="name (opt.)", scale=2)
                                cp = gr.Number(value=1, label="copies", precision=0,
                                               minimum=1, scale=1)
                                rm = gr.Button("\u2715", min_width=40, scale=0)
                            rows_b.append(row)
                            ctrls_b.append([tx, cc, nm, cp])
                            rm_b.append(rm)
                gr.Markdown(
                    "Each row is one component; every component of List A is crossed with "
                    "every component of List B (A₁×B₁, A₁×B₂, …). `ccd / mod` and `copies` "
                    "as in single mode."
                )

            mode.change(_toggle_mode, mode, [single_block, combo_block])

            # slot visibility state
            single_flags = gr.State([True] + [False] * (MAX_SINGLE - 1))
            a_flags = gr.State([False] * MAX_COMBO)
            b_flags = gr.State([False] * MAX_COMBO)

            # ---------- MSA / folding / output controls ----------
            with gr.Row():
                use_msa = gr.Checkbox(value=False, label="Use MSA",
                                      info="ESMFold2 only; ignored by ESMFold2-Fast")
                msa = gr.File(file_count="multiple", file_types=[".a3m", ".aln", ".fasta"],
                              label="MSA files", height=80, scale=4)
                with gr.Column(scale=2):
                    with gr.Row():
                        loops = gr.Number(value=20, label="loops", info="refinement loops")
                        steps = gr.Number(value=100, label="steps", info="diffusion steps")
                    outdir = gr.Textbox(value=os.path.join(BASE_DIR, "esmfold2_runs"),
                                        label="Output directory")

            with gr.Row():
                run_btn = gr.Button("Run", variant="primary")
                stop_btn = gr.Button("Stop")
                preview_btn = gr.Button("Preview batch list")
            status = gr.Markdown("Ready.")
            preview_tbl = gr.Dataframe(label="Batch preview", interactive=False)

            # add / remove wiring (states must exist first)
            add_s.click(lambda f: on_add_slot(f, MAX_SINGLE), [single_flags],
                        [single_flags, *single_rows, status])
            add_a.click(lambda f: on_add_slot(f, MAX_COMBO), [a_flags],
                        [a_flags, *rows_a, status])
            add_b.click(lambda f: on_add_slot(f, MAX_COMBO), [b_flags],
                        [b_flags, *rows_b, status])
            for k, rm in enumerate(single_rm):
                rm.click(lambda i=k: on_remove_slot(i, single_flags, MAX_SINGLE),
                         [single_flags], [single_flags, *single_rows, status])
            for k, rm in enumerate(rm_a):
                rm.click(lambda i=k: on_remove_slot(i, a_flags, MAX_COMBO),
                         [a_flags], [a_flags, *rows_a, status])
            for k, rm in enumerate(rm_b):
                rm.click(lambda i=k: on_remove_slot(i, b_flags, MAX_COMBO),
                         [b_flags], [b_flags, *rows_b, status])

            # run / preview inputs: slot values flattened in slot order
            single_vals = [c for slot in single_ctrls for c in slot]
            a_vals = [c for slot in ctrls_a for c in slot]
            b_vals = [c for slot in ctrls_b for c in slot]

            run_event = run_btn.click(on_run_click,
                          [mode, model, use_msa, msa, loops, steps, outdir,
                           *single_vals, single_flags,
                           type_a, *a_vals, a_flags,
                           type_b, *b_vals, b_flags],
                          status)
            stop_btn.click(stop_jobs, None, status)
            preview_btn.click(on_preview_click,
                              [mode, model, use_msa, msa, loops, steps, outdir,
                               *single_vals, single_flags,
                               type_a, *a_vals, a_flags,
                               type_b, *b_vals, b_flags],
                              preview_tbl)

        with gr.Tab("Results"):
            with gr.Row():
                refresh_btn = gr.Button("Refresh jobs")
                job_sel = gr.Dropdown(label="Job", choices=[])
            results_state = gr.State([])
            job_status = gr.Markdown("")
            viewer = gr.HTML(label="3D structure (pLDDT-colored)")
            metrics = gr.Dataframe(label="Metrics")
            with gr.Row():
                pae_plot = gr.Image(label="PAE heatmap")
                plddt_plot = gr.Image(label="mean pLDDT per chain")
            job_zip = gr.DownloadButton("Download job zip", variant="secondary")

            refresh_btn.click(on_refresh, outdir, [job_sel, results_state, job_status])
            job_sel.change(on_select, [job_sel, results_state],
                           [viewer, pae_plot, plddt_plot, metrics, job_zip, job_status])

        with gr.Tab("Help"):
            gr.Markdown(
                "## Model\n"
                "- **ESMFold2-Fast (default)** — inference-optimized single-sequence "
                "model; several times faster than ESMFold2. MSA files are ignored.\n"
                "- **ESMFold2 (MSA-capable)** — use the `Use MSA` switch and upload "
                "`.a3m` per protein/RNA chain (file stem = chain name).\n"
                "## Single complex\n"
                "`+ Add molecule` opens a new row; every row is one molecule of the SAME "
                "complex.\n"
                "- **type** — `protein` / `dna` / `rna` / `ligand`; the text field then holds "
                "AA / bases / SMILES respectively.\n"
                "- **ccd / mod** — ligand CCD code (e.g. `SAH`, `HEM`) or polymer "
                "modification `CCD@POSITION` (e.g. `SEP@5`; comma-separated).\n"
                "- **name (opt.)** — chain label; blank = auto (`chain1`, …). MSA files must "
                "use this name as the file stem.\n"
                "- **copies** — homomultimer size (≥ 1); chain gets a `_1.._n` suffix.\n"
                "## N×N combination\n"
                "Two boxes, each with its own molecule type. Every row (= one component) of "
                "List A is crossed with every row of List B (A₁×B₁, A₁×B₂, …); each job = "
                "one A component + one B component.\n"
                "## Metrics\n"
                "- Native: per-residue pLDDT, PAE matrix, pTM, ipTM, per-chain mean pLDDT.\n"
                "- pDockQ (reference-free, DockQ package) is temporarily skipped.\n"
                "## GPU / memory\n"
                "Requires a CUDA GPU. The 6B model needs ~13 GB RAM + ~13 GB VRAM — close "
                "other heavy apps (browsers, games) before folding; if loading seems to "
                "freeze it is usually disk swap, give it a few minutes the first time. "
                "Total tokens are hard-capped at 2048; reducing `loops` / `steps` speeds "
                "folding up. Model + CCD are loaded once per program run and reused."
            )

        demo.queue()

        # auto-refresh the Results job list after a run finishes
        run_event.then(on_refresh, outdir, [job_sel, results_state, job_status])
    return demo


def main():
    os.makedirs(os.path.join(BASE_DIR, "esmfold2_runs"), exist_ok=True)
    demo = make_ui()
    # serve the local 3Dmol.js so the viewer works without a CDN.
    # NOTE: gradio rebuilds its FastAPI app on launch(), so mount AFTER launch.
    ensure_3dmol_js()
    demo.launch(inbrowser=True, show_error=True, prevent_thread_lock=True)
    if os.path.isdir(STATIC_DIR):
        from fastapi.staticfiles import StaticFiles
        try:
            # note: gradio already owns "/static"; use a distinct prefix
            demo.app.mount("/3dmol-static", StaticFiles(directory=STATIC_DIR),
                           name="3dmol-static")
            print(f"Static viewer files served from {STATIC_DIR}", flush=True)
        except Exception as exc:
            print(f"[warn] could not mount /3dmol-static: {exc}", flush=True)
    try:
        demo.block_thread()
    except Exception:
        import time
        while True:
            time.sleep(3600)


if __name__ == "__main__":
    main()