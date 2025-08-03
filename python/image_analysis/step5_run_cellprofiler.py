# -*- coding: utf-8 -*-

import argparse
import math
import os
import time
# import itertools
from shutil import copyfile
import in_place
# import random
from datetime import timedelta
import sqlite3
from sqlalchemy import create_engine
from pathlib import Path
import pandas as pd
# import tqdm
# import numpy as np
from natsort import natsorted
# import warnings
# warnings.filterwarnings('ignore', '.*grid_sample and affine_grid behavior.*', )
# warnings.filterwarnings('ignore', '.*Named tensors and all their associated.*', )

from subprocess import run, Popen, PIPE  # , DEVNULL
from multiprocessing import Process


# run cellprofiler on dataloaders
def cellprofiler_single(cp_dataloader_path,
                        cp_project_path,  # load saved checkpoint
                        cp_exec_path  # operetta exported measurement directory
                        ):
    output_dir = os.path.dirname(cp_dataloader_path)
    cmd = f'. ~/miniconda3/etc/profile.d/conda.sh && conda activate cellprofiler && conda info && cellprofiler -c -r -p "{cp_project_path}" -o "{output_dir}" --data-file "{cp_dataloader_path}"'
    # cmd = f'{cp_exec_path} -c -r -p "{cp_project_path}" -o "{output_dir}" --data-file "{cp_dataloader_path}"'
    logfile = open(os.path.join(output_dir, 'cellprofiler.log'), 'w')
    run(cmd, shell=True, stdout=logfile, stderr=logfile)
    # process = Popen(cmd, shell=True, stdout=logfile, stderr=logfile, text=True)
    # out, err = process.communicate()
    return None


def cellprofiler_multiple(
    measurement: str = None,
    cp_project_path: str = None, # full path
    cp_exec_path: str = "/home/hao/miniconda3/envs/cellprofiler/bin/cellprofiler",
    thread_num: int = 16
):
    # get dataloader info
    dataloader_filename = "cp_dataloader.csv"
    dataloader_path = os.path.join(measurement, dataloader_filename)
    if not os.path.exists(dataloader_path):
        print(f"dataloader file {dataloader_path} not found, skip")
        return None
    else:
        dataloader = pd.read_csv(dataloader_path, index_col=False)
        # print(dataloader)
    
    if not os.path.exists(cp_project_path):
        print(f"project file {cp_project_path} not found, skip")
        return None

    cp_result_merged_name = "cp_result.db"
    if os.path.exists(os.path.join(measurement, cp_result_merged_name)):
        # os.remove(cp_result_merged_name)
        print(f"{cp_result_merged_name} found, skip")
        return None

    # split dataloader to multiple processes
    group_size = math.ceil(len(dataloader) / thread_num)
    dataloader_list = [dataloader.iloc[i:i + group_size]
                       for i in range(0, len(dataloader), group_size)]
    dataloader_list_name = [os.path.join(measurement, 'jobs', f'job_{i}', f'job_{i}.csv')
                            for i in range(1, (len(dataloader_list)+1))]
    # save split dataloader
    for idx, item in enumerate(dataloader_list):
        os.makedirs(os.path.dirname(dataloader_list_name[idx]), exist_ok=True)
        item.to_csv(dataloader_list_name[idx], index=False)

    # run cellprofiler
    start_time = time.time()
    print("run cellprofiler subprocess jobs")
    jobs = []
    for dataloader in dataloader_list_name:
        process = Process(
            target=cellprofiler_single,
            args=(dataloader, cp_project_path, cp_exec_path,))
        jobs.append(process)

    [j.start() for j in jobs]
    [j.join() for j in jobs]
    [j.close() for j in jobs]
    print("elapsed time: " + str(timedelta(seconds=int(time.time() - start_time))))

    return None


def merge_db(
        measurement,
        job_dir_name="jobs",
        merged_name="cp_result.db"
):
    
    if os.path.exists(f"{measurement}/{merged_name}"):
        print(f"{merged_name} existed, skip")
        return None
    
    start_time = time.time()
    
    # open conn
    db_conn = create_engine(f'sqlite:///{measurement}/{merged_name}')
    # db_conn = sqlite3.connect(f'{measurement}/{merged_name}')

    # get all split dbs
    job_paths = Path(measurement).joinpath(job_dir_name).rglob('*.db')
    # job_paths = [str(f) for f in job_paths]
    job_paths = natsorted(job_paths)
    N = 0
    for f in job_paths:
        # print(f'processing: {str(f)}')
        con_tmp = create_engine(f'sqlite:///{str(f)}')

        # write extra tables once
        if N == 0:
            # read first db to get tables
            experiment = pd.read_sql('Experiment', con=con_tmp)
            experiment_properties = pd.read_sql('Experiment_Properties', con=con_tmp)
            per_experiment = pd.read_sql('Per_Experiment', con=con_tmp)
            # write to merged db
            experiment.to_sql('Experiment', db_conn, index=False, if_exists='replace')
            experiment_properties.to_sql('Experiment_Properties', db_conn,
                                         index=False, if_exists='replace')
            per_experiment.to_sql('Per_Experiment', db_conn, index=False, if_exists='replace')

            # get object tables
            tmp = sqlite3.connect(f'{str(f)}')
            cursor = tmp.cursor()
            cursor.execute("""SELECT name FROM sqlite_master WHERE type='table';""")
            object_names = cursor.fetchall()
            object_names = [x[0] for x in object_names]
            object_names = [x for x in object_names if x not in [
                "Per_Image", "Experiment", "sqlite_sequence", "Experiment_Properties", "Per_Experiment"]]
            tmp.close()

        # write per_image and per_object tables
        per_image = pd.read_sql('Per_Image', con=con_tmp)
        # modify ImageNumber
        per_image['ImageNumber'] += N
        # write to db
        per_image.to_sql('Per_Image', db_conn, index=False, if_exists='append')

        for name in object_names:
            per_object = pd.read_sql(name, con=con_tmp)
            per_object['ImageNumber'] += N
            per_object.to_sql(name, db_conn, index=False, if_exists='append')
        # counter
        N += len(per_image)

    # modify Experiment_Properties db_sqlite_file value
    experiment_properties = pd.read_sql('Experiment_Properties', con=db_conn)
    experiment_properties.loc[1, 'value'] = f'{measurement}/{merged_name}'
    # experiment_properties['value'][1] = f'{measurement}/{merged_name}'
    experiment_properties.to_sql('Experiment_Properties', db_conn, index=False, if_exists='replace')

    # get a new .propoerty file
    for name in object_names:
        name = name.replace("Per_", "")
        job_1_file = f"{measurement}/jobs/job_1/cp_result_{name}.properties"
        target_file = f"{measurement}/cp_result_{name}.properties"
        if os.path.exists(job_1_file):
            copyfile(job_1_file, target_file)
            with in_place.InPlace(target_file) as file:
                for line in file:
                    line = line.replace("jobs/job_1/", "")
                    file.write(line)
                file.close()
    
    print("elapsed time: " + str(timedelta(seconds=int(time.time() - start_time))))

    return None



if __name__ == "__main__":
    # get parameters
    parser = argparse.ArgumentParser(
        description='run cellprofiler pipeline with save .cpproj file and .csv dataloader')
    parser.add_argument('--measurements', type=str, nargs="+",
                        help='directories (parent or subfolder) exported from operetta')
    parser.add_argument('--cp_project_path', type=str, help='cellprofiler pipeline file path')
    parser.add_argument('--cp_exec_path', type=str,
                        default='/home/hao/miniconda3/envs/cellprofiler/bin/cellprofiler',
                        metavar='/home/hao/miniconda3/envs/cellprofiler/bin/cellprofiler',
                        help='cellprofiler executable program path')
    parser.add_argument('--thread_num', type=int, default=16,
                        metavar=16, help='number of threads to use')
    args = parser.parse_args()

    if not os.path.exists(args.cp_project_path):
        quit(f"cellprofiler proj file {args.cellprofiler_proj} not found, please check")
    if not os.path.exists(args.cp_exec_path):
        quit(f"cellprofiler execute file {args.cellprofiler_exec} not found, please check")

    # dir = ['../AR_operreta/2023-03-04_CellCarrier-96_AR_Cmpds__2023-03-04T17_47_34-Measurement 1',
    #             '../AR_operreta/2023-03-08_CellCarrier-96_AR-V7_Cmpds__2023-03-08T17_59_35-Measurement 1']
    # measurement_dirs = [[d[0] for d in os.walk(m)] for m in args.dir]
    # measurement_dirs = list(itertools.chain(*measurement_dirs))
    # measurement_dirs = [os.path.abspath(d) for d in measurement_dirs if re.search("-Measurement \\d+$", d)]
    measurements = natsorted(args.measurements)
    # print(measurement_dirs)

    # # get all sub-jobs
    # jobs = [j[0] for m in measurement_dirs for j in os.walk(m) if re.search("job_", j[0])]

    # # remove no image jobs
    # jobs_Images = [os.path.join(os.path.dirname(os.path.dirname(j)), "Images") for j in jobs]
    # jobs = [j for idx, j in enumerate(jobs) if os.path.exists(jobs_Images[idx])]

    # # search for killed sub-jobs
    # jobs_killed = [j for j in jobs if not os.path.exists(os.path.join(j, "cp_result.properties"))]

    # # get related directory
    # jobs_killed_dataloader = [os.path.join(j, os.path.basename(j) + ".csv") for j in jobs_killed]

    # # remove sub-jobs result
    # for j in jobs_killed:
    #     job_result = os.path.join(j, "cp_result.db")
    #     if os.path.exists(job_result): os.remove(job_result)

    # # using partial run
    # partial_func = partial(run_cp_pipeline,
    #                        cellprofiler_exec=args.cellprofiler_exec,
    #                        cellprofiler_proj=args.cellprofiler_proj)

    # _ = process_map(partial_func, jobs_killed_dataloader, max_workers=args.thread_num)

    # process measurement one by one
    for measurement in measurements:
        print(">>>>> processing: " + measurement)

        # run cellprofiler
        cellprofiler_multiple(
            measurement,
            args.cp_project_path,
            args.cp_exec_path,
            args.thread_num)

        # merge multiple db
        print("merge cellprofiler subprocess results")
        merge_db(measurement, job_dir_name="jobs", merged_name="cp_result.db")

    print("\n----------------- done -----------------")
