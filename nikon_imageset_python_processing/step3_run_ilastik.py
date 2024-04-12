# -*- coding: utf-8 -*-

import argparse
import os, glob, random, re, time, math
from datetime import timedelta
# import itertools
from natsort import natsorted
import subprocess
from functools import partial
from tqdm.contrib.concurrent import process_map
# import itertools
# import pandas as pd



# run ilastik prediction on images
def ilastik_prediction(images,
                       ilastik_proj,  # load saved checkpoint
                       ilastik_exec="D:/Software/Ilastik/ilastik-1.4.0b21-gpu/ilastik.exe"  # operetta exported measurement directory
                       ):
    # device = torch.device("cuda:" + str(gpu_id) if torch.cuda.is_available() else 'cpu')
    # print('using device:', device)
    # if torch.cuda.is_available():
    #     print('device name:', torch.cuda.get_device_name(gpu_id))
    #     torch.cuda.empty_cache()

    # cmd = f'"{ilastik_exec}" --headless --project="{ilastik_proj}" {images}'
    # cmd = f'nohup "{ilastik_exec}" --headless --project="{ilastik_proj}" {images} &'
    # cmd = f'LAZYFLOW_THREADS=2 LAZYFLOW_TOTAL_RAM_MB=2000 "{ilastik_exec}" --readonly --headless --project="{ilastik_proj}" {images}'
    cmd = f'"{ilastik_exec}" --headless --project="{ilastik_proj}" {images}'
    subprocess.run(cmd, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)




if __name__ == "__main__":
    # get parameters
    parser = argparse.ArgumentParser(description='Predict images exported from operetta using segmentation models')
    parser.add_argument('--dir', type=str, nargs="+",
                        help='directories (parent or subfolder) exported from operetta')
    parser.add_argument('--image_suffix', type=str, default=".tiff", metavar=".tiff", help='image suffix')
    parser.add_argument('--subset_pattern', type=str, default = None, metavar = 'None',
                        help = 'use regrex to subset images, eg. "r0[1]c0[1-2]f.+p.+-" used to subset A1 and A2 wells')
    parser.add_argument('--output_suffix', type=str, default="Probabilities", metavar="Probabilities",
                        help='output suffix of ilastik, please check it before')
    parser.add_argument('--ilastik_proj', type=str, help='trained ilastik project file path')
    parser.add_argument('--ilastik_exec', type=str,
                        default="D:/Software/Ilastik/ilastik-1.4.0b21-gpu/ilastik.exe",
                        metavar="D:/Software/Ilastik/ilastik-1.4.0b21-gpu/ilastik.exe",
                        help='ilastik executable program path')
    parser.add_argument('--batchsize', type=int, default=50, metavar=50,
                        help='split images to batches on single core, if 0, means auto adjust to [thread_num]')
    parser.add_argument('--thread_num', type=int, default=8, metavar=8, help='number of threads to use')
    parser.add_argument('--overwrite_mask', action="store_true", help='overwrite existed masks')
    # parser.add_argument('--gpu_id', type=int, default=0, metavar=0, help='GPU id to use, if avaiable')
    parser.add_argument('--test', action='store_true', help='test on random indicated number of images')
    parser.add_argument('--test_num', type=int, default=4, metavar=4,
                        help='test number of images, only work if --test exists')
    args = parser.parse_args()
    
    
    # measurements = ["D:/Postdoc_Data/2024-01-06_condesate_reporter_LJY"]
    # subdir = "images"
    # # subdir = None
    # subset_pattern = ".*_c2.tiff"
    # # subset_pattern = None
    # ilastik_proj="D:/Postdoc_Data/2023-12-27_condesate_reporter_LJY/MyProject.ilp"
    
    
    measurements = args.dir
    image_suffix = args.image_suffix
    subset_pattern =args.subset_pattern
    output_suffix = args.output_suffix
    ilastik_proj = args.ilastik_proj
    ilastik_exec = args.ilastik_exec
    batchsize = args.batchsize
    thread_num = args.thread_num
    overwrite_mask = args.overwrite_mask
    test = args.test
    test_num = args.test_num
    
    if True:
        measurements=["D:\\Postdoc_Data\\2024-02-19_src_project"]
        image_suffix = ".tiff"
        subset_pattern = ".*_c1.tiff"
        output_suffix = "Probabilities"
        ilastik_proj = "D:\\Postdoc_Data\\2024-02-19_src_project\\MyProject.ilp"
        ilastik_exec = "D:/Software/Ilastik/ilastik-1.4.0b21-gpu/ilastik.exe"
        batchsize = 50
        thread_num = 6
        overwrite_mask = False
        test = False
        test_num = 4

    
    # check parameters
    if measurements is None:
        quit('no valid measurement found, check it!')
    if not os.path.exists(ilastik_proj):
        quit(f"ilastik proj file not found, please check:\n  {ilastik_proj}")
    if not os.path.exists(ilastik_exec):
        quit(f"ilastik excute file not found, please check:\n  {ilastik_exec}")

    # process datasets one by one
    for measurement in measurements:
        print(">>>>> processing: " + measurement)
        
        os.chdir(measurement)
        subdir = "images"
        if subdir is not None:
            os.chdir(subdir)
        
        # get images
        files = glob.glob("*" + image_suffix)
       
        # remove any files containing output suffix
        files = [f for f in files if not re.search(output_suffix + ".tiff", f)]
        print("found %d intensity files" % (len(files)))
        
        # subset_pattern = "r0[1]c0[1-2]f.+p.+-"
        # subset_pattern = None
        if subset_pattern is not None:
            files = [f for f in files if re.search(subset_pattern, f)]
            print("found %d files, filtered by %s" % (len(files), subset_pattern))

        # remove founded segmentated images
        if overwrite_mask:
            print(f"overwrite any existed masks, will process {len(files)} images")
        else:
            files_save_name = [f.replace(".tiff",  "_" + output_suffix + ".tiff") for f in files]
            # print(files_save_name[0:5])
            segmentated_files = []
            for idx, f in enumerate(files_save_name):
                if os.path.isfile(f):
                    segmentated_files.append(files[idx])
            print("ignore %d predicted masks from total %d images contain %s" %
                  (len(segmentated_files), len(files), subset_pattern))
            files = list(set(files) - set(segmentated_files))
        files = natsorted(files)
        # print(files[0:5])
       
        # test mode
        if test is True:
            print("enable test mode, process random 4 files")
            files = random.sample(files, 4)
        
        
        # run in parallel
        partial_func = partial(ilastik_prediction,
                               ilastik_exec=ilastik_exec,
                               ilastik_proj=ilastik_proj)
        # batch in each cpu core
        batch_images = []
        batch_size = math.ceil(len(files) / 6)
        print(f'batch size: {batch_size}')
        for i in range(0, len(files), batch_size):
            batch = files[i:i+batch_size]
            # ' '.join([str(elem) for elem in batch])
            batch_images.append(' '.join(f'"{w}"' for w in batch))
        # print(batch_images)
        
        # run in parallel
        start_time = time.time()
        _ = process_map(partial_func, batch_images, max_workers=thread_num)
        
        # # run all in a batch
        # start_time = time.time()
        # ilastik_prediction(" ".join(files), ilastik_proj)
        
        print("elapsed time: " + str(timedelta(seconds=int(time.time() - start_time))))
       
    print("\n----------------- done -----------------")
    
    
    


#     # find operreta_datasets
#     # dir = ['../AR_operreta/2023-03-04_CellCarrier-96_AR_Cmpds__2023-03-04T17_47_34-Measurement 1',
#     #        '../AR_operreta/2023-03-08_CellCarrier-96_AR-V7_Cmpds__2023-03-08T17_59_35-Measurement 1']
#     measurements = "D:\\Others\\2023-11-08_LJY_SPARK_spot_intensity\\images"

#     # process datasets one by one
#     for measurement in measurements:
#         print(">>>>> processing: " + measurement)
#         # get images
#         files = glob.glob(os.path.join(measurement, "Images", "*" + args.image_suffix))
#         # subset only operetta images
#         files = [f for f in files if re.search('r[0-9]{2}c[0-9]{2}f[0-9]{2,}p[0-9]{2,}-ch[0-9]{1,}sk[0-9]{1,}fk1fl1' + f'({args.image_suffix})', f)]
#         # print(files[0:5])

#         # subset_pattern = "r0[1]c0[1-2]f.+p.+-"
#         if args.subset_pattern is not None:
#             files = [f for f in files if re.search(args.subset_pattern, f)]
#             print("found %d files, filtered by %s" % (len(files), args.subset_pattern))
#             # print(len(files))

#         # remove founded segmentated images
#         if args.overwrite_masks:
#             print(f"overwrite any existed masks, will process {len(files)} images")
#         else:
#             files_save_name = [f.replace(args.image_suffix, "_" + args.output_suffix + args.image_suffix) for f in files]
#             # print(files_save_name[0:5])
#             segmentated_files = []
#             for idx, f in enumerate(files_save_name):
#                 if os.path.isfile(f):
#                     segmentated_files.append(files[idx])
#             print("ignore %d predicted masks from total %d images contain %s" %
#                   (len(segmentated_files), len(files), args.subset_pattern))
#             files = list(set(files) - set(segmentated_files))
#         files = natsorted(files)
#         # print(files[0:5])

#         # test mode
#         if test is True:
#             print(f"enable test mode, process random {args.test_num} files")
#             files = random.sample(files, args.test_num)

#         # using partial
#         partial_func = partial(ilastik_prediction,
#                                ilastik_exec=args.ilastik_exec,
#                                ilastik_proj=args.ilastik_proj)
#         # batch in each cpu core
#         batch_images = []
#         batch_size = math.ceil(len(files) / 6)
#         print(f'batch size: {batch_size}')
#         for i in range(0, len(files), batch_size):
#             batch = files[i:i+batch_size]
#             # ' '.join([str(elem) for elem in batch])
#             batch_images.append(' '.join(f'"{w}"' for w in batch))
#         # print(batch_images)

#         # run in parallel
#         start_time = time.time()
#         _ = process_map(partial_func, batch_images, max_workers=a)

#         print("elapsed time: " + str(timedelta(seconds=int(time.time() - start_time))))

#     print("\n----------------- done -----------------")
