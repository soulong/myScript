#! /bin/bash

# by Hao He, 2025-04-09



conda activate unidock

source /home/hao/Documents/GitHub/myScript/bash_function/ligand_process_helper.sh




# downloaded from https://mcule.com/database/
# subset: ultimate_express1
# Aug. 28, 2022	
# total 486,387
# format 2D sdf
cd '/home/hao/docking_tools/cmpd_library/mcule_ultimate_express1_220828'
# split into small chunks
split_sdf_into_chunks '.' 'mcule_ultimate_express1_220828.sdf' 20000 'chunks'
# generate 3d and minimize
sdf_to_3d_minimized '.' 'chunks' 'sdf_3d' 25
# sdf to pdbqt
convert_to_pdbqt '.' 'sdf_3d' 'pdbqt' 25 'true'
# convert_to_pdbqt '.' 'mcule_ultimate_express1_220828.sdf' 'pdbqt' 25 'true'



cd '/home/hao/docking_tools/cmpd_library/VS_Campaign'
# split into small chunks
split_sdf_into_chunks '.' 'diverse_20k.sdf' 2000 'chunks'
# generate 3d and minimize
sdf_to_3d_minimized '.' 'chunks' 'sdf_3d' 10
# sdf to pdbqt
convert_to_pdbqt '.' 'sdf_3d' 'pdbqt' 10 'true'




