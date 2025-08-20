#! /bin/bash
#Author: Hao He, 2025
#System: Ubuntu-22.04
#Env: gromacs_2024


conda activate gromacs_2024

cd /media/hao/Data/Project_MYC/2025-07-18_AURKA_MYC_MD

input_protein=AURKA_NMYC_clean.pdb
protein_name=AURKA_NMYC # can't be same with input_protein

input_ligand=OTSSP167_diffdock.sdf
# ligand_name=OTS # 3-letter name
# input_ligand="" # if "" will omit ligand


# ---------------- MD preparation ------------------#
## 0. re-docking
# re-dock protein-protein/ligand if nesscerry

## 1. primary clean
# remove unwanted chain or ions, water, or reducntant ligands

## 2. get ligand
# select only ligand
# save selected ligand to sdf
# current export sdf from pymol is fine, the others are not
# add hyrogen to ligand, change ligand to ligand 3-LETTER name
obabel -isdf $input_ligand --title $ligand_name -h -osdf -O ${ligand_name}.sdf

## 3. get protein
# select only protein (multiple chians), saved as pdb
# if have non-standard res, using pdbfixer to fix it first
# ref to: ~/Documents/GitHub/myScript/python_functions/md_helper.py -> prepare_protein()
# proteonatize protein (ref to streadmd mannual)
pdb2pqr $input_protein ${protein_name}.pdb --ffout AMBER --keep-chain


# ---------------- MD stimulation ------------------#
# forcefield: ~/miniconda3/envs/gromacs_2024/share/gromacs/top
if [ "$input_ligand" == "" ]; then
	md_dir=${protein_name}
else
	md_dir=${protein_name}_${ligand_name}
fi

run_md \
	--wdir $md_dir \
	--nvt_time 1000 --npt_time 1000 --md_time 200 \
	--device "gpu" \
	--protein ${protein_name}.pdb \
	--ligand ${ligand_name}.sdf


## post-processing trajectory

md_res_dir=${md_dir}/md_files/md_run/${md_dir}
# md_res_dir=/media/hao/Data/Project_MYC/2025-07-18_AURKA_MYC_MD/5G1X/md_files/md_run/prot

if [ "$input_ligand" == "" ]; then
	# check index
	# echo -e "q" | gmx make_ndx -n $md_res_dir/index.ndx -o $md_res_dir/index.ndx
	keep_group=1
else
	# make new index only for [protein_ligand] if not created
	# echo -e "1 | 13\nname 20 protein_ligand\nq" | \
	# 	gmx make_ndx -f $md_res_dir/md_out.tpr -n $md_res_dir/index.ndx -o $md_res_dir/index.ndx
	keep_group=20
fi
# remove pbc and water_ions
echo -e "$keep_group" | gmx trjconv -f $md_res_dir/md_out.xtc -s $md_res_dir/md_out.tpr \
	-n $md_res_dir/index.ndx -o $md_res_dir/md_out_nojump.xtc -pbc nojump
echo -e "4\n$keep_group" | gmx trjconv -f $md_res_dir/md_out_nojump.xtc -s $md_res_dir/md_out.tpr \
	-n $md_res_dir/index.ndx -o $md_res_dir/md_out_nojump_fit.xtc -fit rot+trans
echo -e "1\n$keep_group" | gmx trjconv -f $md_res_dir/md_out_nojump_fit.xtc -s $md_res_dir/md_out.tpr \
	-n $md_res_dir/index.ndx -o $md_res_dir/md_out_nojump_fit_center.xtc -pbc mol -center
echo -e "$keep_group" | gmx trjconv -f $md_res_dir/md_out_nojump_fit_center.xtc -s $md_res_dir/md_out.tpr \
	-n $md_res_dir/index.ndx -o $md_res_dir/frame_clean.pdb -dump 0



# remove redundant files
rm ${ligand_name}.sdf, ${protein_name}.pdb, ${protein_name}.log







# ---------------- MD analysis ------------------#
## ProLIF analysis
run_prolif --wdir_to_run md_files/md_run/myc_aurka_OTS \
	--wdir md_files/md_run/myc_aurka_OTS/prolif \
	--ligand 'UNL' \
	--protein_selection 'protein' \
	--step 1


## MM-PBSA/GBSA analysis
conda activate gmx_mmpbsa
run_gbsa --wdir_to_run md_files/md_run/protein_h_OTS \
	--wdir md_files/md_run/protein_h_OTS/mm_pbsa \
	--ligand_id 'UNL'


# # add new index for protein-protein analysis
# cd ~/Downloads/results_MYC_1x143__AURKA_122x403__OTS_f4d50/streamd/md_files/md_run/myc_aurka_OTS
# gmx make_ndx -f md_out.tpr -n index.ndx -o index_protein_chains.ndx

index_protein_chains




# ----------------  myc_aurka ------------------#

## MD stimulation
run_md --wdir . \
	--protein myc_aurka.pdb \
	--ncpu 20 --device auto --ntmpi_per_gpu 2 \
	--nvt_time 500 --npt_time 1000 --md_time 200 \
	--save_traj_without_water
	# --out_suffix myc_aurka \

## ProLIF analysis
run_prolif --wdir_to_run md_files/md_run/myc_aurka \
	--wdir md_files/md_run/myc_aurka/prolif \
	--protein_selection 'protein' \
	--step 1


## MM-PBSA/GBSA analysis
conda activate gmx_mmpbsa
run_gbsa --wdir_to_run md_files/md_run/myc_aurka \
	--wdir md_files/md_run/myc_aurka/mm_pbsa


# add new index for protein-protein analysis




