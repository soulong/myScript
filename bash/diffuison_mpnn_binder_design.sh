#! /bin/bash

BASEDIR="/media/hao/Data/Project_kinase_substrate_design/2025-05-12_ALK_sub"
cd ~/RFdiffusion


## step 1
## diffusion backbone
conda activate SE3nv
diffusion_output_dir=$BASEDIR/diffusion

scripts/run_inference.py \
	inference.input_pdb=${BASEDIR}/ALK-KD_ATP_Mg.pdb \
	inference.output_prefix=$diffusion_output_dir/sub \
	'contigmap.contigs=[A1-277/0 11]' \
	'ppi.hotspot_res=[A134,A138,A174,A176,A178]' \
	inference.num_designs=400 \
	denoiser.noise_scale_ca=0.2 \
	denoiser.noise_scale_frame=0.2

## add constrain (optional)
conda activate proteinmpnn_binder_design
constrained_pdb_dir=$BASEDIR/diffusion_constrained

# modify pdb to add constrain residue annd pdb label
python '/home/hao/Documents/GitHub/myScript/python_functions/constrain_residues.py' \
	--constraints "A:M0:Y" \
	--output_dir $constrained_pdb_dir \
	$diffusion_output_dir



## step 2
## protein-mpnn sequence
conda activate proteinmpnn_binder_design
cd ~/dl_binder_design

protein_mpnn_dir=$BASEDIR/protein_mpnn
mpnn_fr/dl_interface_design.py \
	-pdbdir $constrained_pdb_dir \
	-outpdbdir $protein_mpnn_dir \
	-checkpoint_path mpnn_fr/ProteinMPNN/soluble_model_weights/v_48_020.pt \
	-relax_cycles 0 \
	-seqs_per_struct 2 \
	-temperature 0.0001 \
	-debug



## step 3
## AF2 prediction
conda activate af2_binder_design
cd ~/dl_binder_design
af2_dir=$BASEDIR/af2

af2_initial_guess/predict.py \
	-pdbdir $protein_mpnn_dir \
	-outpdbdir $af2_dir \
	-scorefilename $BASEDIR/score.sc \
	-debug
## organize
python ~/Documents/GitHub/myScript/python_functions/tidy_mpnn_result.py \
	$af2_dir $BASEDIR/score.sc $BASEDIR/af2_score.csv



## step 3 (optional method, AlphaPulldown prediction)
# pdb to fasta
alphapulldown_dir=$BASEDIR/alphapulldown
mkdir -p $alphapulldown_dir

python ~/Documents/GitHub/myScript/python_functions/extract_pdb_sequence.py \
	$BASEDIR/ALK-KD_ATP_Mg.pdb -o $BASEDIR/candidate.fasta --chain A
python ~/Documents/GitHub/myScript/python_functions/extract_pdb_sequence.py \
	$BASEDIR/protein_mpnn -o $BASEDIR/bait.fasta --chain A
# extract id to txt
python ~/Documents/GitHub/myScript/python_functions/extract_fasta_ids.py \
	$BASEDIR/candidate.fasta
python ~/Documents/GitHub/myScript/python_functions/extract_fasta_ids.py \
	$BASEDIR/bait.fasta

conda activate alphapulldown
create_individual_features.py \
  --fasta_paths="$BASEDIR/candidate.fasta,$BASEDIR/bait.fasta" \
  --data_dir="/media/hao/Data1/alphafold_data_reduced/" \
  --output_dir=$BASEDIR/msa \
  --save_msa_files=False \
  --skip_existing=True \
  --use_mmseqs2=True \
  --max_template_date=2050-01-01 \
  --db_preset reduced_dbs
 run_multimer_jobs.py \
  --mode=pulldown \
  --num_cycle=3 \
  --num_predictions_per_model=1 \
  --model_names=model_5_multimer_v3 \
  --output_path=$alphapulldown_dir \
  --data_dir="/media/hao/Data1/alphafold_data_reduced/" \
  --protein_lists="$BASEDIR/candidate.txt,$BASEDIR/bait.txt" \
  --monomer_objects_dir=$BASEDIR/msa \
  --compress_result_pickles=True
   # --models_to_relax=Best \
singularity exec \
    --no-home \
    --bind $alphapulldown_dir/:/mnt \
    ~/alpha-analysis_jax_0.4.sif \
    run_get_good_pae.sh \
    --output_dir=/mnt \
    --cutoff=200


