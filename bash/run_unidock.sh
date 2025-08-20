#! /bin/bash
# env: unidock


# Setup -------------------
working_dir="/home/hao/docking_tools/VS_HPK1"
cd $working_dir

recepor_file="7M0M.pdbqt"
center_x=-2.93
center_y=-6.681
center_z=-9.034
size_x=15.0
size_y=15.0
size_z=15.0

# ligand_dir="HPK1_inhitors_pdbqt"
# output_dir="vina_result_HPK1_inhitors"

# ligand_dir="/home/hao/docking_tools/cmpd_library/VS_Campaign/pdbqt"
# output_dir="vina_result_diverse_20k"

ligand_dir="/home/hao/docking_tools/cmpd_library/mcule_ultimate_express1_220828/pdbqt"
output_dir="vina_result_mcule_ue1"


batch_size=1000 # only works with batch

cd $working_dir
mkdir -p $output_dir


# ## process active cmpds -------------------
# # run only once
# source /home/hao/Documents/GitHub/myScript/bash_function/ligand_process_helper.sh
# sdf_to_3d_minimized '.' 'HPK1_inhitors_sdf' 'HPK1_inhitors_sdf_3d' 24
# convert_to_pdbqt '.' 'HPK1_inhitors_sdf_3d' 'HPK1_inhitors_pdbqt' 24 'false'


## Option 1: with progress bar -------------------
progress_bar() {
    local progress=$1
    local total=$2
    local width=50  # width of progress bar
    local percent=$(( 100 * progress / total ))
    local filled=$(( width * progress / total ))
    local empty=$(( width - filled ))
    local bar=$(printf "%${filled}s" | tr ' ' '#')
    bar+=$(printf "%${empty}s" | tr ' ' '-')
    echo -ne "\r[${bar}] ${percent}% ($progress/$total)"
}

# Get a list of all ligand files
ligand_files=(${ligand_dir}/*.pdbqt)
total=${#ligand_files[@]}
processed=0


## Option 1.1: with batch
echo "Starting Uni-Dock docking with $total ligands (batch size = $batch_size)..."
while [ $processed -lt $total ]; do
    batch=("${ligand_files[@]:$processed:$batch_size}")

    # Construct the batch ligand argument string
    gpu_batch_args=()
    for lig in "${batch[@]}"; do
        gpu_batch_args+=("$lig")
    done

    # Run Uni-Dock on this batch
    unidock \
        --receptor "$recepor_file" \
        --gpu_batch "${gpu_batch_args[@]}" \
        --center_x "$center_x" \
        --center_y "$center_y" \
        --center_z "$center_z" \
        --size_x "$size_x" \
        --size_y "$size_y" \
        --size_z "$size_z" \
        --search_mode balance \
        --scoring vina \
        --num_modes 1 \
        --dir "$output_dir" \
        --verbosity 0 \
        --seed 42 \
        >> "docking.log" 2>&1

    processed=$((processed + ${#batch[@]}))
    progress_bar "$processed" "$total"
done
echo -e "\n✅ Docking complete!"




## Option 1.2: without batch
echo "Starting Uni-Dock docking with $total ligands..."
for lig in "${ligand_files[@]}"
do
	unidock \
		--receptor "$recepor_file" \
		--gpu_batch "$lig" \
		--search_mode balance \
		--scoring vina \
		--center_x "$center_x" \
		--center_y "$center_y" \
		--center_z "$center_z" \
		--size_x "$size_x" \
		--size_y "$size_y" \
		--size_z "$size_z" \
		--num_modes 1 \
		--dir "$output_dir" \
		--verbosity 0 \
		>> "docking.log" 2>&1

	((processed++))
	progress_bar "$processed" "$total"
done
echo -e "\n✅ Docking complete!"





## Option 2: no progress bar -------------------
for lig in ${ligand_dir}/*.pdbqt
do
	echo $lig
	unidock \
		--receptor $recepor_file \
		--gpu_batch $lig \
		--search_mode balance \
		--scoring vina \
        --center_x $center_x \
        --center_y $center_y \
        --center_z $center_z \
        --size_x $size_x \
        --size_y $size_y \
        --size_z $size_z \
		--num_modes 1 \
		--dir $output_dir \
		--verbosity 0 \
		>> "docking.log" 2>&1
done
