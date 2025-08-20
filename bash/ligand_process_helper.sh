#!/bin/bash


# split_sdf_into_chunks '.' 'mcule_ultimate_express1_220828.sdf' 20000 'chunks'
split_sdf_into_chunks() {
    local work_dir=$1
    local input_sdf=$2
    local chunk_size=$3
    local output_dir=$4

    local basename_unit="$((chunk_size / 1000))k"
    local basename="chunk_${basename_unit}"

    cd "$work_dir" || exit
    mkdir -p "$output_dir"

    # Count number of molecules
    local total_mols
    total_mols=$(obabel "$input_sdf" -omol -O /dev/null 2>&1 | grep -c "1 molecule converted")
    if [ "$total_mols" -eq 0 ]; then
        total_mols=$(grep -c "^\\$\\$\\$\\$" "$input_sdf")
    fi
    echo "Total molecules: $total_mols"

    local start=1
    local chunk_index=1
    while [ "$start" -le "$total_mols" ]; do
        local end=$((start + chunk_size - 1))
        if [ "$end" -gt "$total_mols" ]; then
            end=$total_mols
        fi
        local output_file="${basename}_${chunk_index}.sdf"
        echo "Creating $output_file with molecules $start to $end..."
        obabel "$input_sdf" -f "$start" -l "$end" -O "${output_dir}/$output_file"
        start=$((end + 1))
        chunk_index=$((chunk_index + 1))
    done

    echo "Done splitting into chunks of $chunk_size."
}




# sdf_to_3d_minimized '.' 'chunks' 'chunks_3d' 24
sdf_to_3d_minimized() {
    local working_dir=$1
    local input_dir=$2
    local output_dir=$3
    local thread=$4

    cd "$working_dir" || exit
    mkdir -p "$output_dir"

    # Export output_dir to be available in the subprocess
    export output_dir

    parallel -j $thread ' 
        input_file={}
        base_name=$(basename "${input_file%.*}")
        output_file="${output_dir}/${base_name}_3d.sdf"
        echo "Processing $input_file -> $output_file"
        obabel "$input_file" -O "$output_file" \
            --gen3d med \
            -h \
            --minimize --ff MMFF94 --steps 1000 \
            --errorlevel 1
    ' ::: ${input_dir}/*.sdf

    echo "Processed into 3D SDFs with hydrogens and minimized."
}





# Function to convert SDF files to PDBQT using unidocktools

# # Example usage of the function:
# working_dir="/home/hao/docking_tools/cmpd_library/VS_Campaign"
# library_sdf="diverse_20k.sdf"   # Can be a file or directory
# save_pdbqt_dir="diverse_20k_pdbqt"
# use_file_name="true"  # Set to "true" to use --use_file_name in the command

# # Call the function
# convert_to_pdbqt $working_dir $library_sdf $save_pdbqt_dir 24 'true'

convert_to_pdbqt() {
    local working_dir="$1"        # Directory containing the SDF file or folder
    local library_sdf="$2"        # SDF file or directory of SDF files
    local save_pdbqt_dir="$3"     # Directory to save converted PDBQT files
    local thread="$4"
    local use_file_name="$5"      # Boolean flag to use file name in the output (true/false)

    cd "$working_dir" || exit
    mkdir -p "$save_pdbqt_dir"  # Create the save directory if it doesn't exist


    # Determine whether to include --use_file_name in the command
    local use_file_name_flag=""
    if [ "$use_file_name" == "true" ]; then
        use_file_name_flag="--use_file_name"
    fi

    # Check if library_sdf is a file or directory
    if [ -f "$library_sdf" ]; then
        # If it's a file, just convert it to PDBQT
        echo "Converting single file: $library_sdf to PDBQT"
        unidocktools ligandprep \
            --ligands "$library_sdf" \
            --savedir "$save_pdbqt_dir" \
            --save_format pdbqt \
            $use_file_name_flag

    elif [ -d "$library_sdf" ]; then
        # If it's a directory, process all SDF files in the directory in parallel
        echo "Converting all SDF files in directory: $library_sdf to PDBQT"

        # Export output_dir to be available in the subprocess
        export save_pdbqt_dir

        # Find all SDF files in the directory and process them in parallel
        find "$library_sdf" -type f -name "*.sdf" | parallel -j $thread '
            input_file={}
            base_name=$(basename "${input_file%.*}")
            output_file="${save_pdbqt_dir}/${base_name}.pdbqt"
            echo "Converting $input_file -> $output_file"
            unidocktools ligandprep \
                --ligands "$input_file" \
                --savedir "$save_pdbqt_dir" \
                --save_format pdbqt \
                '$use_file_name_flag'
        '
    else
        echo "Error: $library_sdf is neither a valid file nor a directory."
        exit 1
    fi

    echo "Conversion complete."
}


