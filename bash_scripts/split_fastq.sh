#!/bin/bash

# conda activate ngs
# Split FASTQ data into multiple non-overlapping replicates using BBTools reformat.sh
# Supports both single-end and paired-end sequencing data
# Usage: 
#   Single-end: ./split_fastq_bbtools.sh sample.fq.gz --split_ratio 0.5 0.5
#   Paired-end: ./split_fastq_bbtools.sh sample_1.fq.gz sample_2.fq.gz --split_ratio 0.33 0.33 0.34

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <fastq1> [fastq2] --split_ratio <ratio1> <ratio2> [ratio3] ... [--seed <seed>]"
    echo ""
    echo "Arguments:"
    echo "  fastq1              First FASTQ file (or single-end file)"
    echo "  fastq2              Second FASTQ file (for paired-end, optional)"
    echo "  --split_ratio       Split ratios for each replicate (must sum to 1.0)"
    echo "  --seed              Random seed for reproducibility (default: 42)"
    echo ""
    echo "Examples:"
    echo "  Single-end:  $0 sample.fq.gz --split_ratio 0.5 0.5"
    echo "  Paired-end:  $0 sample_1.fq.gz sample_2.fq.gz --split_ratio 0.33 0.33 0.34"
    echo "  Three reps:  $0 sample.fq.gz --split_ratio 0.4 0.3 0.3 --seed 123"
    exit 1
}

# Function to check if a command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is not installed or not in PATH"
        echo "Please install BBTools and ensure reformat.sh is accessible"
        exit 1
    fi
}

# Function to validate split ratios
validate_split_ratios() {
    local ratios=("$@")
    local sum=0
    local ratio
    
    # Check if any ratios provided
    if [ ${#ratios[@]} -eq 0 ]; then
        echo "Error: No split ratios provided"
        return 1
    fi
    
    # Check each ratio and calculate sum
    for ratio in "${ratios[@]}"; do
        # Check if ratio is a valid number
        if ! [[ "$ratio" =~ ^[0-9]*\.?[0-9]+$ ]]; then
            echo "Error: Invalid ratio '$ratio'. Must be a positive number."
            return 1
        fi
        
        # Check if ratio is positive
        if (( $(echo "$ratio <= 0" | bc -l) )); then
            echo "Error: All split ratios must be positive. Got: $ratio"
            return 1
        fi
        
        sum=$(echo "$sum + $ratio" | bc -l)
    done
    
    # Check if sum equals 1.0 (with small tolerance for floating point errors)
    if (( $(echo "$sum < 0.999999 || $sum > 1.000001" | bc -l) )); then
        echo "Error: Split ratios must sum to 1.0, got $sum"
        return 1
    fi
    
    return 0
}

# Function to get base name from FASTQ filename
get_base_name() {
    local filename="$1"
    local basename
    
    # Get just the filename without path
    basename=$(basename "$filename")
    
    # Remove .gz if present
    basename="${basename%.gz}"
    
    # Remove common FASTQ extensions
    basename="${basename%.fastq}"
    basename="${basename%.fq}"
    
    echo "$basename"
}

# Function to generate output filename
get_output_filename() {
    local input_file="$1"
    local replicate_num="$2"
    local pair_num="$3"  # Optional, for paired-end
    
    local base_name
    base_name=$(get_base_name "$input_file")
    
    # Handle paired-end naming
    if [ -n "$pair_num" ]; then
        # Remove existing pair indicators
        base_name=$(echo "$base_name" | sed -E 's/_[12]$|_R[12]$|\.R[12]$//')
        output_name="${base_name}_${pair_num}_rep${replicate_num}"
    else
        output_name="${base_name}_rep${replicate_num}"
    fi
    
    # Add appropriate extension
    if [[ "$input_file" == *.gz ]]; then
        output_name="${output_name}.fq.gz"
    else
        output_name="${output_name}.fq"
    fi
    
    echo "$output_name"
}

# Function to detect if files are paired-end
detect_paired_end() {
    local file1="$1"
    local file2="$2"
    
    if [ -z "$file2" ]; then
        return 1  # Not paired-end
    fi
    
    local base1 base2
    base1=$(get_base_name "$file1")
    base2=$(get_base_name "$file2")
    
    # Remove pair indicators and check if base names match
    local clean_base1 clean_base2
    clean_base1=$(echo "$base1" | sed -E 's/_[12]$|_R[12]$|\.R[12]$//')
    clean_base2=$(echo "$base2" | sed -E 's/_[12]$|_R[12]$|\.R[12]$//')
    
    if [ "$clean_base1" = "$clean_base2" ]; then
        return 0  # Paired-end detected
    else
        return 1  # Not paired-end
    fi
}

# Function to split single-end FASTQ
split_single_end() {
    local input_file="$1"
    local seed="$2"
    shift 2
    local ratios=("$@")
    
    local num_replicates=${#ratios[@]}
    local cumulative_ratio=0
    
    echo "Splitting single-end file: $input_file"
    echo "Split ratios: ${ratios[*]}"
    echo "Random seed: $seed"
    echo ""
    
    # Process each replicate
    for i in $(seq 0 $((num_replicates - 1))); do
        local ratio=${ratios[$i]}
        local replicate_num=$((i + 1))
        local output_file
        output_file=$(get_output_filename "$input_file" "$replicate_num" "")
        
        echo "Creating replicate $replicate_num: $output_file (${ratio} = $(echo "$ratio * 100" | bc -l | xargs printf "%.1f")%)"
        
        if [ $i -eq $((num_replicates - 1)) ]; then
            # Last replicate gets all remaining reads
            local skip_ratio
            skip_ratio=$(echo "$cumulative_ratio" | bc -l)
            
            reformat.sh \
                in="$input_file" \
                out="$output_file" \
                skipreads="$skip_ratio" \
                -Xmx4g \
                overwrite=true
        else
            # Calculate the sample ratio for this replicate
            local sample_ratio
            sample_ratio=$(echo "scale=10; $ratio / (1 - $cumulative_ratio)" | bc -l)
            
            local skip_ratio
            skip_ratio=$(echo "$cumulative_ratio" | bc -l)
            
            reformat.sh \
                in="$input_file" \
                out="$output_file" \
                skipreads="$skip_ratio" \
                samplerate="$sample_ratio" \
                seed="$seed" \
                -Xmx4g \
                overwrite=true
        fi
        
        cumulative_ratio=$(echo "$cumulative_ratio + $ratio" | bc -l)
    done
}

# Function to split paired-end FASTQ
split_paired_end() {
    local input_file1="$1"
    local input_file2="$2"
    local seed="$3"
    shift 3
    local ratios=("$@")
    
    local num_replicates=${#ratios[@]}
    local cumulative_ratio=0
    
    echo "Splitting paired-end files: $input_file1, $input_file2"
    echo "Split ratios: ${ratios[*]}"
    echo "Random seed: $seed"
    echo ""
    
    # Process each replicate
    for i in $(seq 0 $((num_replicates - 1))); do
        local ratio=${ratios[$i]}
        local replicate_num=$((i + 1))
        local output_file1 output_file2
        output_file1=$(get_output_filename "$input_file1" "$replicate_num" "1")
        output_file2=$(get_output_filename "$input_file2" "$replicate_num" "2")
        
        echo "Creating replicate $replicate_num: $output_file1, $output_file2 (${ratio} = $(echo "$ratio * 100" | bc -l | xargs printf "%.1f")%)"
        
        if [ $i -eq $((num_replicates - 1)) ]; then
            # Last replicate gets all remaining reads
            local skip_ratio
            skip_ratio=$(echo "$cumulative_ratio" | bc -l)
            
            reformat.sh \
                in1="$input_file1" \
                in2="$input_file2" \
                out1="$output_file1" \
                out2="$output_file2" \
                skipreads="$skip_ratio" \
                -Xmx4g \
                overwrite=true
        else
            # Calculate the sample ratio for this replicate
            local sample_ratio
            sample_ratio=$(echo "scale=10; $ratio / (1 - $cumulative_ratio)" | bc -l)
            
            local skip_ratio
            skip_ratio=$(echo "$cumulative_ratio" | bc -l)
            
            reformat.sh \
                in1="$input_file1" \
                in2="$input_file2" \
                out1="$output_file1" \
                out2="$output_file2" \
                skipreads="$skip_ratio" \
                samplerate="$sample_ratio" \
                seed="$seed" \
                -Xmx4g \
                overwrite=true
        fi
        
        cumulative_ratio=$(echo "$cumulative_ratio + $ratio" | bc -l)
    done
}

# Main function
main() {
    # Check if BBTools is available
    check_command "reformat.sh"
    check_command "bc"
    
    # Parse arguments
    local fastq1="" fastq2="" split_ratios=() seed=42
    local parse_ratios=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            --split_ratio)
                parse_ratios=true
                shift
                ;;
            --seed)
                seed="$2"
                shift 2
                ;;
            -h|--help)
                usage
                ;;
            *)
                if [ "$parse_ratios" = true ]; then
                    # Check if this looks like a ratio (number with optional decimal)
                    if [[ "$1" =~ ^[0-9]*\.?[0-9]+$ ]]; then
                        split_ratios+=("$1")
                    else
                        parse_ratios=false
                        # This might be fastq2, handle it below
                    fi
                fi
                
                if [ "$parse_ratios" = false ]; then
                    if [ -z "$fastq1" ]; then
                        fastq1="$1"
                    elif [ -z "$fastq2" ]; then
                        fastq2="$1"
                    else
                        echo "Error: Too many input files specified"
                        usage
                    fi
                fi
                shift
                ;;
        esac
    done
    
    # Validate inputs
    if [ -z "$fastq1" ]; then
        echo "Error: No input files specified"
        usage
    fi
    
    if [ ! -f "$fastq1" ]; then
        echo "Error: Input file '$fastq1' not found"
        exit 1
    fi
    
    if [ -n "$fastq2" ] && [ ! -f "$fastq2" ]; then
        echo "Error: Input file '$fastq2' not found"
        exit 1
    fi
    
    # Set default split ratios if none provided
    if [ ${#split_ratios[@]} -eq 0 ]; then
        split_ratios=(0.5 0.5)
        echo "No split ratios specified, using default: ${split_ratios[*]}"
    fi
    
    # Validate split ratios
    if ! validate_split_ratios "${split_ratios[@]}"; then
        exit 1
    fi
    
    local num_replicates=${#split_ratios[@]}
    echo "Splitting into $num_replicates replicates"
    
    # Determine if single-end or paired-end and process accordingly
    if [ -z "$fastq2" ]; then
        echo "Detected single-end sequencing data"
        split_single_end "$fastq1" "$seed" "${split_ratios[@]}"
    else
        if detect_paired_end "$fastq1" "$fastq2"; then
            echo "Detected paired-end sequencing data"
        else
            echo "Warning: File names don't follow typical paired-end naming convention"
            echo "Proceeding as paired-end data..."
        fi
        split_paired_end "$fastq1" "$fastq2" "$seed" "${split_ratios[@]}"
    fi
    
    echo ""
    echo "Splitting completed successfully!"
    echo "Output files created in the current directory."
}

# Run main function with all arguments
main "$@"