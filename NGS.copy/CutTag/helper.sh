#!/usr/bin/env bash
# helper.sh - Common functions for genomic pipelines
# Usage: source this file in your main pipeline script

set -euo pipefail

# ============================= Logging Functions =============================

log() {
    # Print timestamped log message to stderr
    # Usage: log "message"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

log_info() {
    # Print info message
    # Usage: log_info "message"
    log "[INFO] $*"
}

log_warn() {
    # Print warning message
    # Usage: log_warn "message"
    log "[WARN] $*"
}

log_error() {
    # Print error message
    # Usage: log_error "message"
    log "[ERROR] $*"
}

die() {
    # Print error and exit
    # Usage: die "error message"
    log_error "$*"
    exit 1
}

# ============================= File Validation =============================

require_file() {
    # Check if file exists and is not empty
    # Usage: require_file "/path/to/file" ["custom error message"]
    local file="$1"
    local msg="${2:-Missing or empty required file: $file}"
    
    if [[ ! -f "$file" ]]; then
        die "$msg (file not found)"
    elif [[ ! -s "$file" ]]; then
        die "$msg (file is empty)"
    fi
}

require_dir() {
    # Check if directory exists
    # Usage: require_dir "/path/to/dir" ["custom error message"]
    local dir="$1"
    local msg="${2:-Missing required directory: $dir}"
    
    [[ -d "$dir" ]] || die "$msg"
}

require_command() {
    # Check if command is available
    # Usage: require_command "samtools" ["custom error message"]
    local cmd="$1"
    local msg="${2:-Required command not found: $cmd}"
    
    if ! command -v "$cmd" &>/dev/null; then
        die "$msg"
    fi
}

# ============================= File Processing =============================

safe_read_csv() {
    # Safely read CSV file, handling DOS line endings and headers
    # Usage: safe_read_csv <csv_file> <callback_function>
    # Callback receives: sample group control target fq1 fq2
    local file="$1"
    local callback="$2"
    
    require_file "$file"
    
    while IFS=, read -r sample group control target fq1 fq2 || [[ -n "${sample:-}" ]]; do
        # Skip header
        [[ "$sample" == "sample" ]] && continue
        # Skip empty lines
        [[ -z "$sample" ]] && continue
        
        "$callback" "$sample" "$group" "$control" "$target" "$fq1" "$fq2"
    done < <(tail -n +2 "$file" | tr -d '\r')
}

get_base_name() {
    # Extract base name from FASTQ file, removing common extensions
    # Usage: get_base_name "sample_R1.fastq.gz"
    # Returns: sample_R1
    local file="$1"
    local name
    name=$(basename "$file")
    name=${name%.fastq.gz}
    name=${name%.fq.gz}
    name=${name%.fastq}
    name=${name%.fq}
    echo "$name"
}

# ============================= Statistics Extraction =============================

get_flagstat_value() {
    # Extract value from samtools flagstat output
    # Usage: get_flagstat_value <flagstat_file> <pattern> [column]
    # Example: get_flagstat_value "sample.flagstat" "mapped" 1
    local file="$1"
    local pattern="$2"
    local col="${3:-1}"
    
    if [[ ! -f "$file" ]]; then
        log_warn "Flagstat file not found: $file"
        echo "0"
        return 1
    fi
    
    awk -v p="$pattern" -v c="$col" '
        $0 ~ p { print $c; exit }
    ' "$file"
}

# ============================= FASTQ Generation =============================

generate_samplesheet() {
    # Auto-generate samplesheet by scanning for paired-end FASTQ files
    # Usage: generate_samplesheet [input_dir] [output_csv]
    # Assumes R1/R2 naming convention
    local input_dir="${1:-.}"
    local output_csv="${2:-samplesheet.csv}"
    
    if [[ ! -d "$input_dir" ]]; then
        log_error "Directory does not exist: $input_dir"
        return 1
    fi
    
    log_info "Generating samplesheet from: $input_dir"
    
    # Write header
    printf "sample,group,control,target,fq1,fq2\n" > "$output_csv"
    
    # Find R1 files and infer R2
    find "$input_dir" -type f -regextype posix-extended \
        -regex '.*_[rR]?1(_.*)?\.f(ast)?q\.gz' | sort | while IFS= read -r fq1; do
        
        local fq1_dir fq1_base fq2_base fq2
        fq1_dir=$(dirname "$fq1")
        fq1_base=$(basename "$fq1")
        
        # Convert R1 to R2
        fq2_base=$(echo "$fq1_base" | sed -E 's/^(.*)(_[rR]?)1((_.*)?\.f(ast)?q\.gz)$/\1\22\3/')
        fq2="$fq1_dir/$fq2_base"
        
        if [[ ! -f "$fq2" ]]; then
            log_warn "Missing mate pair: $fq2 (derived from $fq1)"
            fq2="NA"
        fi
        
        printf ",,,,%s,%s\n" "$fq1" "$fq2" >> "$output_csv"
    done
    
    log_info "Generated samplesheet: $output_csv"
    log_info "Please fill in sample, group, control, and target columns manually"
}

# ============================= Peak Processing =============================

merge_peaks_to_consensus() {
    # Merge multiple peak files into consensus BED
    # Usage: merge_peaks_to_consensus -i "peak1.bed peak2.bed" -o output.bed -g genome.sizes
    local input_peaks=""
    local output_bed=""
    local genome_file=""
    
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -i|--input)
                input_peaks="$2"
                shift 2
                ;;
            -o|--output)
                output_bed="$2"
                shift 2
                ;;
            -g|--genome)
                genome_file="$2"
                shift 2
                ;;
            -h|--help)
                echo "Usage: merge_peaks_to_consensus -i 'file1 file2' -o output.bed -g genome.sizes"
                echo "  -i, --input    Space-separated list of peak files"
                echo "  -o, --output   Output consensus BED file"
                echo "  -g, --genome   Genome sizes file for sorting"
                return 0
                ;;
            *)
                log_error "Unknown option: $1"
                return 1
                ;;
        esac
    done
    
    # Validate inputs
    [[ -z "$input_peaks" ]] && { log_error "Input peaks required (-i)"; return 1; }
    [[ -z "$output_bed" ]] && { log_error "Output file required (-o)"; return 1; }
    [[ -z "$genome_file" ]] && { log_error "Genome file required (-g)"; return 1; }
    require_file "$genome_file"
    
    # Convert space-separated string to array
    local -a peak_files
    read -ra peak_files <<< "$input_peaks"
    
    if [[ ${#peak_files[@]} -eq 0 ]]; then
        log_error "No peak files provided"
        return 1
    fi
    
    # Validate all peak files exist
    local missing=0
    for pf in "${peak_files[@]}"; do
        if [[ ! -f "$pf" ]]; then
            log_warn "Peak file not found: $pf"
            ((missing++))
        fi
    done
    
    if [[ $missing -gt 0 ]]; then
        log_error "Some peak files are missing"
        return 1
    fi
    
    log_info "Merging ${#peak_files[@]} peak files into $output_bed"
    
    # Merge peaks: concatenate -> extract coords -> sort -> merge -> sort by genome
    cat "${peak_files[@]}" | \
        cut -f1-3 | \
        bedtools sort 2>/dev/null | \
        bedtools merge 2>/dev/null | \
        bedtools sort -g "$genome_file" 2>/dev/null > "$output_bed"
    
    if [[ -s "$output_bed" ]]; then
        local count
        count=$(wc -l < "$output_bed")
        log_info "Consensus peaks created: $output_bed ($count regions)"
        return 0
    else
        log_error "Failed to generate consensus peaks"
        rm -f "$output_bed"
        return 1
    fi
}

# ============================= Normalization =============================

calculate_scale_factor() {
    # Calculate normalization scale factor
    # Usage: calculate_scale_factor <method> <filtered_reads> [spike_reads] [sample] [spike_free_file]
    # Returns: scale factor as float
    local method="$1"
    local filtered_reads="$2"
    local spike_reads="${3:-0}"
    local sample="${4:-unknown}"
    local spike_free_file="${5:-}"
    
    local scale_factor="1.0"
    
    case "$method" in
        CPM)
            if [[ $filtered_reads -gt 0 ]]; then
                scale_factor=$(awk "BEGIN{printf \"%.6f\", 1e6 / $filtered_reads}")
            else
                log_warn "Filtered reads = 0, using scale factor 1.0"
            fi
            ;;
        Spike)
            if [[ $spike_reads -gt 0 ]]; then
                scale_factor=$(awk "BEGIN{printf \"%.6f\", 1e4 / $spike_reads}")
            else
                log_warn "Spike reads = 0, using scale factor 1.0"
            fi
            ;;
        SpikeFree)
            if [[ -z "$spike_free_file" || ! -f "$spike_free_file" ]]; then
                die "SpikeFree method requires spike-free scale factor file"
            fi
            
            local sf_raw
            sf_raw=$(awk -F'\t' -v s="$sample" '
                NR>1 {
                    gsub(/\.filtered\.bam$/, "", $1)
                    if ($1 == s) { print $7; exit }
                }
            ' "$spike_free_file")
            
            if [[ -z "$sf_raw" ]]; then
                die "Sample $sample not found in $spike_free_file"
            fi
            
            scale_factor=$(awk -v fr="$filtered_reads" -v sf="$sf_raw" \
                'BEGIN{printf "%.6f", 1e6 / (fr * sf)}')
            ;;
        *)
            log_error "Unknown normalization method: $method"
            return 1
            ;;
    esac
    
    echo "$scale_factor"
}

# ============================= Index Creation =============================

ensure_fasta_index() {
    # Create FASTA index if it doesn't exist
    # Usage: ensure_fasta_index <fasta_file>
    local fasta="$1"
    local fai="${fasta}.fai"
    
    require_file "$fasta"
    
    if [[ ! -f "$fai" ]]; then
        log_info "Creating FASTA index: $fai"
        samtools faidx "$fasta" || die "Failed to create FASTA index"
    fi
    
    echo "$fai"
}

# ============================= Validation =============================

validate_species() {
    # Validate species code
    # Usage: validate_species "hs" || die "Invalid species"
    local species="$1"
    
    case "$species" in
        hs|chm13|mm)
            return 0
            ;;
        *)
            log_error "Invalid species: $species (must be: hs, chm13, mm)"
            return 1
            ;;
    esac
}

validate_peak_type() {
    # Validate peak type
    # Usage: validate_peak_type "narrow" || die "Invalid peak type"
    local peak_type="$1"
    
    case "$peak_type" in
        narrow|broad)
            return 0
            ;;
        *)
            log_error "Invalid peak type: $peak_type (must be: narrow, broad)"
            return 1
            ;;
    esac
}

# ============================= Array Utilities =============================

array_contains() {
    # Check if array contains element
    # Usage: array_contains "element" "${array[@]}" && echo "found"
    local seeking="$1"
    shift
    local element
    
    for element; do
        [[ "$element" == "$seeking" ]] && return 0
    done
    return 1
}

get_unique_values() {
    # Get unique values from array
    # Usage: unique_vals=($(get_unique_values "${all_vals[@]}"))
    printf '%s\n' "$@" | sort -u
}

# ============================= Export Functions =============================

# Make functions available to subshells if needed
export -f log log_info log_warn log_error die
export -f require_file require_dir require_command
export -f get_base_name get_flagstat_value
export -f validate_species validate_peak_type