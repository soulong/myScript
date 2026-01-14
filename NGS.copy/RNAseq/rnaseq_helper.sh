#!/usr/bin/env bash
# rnaseq_helper.sh - Helper functions for RNA-seq pipeline
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

log_duration() {
    # Log elapsed time
    # Usage: log_duration $SECONDS "Task name"
    local seconds="$1"
    local task_name="${2:-Task}"
    local minutes=$((seconds / 60))
    local remaining=$((seconds % 60))
    log_info "[$task_name] Completed in ${minutes}m ${remaining}s"
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

check_dependencies() {
    # Check all required tools are installed
    # Usage: check_dependencies
    log_info "Checking dependencies..."
    
    local missing_tools=()
    local required_tools=(
        "fastqc" "fastp" "salmon" "STAR" "samtools" 
        "featureCounts" "multiqc"
    )
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &>/dev/null; then
            missing_tools+=("$tool")
            log_warn "  ✗ $tool (not found)"
        else
            log_info "  ✓ $tool"
        fi
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        log_error "Please install missing tools before running the pipeline"
        return 1
    fi
    
    log_info "All dependencies satisfied"
    return 0
}

# ============================= FASTQ Processing =============================

detect_fastq_suffix() {
    # Detect R1/R2 suffixes from FASTQ files
    # Usage: detect_fastq_suffix <fastq_dir>
    # Returns: "suffix_r1 suffix_r2"
    local fastq_dir="$1"
    
    require_dir "$fastq_dir"
    
    # Find all gzipped FASTQ files
    local -a fastq_files
    mapfile -t fastq_files < <(find "$fastq_dir" -type f -name "*.gz" | sort -u)
    
    if [[ ${#fastq_files[@]} -eq 0 ]]; then
        log_error "No gzipped FASTQ files found in: $fastq_dir"
        return 1
    fi
    
    # Detect R1 and R2 suffixes
    local suffix_r1 suffix_r2
    suffix_r1=$(printf '%s\n' "${fastq_files[@]}" | \
        grep -oE "_([rR])?1\.f(ast)?q\.gz" | sort -u | head -n1)
    suffix_r2=$(printf '%s\n' "${fastq_files[@]}" | \
        grep -oE "_([rR])?2\.f(ast)?q\.gz" | sort -u | head -n1)
    
    if [[ -z "$suffix_r1" || -z "$suffix_r2" ]]; then
        log_error "Could not detect R1/R2 suffixes in FASTQ files"
        return 1
    fi
    
    echo "$suffix_r1 $suffix_r2"
}

collect_paired_fastq() {
    # Collect paired FASTQ files and extract sample names
    # Usage: collect_paired_fastq <fastq_dir> <suffix_r1> <suffix_r2>
    # Populates global arrays: FASTQ_R1, FASTQ_R2, SAMPLE_NAMES
    local fastq_dir="$1"
    local suffix_r1="$2"
    local suffix_r2="$3"
    
    declare -ga FASTQ_R1=() FASTQ_R2=() SAMPLE_NAMES=()
    
    # Find all R1 files
    local -a r1_files
    mapfile -t r1_files < <(find "$fastq_dir" -type f -name "*${suffix_r1}" | sort)
    
    if [[ ${#r1_files[@]} -eq 0 ]]; then
        log_error "No R1 files found with suffix: $suffix_r1"
        return 1
    fi
    
    # For each R1, find corresponding R2 and extract sample name
    for r1 in "${r1_files[@]}"; do
        local r2="${r1/${suffix_r1}/${suffix_r2}}"
        
        if [[ ! -f "$r2" ]]; then
            log_warn "Missing R2 file for: $r1"
            continue
        fi
        
        # Extract sample name
        local basename_file
        basename_file=$(basename "$r1")
        local sample="${basename_file%${suffix_r1}}"
        
        FASTQ_R1+=("$r1")
        FASTQ_R2+=("$r2")
        SAMPLE_NAMES+=("$sample")
    done
    
    if [[ ${#SAMPLE_NAMES[@]} -eq 0 ]]; then
        log_error "No valid paired-end samples found"
        return 1
    fi
    
    log_info "Found ${#SAMPLE_NAMES[@]} paired-end samples"
    return 0
}

# ============================= File Utilities =============================

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

safe_create_dir() {
    # Safely create directory with logging
    # Usage: safe_create_dir <directory>
    local dir="$1"
    
    if [[ ! -d "$dir" ]]; then
        mkdir -p "$dir" || die "Failed to create directory: $dir"
        log_info "Created directory: $dir"
    fi
}

# ============================= Process Management =============================

setup_trap() {
    # Setup trap for interrupt signals with cleanup
    # Usage: setup_trap <cleanup_pattern>
    # Example: setup_trap "output_dir/sample*"
    local cleanup_pattern="$1"
    
    trap "echo 'Interrupted! Cleaning up...'; rm -rf $cleanup_pattern; exit 1" SIGINT SIGTERM
}

clear_trap() {
    # Clear interrupt traps
    # Usage: clear_trap
    trap '' SIGINT SIGTERM
}

# ============================= Statistics Functions =============================

extract_salmon_stats() {
    # Extract mapping statistics from Salmon log
    # Usage: extract_salmon_stats <salmon_dir/sample>
    # Returns: "total_reads mapped_reads mapping_rate"
    local salmon_dir="$1"
    local log_file="$salmon_dir/logs/salmon_quant.log"
    
    if [[ ! -f "$log_file" ]]; then
        echo "0 0 0.0"
        return 1
    fi
    
    local total_reads mapped_reads mapping_rate
    
    # Extract total fragments processed
    total_reads=$(grep -oP "Observed \K[0-9,]+" "$log_file" | tr -d ',' || echo "0")
    
    # Extract mapping rate
    mapping_rate=$(grep -oP "Mapping rate = \K[0-9.]+%" "$log_file" | tr -d '%' || echo "0.0")
    
    # Calculate mapped reads
    mapped_reads=$(awk -v total="$total_reads" -v rate="$mapping_rate" \
        'BEGIN{printf "%.0f", total * rate / 100}')
    
    echo "$total_reads $mapped_reads $mapping_rate"
}

extract_star_stats() {
    # Extract mapping statistics from STAR log
    # Usage: extract_star_stats <star_log_file>
    # Returns: "total_reads uniquely_mapped multi_mapped unmapped"
    local log_file="$1"
    
    if [[ ! -f "$log_file" ]]; then
        echo "0 0 0 0"
        return 1
    fi
    
    local total uniquely multi unmapped
    
    total=$(grep "Number of input reads" "$log_file" | awk '{print $NF}')
    uniquely=$(grep "Uniquely mapped reads number" "$log_file" | awk '{print $NF}')
    multi=$(grep "Number of reads mapped to multiple loci" "$log_file" | awk '{print $NF}')
    unmapped=$(grep "Number of reads unmapped: too short" "$log_file" | awk '{print $NF}')
    
    echo "${total:-0} ${uniquely:-0} ${multi:-0} ${unmapped:-0}"
}

# ============================= Index Validation =============================

validate_salmon_index() {
    # Check if Salmon index exists and is valid
    # Usage: validate_salmon_index <index_path>
    local index_path="$1"
    
    if [[ ! -d "$index_path" ]]; then
        log_error "Salmon index not found: $index_path"
        return 1
    fi
    
    # Check for required files
    local required_files=("hash.bin" "pos.bin" "seq.bin" "info.json")
    for file in "${required_files[@]}"; do
        if [[ ! -f "$index_path/$file" ]]; then
            log_error "Incomplete Salmon index (missing $file): $index_path"
            return 1
        fi
    done
    
    log_info "Salmon index validated: $index_path"
    return 0
}

validate_star_index() {
    # Check if STAR index exists and is valid
    # Usage: validate_star_index <index_path>
    local index_path="$1"
    
    if [[ ! -d "$index_path" ]]; then
        log_error "STAR index not found: $index_path"
        return 1
    fi
    
    # Check for required files
    local required_files=("SA" "SAindex" "Genome")
    for file in "${required_files[@]}"; do
        if [[ ! -f "$index_path/$file" ]]; then
            log_error "Incomplete STAR index (missing $file): $index_path"
            return 1
        fi
    done
    
    log_info "STAR index validated: $index_path"
    return 0
}

validate_gtf() {
    # Check if GTF annotation file exists
    # Usage: validate_gtf <gtf_file>
    local gtf_file="$1"
    
    require_file "$gtf_file" "GTF annotation file not found"
    
    # Check if GTF is not empty and has proper format
    if ! head -n1 "$gtf_file" | grep -qE "^#|^chr|^[0-9]"; then
        log_error "Invalid GTF file format: $gtf_file"
        return 1
    fi
    
    log_info "GTF file validated: $gtf_file"
    return 0
}

# ============================= Output File Checking =============================

check_output_exists() {
    # Check if output file/directory exists (for resume functionality)
    # Usage: check_output_exists <path> [<required_size_bytes>]
    # Returns: 0 if exists and valid, 1 otherwise
    local path="$1"
    local min_size="${2:-1}"  # Minimum 1 byte by default
    
    if [[ -e "$path" ]]; then
        if [[ -f "$path" ]]; then
            local size
            size=$(stat -f%z "$path" 2>/dev/null || stat -c%s "$path" 2>/dev/null || echo "0")
            if [[ $size -ge $min_size ]]; then
                return 0
            fi
        elif [[ -d "$path" ]]; then
            return 0
        fi
    fi
    
    return 1
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

# ============================= Validation Functions =============================

validate_species() {
    # Validate species code
    # Usage: validate_species "hs" || die "Invalid species"
    local species="$1"
    
    case "$species" in
        hs|mm|dm|ce)
            return 0
            ;;
        *)
            log_error "Invalid species: $species (must be: hs, mm, dm, ce)"
            return 1
            ;;
    esac
}

validate_library_type() {
    # Validate library strandedness
    # Usage: validate_library_type "IU" || die "Invalid library type"
    local lib_type="$1"
    
    case "$lib_type" in
        IU|ISF|ISR)
            return 0
            ;;
        *)
            log_error "Invalid library type: $lib_type (must be: IU, ISF, ISR)"
            return 1
            ;;
    esac
}

# ============================= Export Functions =============================

# Make functions available to subshells if needed
export -f log log_info log_warn log_error die log_duration
export -f require_file require_dir require_command
export -f get_base_name check_output_exists
export -f validate_species validate_library_type