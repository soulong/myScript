#!/usr/bin/env bash
# RNA-seq Analysis Pipeline
# Performs: QC, trimming, pseudo-alignment (Salmon), alignment (STAR), quantification
# Usage: ./run_rnaseq.sh [config.yml]

set -euo pipefail

# ============================= Script Setup =============================

# Get script directory (works with symlinks)
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"

# Source helper functions
HELPER_SCRIPT="${SCRIPT_DIR}/rnaseq_helper.sh"
if [[ ! -f "$HELPER_SCRIPT" ]]; then
    echo "ERROR: Helper script not found: $HELPER_SCRIPT" >&2
    exit 1
fi
source "$HELPER_SCRIPT"

# ============================= Configuration Loading =============================

load_config() {
    local config_file="${1:-config.yml}"
    
    # Default configuration
    SPECIES="hs"
    FASTQ_DIR="01.RawData"
    OUTDIR="result"
    THREADS=28
    
    # Analysis options
    RUN_STAR=false
    RUN_FEATURECOUNTS=false
    
    # Index paths
    INDEX_ROOTDIR="/mnt/f/index"
    
    # Library type for featureCounts: IU=0 (unstranded), ISF=1 (forward), ISR=2 (reverse)
    LIBRARY_TYPE="IU"
    
    # If config file exists, parse it
    if [[ -f "$config_file" ]]; then
        log_info "Loading configuration from: $config_file"
        parse_yaml_config "$config_file"
    else
        log_info "Using default configuration (no config file provided)"
    fi
    
    # Validate configuration
    validate_config
}

parse_yaml_config() {
    local config="$1"
    
    # Simple YAML parser
    while IFS=': ' read -r key value || [[ -n "$key" ]]; do
        # Skip comments and empty lines
        [[ "$key" =~ ^[[:space:]]*# ]] && continue
        [[ -z "$key" ]] && continue
        
        # Remove leading/trailing whitespace and quotes
        key=$(echo "$key" | xargs)
        value=$(echo "$value" | xargs | sed 's/^["'\'']\|["'\'']$//g')
        
        # Skip empty values
        [[ -z "$value" ]] && continue
        
        # Map config keys to variables
        case "$key" in
            species) SPECIES="$value" ;;
            fastq_dir) FASTQ_DIR="$value" ;;
            outdir) OUTDIR="$value" ;;
            threads) THREADS="$value" ;;
            run_star) RUN_STAR="$value" ;;
            run_featurecounts) RUN_FEATURECOUNTS="$value" ;;
            index_rootdir) INDEX_ROOTDIR="$value" ;;
            library_type) LIBRARY_TYPE="$value" ;;
        esac
    done < "$config"
}

validate_config() {
    log_info "Validating configuration..."
    
    # Validate species
    validate_species "$SPECIES" || die "Invalid species: $SPECIES"
    
    # Validate library type
    validate_library_type "$LIBRARY_TYPE" || die "Invalid library type: $LIBRARY_TYPE"
    
    # Validate numeric values
    [[ "$THREADS" -gt 0 ]] || die "Invalid threads: $THREADS"
    
    # Check FASTQ directory exists
    require_dir "$FASTQ_DIR" "FASTQ directory not found: $FASTQ_DIR"
}

setup_indices() {
    log_info "Setting up genome indices for species: $SPECIES"
    
    case "$SPECIES" in
        hs)
            INDEX_DIR="$INDEX_ROOTDIR/hs/v49"
            SALMON_INDEX="$INDEX_DIR/salmon"
            STAR_INDEX="$INDEX_DIR/star"
            GTF_FILE="$INDEX_DIR/gencode.v49.annotation.gtf"
            ;;
        mm)
            INDEX_DIR="$INDEX_ROOTDIR/mm/vM38"
            SALMON_INDEX="$INDEX_DIR/salmon"
            STAR_INDEX="$INDEX_DIR/star"
            GTF_FILE="$INDEX_DIR/gencode.vM38.annotation.gtf"
            ;;
    esac
    
    # Validate Salmon index
    validate_salmon_index "$SALMON_INDEX" || die "Invalid Salmon index"
    
    # Validate STAR index if needed
    if [[ "$RUN_STAR" == true ]]; then
        validate_star_index "$STAR_INDEX" || die "Invalid STAR index"
    fi
    
    # Validate GTF if featureCounts is needed
    if [[ "$RUN_FEATURECOUNTS" == true ]]; then
        validate_gtf "$GTF_FILE" || die "Invalid GTF file"
    fi
}

setup_directories() {
    log_info "Creating output directories..."
    
    FASTQC_DIR="$OUTDIR/01_fastqc"
    FASTP_DIR="$OUTDIR/02_fastp"
    SALMON_DIR="$OUTDIR/03_salmon"
    STAR_DIR="$OUTDIR/04_star"
    COUNTS_DIR="$OUTDIR/05_counts"
    MULTIQC_DIR="$OUTDIR/06_multiqc"
    
    safe_create_dir "$OUTDIR"
    safe_create_dir "$FASTQC_DIR"
    safe_create_dir "$FASTP_DIR"
    safe_create_dir "$SALMON_DIR"
    
    if [[ "$RUN_STAR" == true ]]; then
        safe_create_dir "$STAR_DIR"
    fi
    
    if [[ "$RUN_FEATURECOUNTS" == true ]]; then
        safe_create_dir "$COUNTS_DIR"
    fi
    
    safe_create_dir "$MULTIQC_DIR"
    
    # Initialize statistics CSV
    STATS_CSV="$OUTDIR/statistics.csv"
    if [[ ! -f "$STATS_CSV" ]]; then
        log_info "Initializing statistics file: $STATS_CSV"
        cat > "$STATS_CSV" <<EOF
sample,raw_r1,raw_r2,clean_r1,clean_r2,salmon_total,salmon_mapped,salmon_rate,star_total,star_unique,star_multi,star_unmapped
EOF
    fi
}

# ============================= Sample Processing Functions =============================

run_fastqc_analysis() {
    local sample="$1"
    local r1="$2"
    local r2="$3"
    
    local qc1="$FASTQC_DIR/$(get_base_name "$r1")_fastqc.zip"
    local qc2="$FASTQC_DIR/$(get_base_name "$r2")_fastqc.zip"
    
    if check_output_exists "$qc1" && check_output_exists "$qc2"; then
        log_info "[$sample] FastQC outputs exist, skipping"
        return 0
    fi
    
    log_info "[$sample] Running FastQC"
    local start_time=$SECONDS
    
    fastqc -t "$THREADS" -o "$FASTQC_DIR" --quiet "$r1" "$r2" 2>/dev/null || {
        log_error "[$sample] FastQC failed"
        return 1
    }
    
    log_duration $((SECONDS - start_time)) "$sample FastQC"
}

run_fastp_trimming() {
    local sample="$1"
    local r1="$2"
    local r2="$3"
    
    local clean_r1="$FASTP_DIR/${sample}_R1.fastq.gz"
    local clean_r2="$FASTP_DIR/${sample}_R2.fastq.gz"
    
    if check_output_exists "$clean_r1" 1000 && check_output_exists "$clean_r2" 1000; then
        log_info "[$sample] Trimmed FASTQ files exist, skipping"
        echo "$clean_r1 $clean_r2"
        return 0
    fi
    
    log_info "[$sample] Running fastp trimming"
    local start_time=$SECONDS
    
    setup_trap "$FASTP_DIR/${sample}*"
    
    fastp -w "$THREADS" \
        -i "$r1" -I "$r2" \
        -o "$clean_r1" -O "$clean_r2" \
        -j "$FASTP_DIR/${sample}.json" \
        -h "$FASTP_DIR/${sample}.html" \
        -R "$sample" \
        2>> "$FASTP_DIR/fastp.log" || {
        log_error "[$sample] fastp failed"
        clear_trap
        return 1
    }
    
    clear_trap
    log_duration $((SECONDS - start_time)) "$sample fastp"
    
    echo "$clean_r1 $clean_r2"
}

run_salmon_quant() {
    local sample="$1"
    local clean_r1="$2"
    local clean_r2="$3"
    
    local salmon_out="$SALMON_DIR/$sample"
    
    if check_output_exists "$salmon_out/quant.sf"; then
        log_info "[$sample] Salmon quantification exists, skipping"
        return 0
    fi
    
    log_info "[$sample] Running Salmon quantification"
    local start_time=$SECONDS
    
    setup_trap "$salmon_out*"
    
    salmon quant -i "$SALMON_INDEX" -l A -p "$THREADS" \
        -1 "$clean_r1" -2 "$clean_r2" \
        --seqBias --gcBias --validateMappings \
        -o "$salmon_out" \
        2>> "$SALMON_DIR/salmon.log" || {
        log_error "[$sample] Salmon failed"
        clear_trap
        return 1
    }
    
    clear_trap
    log_duration $((SECONDS - start_time)) "$sample Salmon"
}

run_star_alignment() {
    local sample="$1"
    local clean_r1="$2"
    local clean_r2="$3"
    
    local star_prefix="$STAR_DIR/${sample}_"
    local star_bam="${star_prefix}Aligned.sortedByCoord.out.bam"
    
    if check_output_exists "$star_bam" 1000000; then
        log_info "[$sample] STAR alignment exists, skipping"
        return 0
    fi
    
    log_info "[$sample] Running STAR alignment"
    local start_time=$SECONDS
    
    setup_trap "${star_prefix}*"
    
    STAR --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesCommand zcat \
        --readFilesIn "$clean_r1" "$clean_r2" \
        --outFileNamePrefix "$star_prefix" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        2>> "$STAR_DIR/star.log" || {
        log_error "[$sample] STAR failed"
        clear_trap
        return 1
    }
    
    # Index BAM file
    if [[ -f "$star_bam" ]]; then
        log_info "[$sample] Indexing BAM file"
        samtools index -@ "$THREADS" "$star_bam" || log_warn "[$sample] BAM indexing failed"
    fi
    
    clear_trap
    log_duration $((SECONDS - start_time)) "$sample STAR"
}

extract_sample_statistics() {
    local sample="$1"
    local r1="$2"
    local r2="$3"
    local clean_r1="$4"
    local clean_r2="$5"
    
    # Count raw reads
    local raw_r1=0 raw_r2=0
    if [[ -f "$r1" ]]; then
        raw_r1=$(($(zcat "$r1" | wc -l) / 4))
    fi
    if [[ -f "$r2" ]]; then
        raw_r2=$(($(zcat "$r2" | wc -l) / 4))
    fi
    
    # Count clean reads
    local clean_r1_count=0 clean_r2_count=0
    if [[ -f "$clean_r1" ]]; then
        clean_r1_count=$(($(zcat "$clean_r1" | wc -l) / 4))
    fi
    if [[ -f "$clean_r2" ]]; then
        clean_r2_count=$(($(zcat "$clean_r2" | wc -l) / 4))
    fi
    
    # Extract Salmon statistics
    local salmon_stats
    read salmon_total salmon_mapped salmon_rate <<< \
        $(extract_salmon_stats "$SALMON_DIR/$sample")
    
    # Extract STAR statistics if available
    local star_total=0 star_unique=0 star_multi=0 star_unmapped=0
    if [[ "$RUN_STAR" == true ]]; then
        local star_log="$STAR_DIR/${sample}_Log.final.out"
        if [[ -f "$star_log" ]]; then
            read star_total star_unique star_multi star_unmapped <<< \
                $(extract_star_stats "$star_log")
        fi
    fi
    
    # Append to statistics CSV
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$sample" "$raw_r1" "$raw_r2" \
        "$clean_r1_count" "$clean_r2_count" \
        "$salmon_total" "$salmon_mapped" "$salmon_rate" \
        "$star_total" "$star_unique" "$star_multi" "$star_unmapped" \
        >> "$STATS_CSV"
}

process_sample() {
    local sample="$1"
    local r1="$2"
    local r2="$3"
    
    log_info ""
    log_info "========== Processing sample: $sample =========="
    
    # Skip if already in stats
    if grep -q "^${sample}," "$STATS_CSV"; then
        log_info "[$sample] Already processed (found in $STATS_CSV)"
        return 0
    fi
    
    # Step 1: FastQC on raw reads
    run_fastqc_analysis "$sample" "$r1" "$r2"
    
    # Step 2: Trimming with fastp
    local clean_r1 clean_r2
    read clean_r1 clean_r2 <<< $(run_fastp_trimming "$sample" "$r1" "$r2")
    
    if [[ ! -f "$clean_r1" || ! -f "$clean_r2" ]]; then
        log_error "[$sample] Trimmed files not found, skipping remaining steps"
        return 1
    fi
    
    # Step 3: Salmon pseudo-alignment
    run_salmon_quant "$sample" "$clean_r1" "$clean_r2"
    
    # Step 4: STAR alignment (optional)
    if [[ "$RUN_STAR" == true ]]; then
        run_star_alignment "$sample" "$clean_r1" "$clean_r2"
    fi
    
    # Extract and save statistics
    extract_sample_statistics "$sample" "$r1" "$r2" "$clean_r1" "$clean_r2"
    
    log_info "[$sample] Processing completed"
}

# ============================= FeatureCounts Quantification =============================

run_featurecounts() {
    if [[ "$RUN_FEATURECOUNTS" != true ]]; then
        return 0
    fi
    
    log_info ""
    log_info "========== Running featureCounts =========="
    
    # Collect all BAM files
    local -a bam_files
    mapfile -t bam_files < <(find "$STAR_DIR" -name "*Aligned.sortedByCoord.out.bam" | sort)
    
    if [[ ${#bam_files[@]} -eq 0 ]]; then
        log_warn "No BAM files found for featureCounts"
        return 0
    fi
    
    log_info "Found ${#bam_files[@]} BAM files"
    
    local counts_file="$COUNTS_DIR/gene_counts.txt"
    
    if check_output_exists "$counts_file"; then
        log_info "featureCounts output exists, skipping"
        return 0
    fi
    
    # Convert library type to featureCounts strandness parameter
    local strandness=0
    case "$LIBRARY_TYPE" in
        IU) strandness=0 ;;
        ISF) strandness=1 ;;
        ISR) strandness=2 ;;
    esac
    
    log_info "Running featureCounts (strandness: $strandness)"
    local start_time=$SECONDS
    
    featureCounts -T "$THREADS" -s "$strandness" \
        -t exon -g gene_id -Q 20 -p -C --donotsort \
        -a "$GTF_FILE" \
        -o "$counts_file" "${bam_files[@]}" \
        2>> "$COUNTS_DIR/featurecounts.log" || {
        log_error "featureCounts failed"
        return 1
    }
    
    log_duration $((SECONDS - start_time)) "featureCounts"
}

# ============================= MultiQC Report =============================

run_multiqc() {
    log_info ""
    log_info "========== Running MultiQC =========="
    
    if check_output_exists "$MULTIQC_DIR/multiqc_report.html"; then
        log_info "MultiQC report exists, regenerating..."
    fi
    
    local start_time=$SECONDS
    
    multiqc -f "$OUTDIR" -o "$MULTIQC_DIR" \
        --title "RNA-seq Analysis Report" \
        2>/dev/null || {
        log_warn "MultiQC failed (non-critical)"
        return 0
    }
    
    log_duration $((SECONDS - start_time)) "MultiQC"
}

# ============================= Main Pipeline =============================

main() {
    local config_file="${1:-}"
    
    # Setup logging
    local log_file="./rnaseq_pipeline.log"
    exec > >(tee -a "$log_file") 2>&1
    
    log_info "============================================================"
    log_info "RNA-seq Analysis Pipeline"
    log_info "============================================================"
    
    # Check dependencies
    check_dependencies || die "Please install missing dependencies"
    
    # Load configuration
    load_config "$config_file"
    
    # Display configuration
    log_info ""
    log_info "Configuration:"
    log_info "  Species: $SPECIES"
    log_info "  FASTQ directory: $FASTQ_DIR"
    log_info "  Output directory: $OUTDIR"
    log_info "  Threads: $THREADS"
    log_info "  Library type: $LIBRARY_TYPE"
    log_info "  Run STAR: $RUN_STAR"
    log_info "  Run featureCounts: $RUN_FEATURECOUNTS"
    log_info "============================================================"
    
    # Setup indices
    setup_indices
    
    # Setup directories
    setup_directories
    
    # Detect FASTQ suffixes
    log_info ""
    log_info "========== Detecting FASTQ Files =========="
    local suffix_r1 suffix_r2
    read suffix_r1 suffix_r2 <<< $(detect_fastq_suffix "$FASTQ_DIR")
    log_info "Detected R1 suffix: $suffix_r1"
    log_info "Detected R2 suffix: $suffix_r2"
    
    # Collect paired FASTQ files
    collect_paired_fastq "$FASTQ_DIR" "$suffix_r1" "$suffix_r2"
    
    log_info ""
    log_info "Paired-end samples detected:"
    printf '  - %s\n' "${SAMPLE_NAMES[@]}"
    log_info ""
    
    # Process each sample
    for i in "${!SAMPLE_NAMES[@]}"; do
        process_sample "${SAMPLE_NAMES[$i]}" "${FASTQ_R1[$i]}" "${FASTQ_R2[$i]}"
    done
    
    # Run featureCounts if requested
    run_featurecounts
    
    # Generate MultiQC report
    run_multiqc
    
    # Summary
    log_info ""
    log_info "============================================================"
    log_info "Pipeline completed successfully!"
    log_info "============================================================"
    log_info "Results:"
    log_info "  Statistics → $STATS_CSV"
    log_info "  FastQC → $FASTQC_DIR/"
    log_info "  Trimmed reads → $FASTP_DIR/"
    log_info "  Salmon quantification → $SALMON_DIR/"
    if [[ "$RUN_STAR" == true ]]; then
        log_info "  STAR alignments → $STAR_DIR/"
    fi
    if [[ "$RUN_FEATURECOUNTS" == true ]]; then
        log_info "  Gene counts → $COUNTS_DIR/gene_counts.txt"
    fi
    log_info "  MultiQC report → $MULTIQC_DIR/multiqc_report.html"
    log_info "  Log file → $log_file"
    log_info "============================================================"
}

# ============================= Entry Point =============================

# Run main pipeline
main "$@"