#!/usr/bin/env bash
# CUT&Tag Analysis Pipeline
# Author: Refactored for modularity and target-specific analysis
# Usage: ./run_cuttag.sh [config.yml]

set -euo pipefail

# ============================= Script Setup =============================

# Get script directory (works with symlinks)
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"

# Source helper functions
HELPER_SCRIPT="${SCRIPT_DIR}/helper.sh"
if [[ ! -f "$HELPER_SCRIPT" ]]; then
    echo "ERROR: Helper script not found: $HELPER_SCRIPT" >&2
    exit 1
fi
source "$HELPER_SCRIPT"

# ============================= Configuration Loading =============================

load_config() {
    local config_file="${1:-config.yml}"
    
    # Default configuration
    SPECIES="chm13"
    SAMPLESHEET="samplesheet.csv"
    OUTDIR="result"
    THREADS=18
    SORT_MEM_LIMIT="16G"
    MAX_FRAG_LENGTH=1000
    
    CALL_PEAK=true
    PEAK_TYPE="narrow"
    RUN_FRIP=true
    RUN_SPIKE=false
    NORM_METHOD="CPM"
    
    INDEX_ROOTDIR="/mnt/f/index"
    
    # Initialize array for custom BED regions
    declare -ga CUSTOM_BED_REGIONS=()
    
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
    
    # Simple YAML parser (handles key: value pairs)
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
            samplesheet) SAMPLESHEET="$value" ;;
            outdir) OUTDIR="$value" ;;
            threads) THREADS="$value" ;;
            sort_mem_limit) SORT_MEM_LIMIT="$value" ;;
            max_frag_length) MAX_FRAG_LENGTH="$value" ;;
            call_peak) CALL_PEAK="$value" ;;
            peak_type) PEAK_TYPE="$value" ;;
            run_frip) RUN_FRIP="$value" ;;
            run_spike) RUN_SPIKE="$value" ;;
            norm_method) NORM_METHOD="$value" ;;
            index_rootdir) INDEX_ROOTDIR="$value" ;;
            bed_region) 
                # Add custom BED regions to array
                if [[ -f "$value" ]]; then
                    CUSTOM_BED_REGIONS+=("$value")
                    log_info "Added custom BED region: $value"
                else
                    log_warn "Custom BED file not found: $value"
                fi
                ;;
        esac
    done < "$config"
}

validate_config() {
    log_info "Validating configuration..."
    
    # Validate species
    validate_species "$SPECIES" || die "Invalid species: $SPECIES"
    
    # Validate peak type
    validate_peak_type "$PEAK_TYPE" || die "Invalid peak type: $PEAK_TYPE"
    
    # Validate normalization method
    case "$NORM_METHOD" in
        CPM|Spike|SpikeFree) ;;
        *) die "Invalid normalization method: $NORM_METHOD (must be: CPM, Spike, SpikeFree)" ;;
    esac
    
    # Validate numeric values
    [[ "$THREADS" -gt 0 ]] || die "Invalid threads: $THREADS"
    [[ "$MAX_FRAG_LENGTH" -gt 0 ]] || die "Invalid max fragment length: $MAX_FRAG_LENGTH"
}

setup_genome_indices() {
    log_info "Setting up genome indices for species: $SPECIES"
    
    local hs_dir="$INDEX_ROOTDIR/hs/v49"
    local chm13_dir="$INDEX_ROOTDIR/hs/chm13"
    local mm_dir="$INDEX_ROOTDIR/mm/vM38"
    local spike_dir="$INDEX_ROOTDIR/Ecoli_novoprotein"
    
    case "$SPECIES" in
        chm13)
            INDEX_DIR="$chm13_dir"
            GSIZE=2913022398
            FA_NAME="chm13v2.0.fa.gz"
            ;;
        hs)
            INDEX_DIR="$hs_dir"
            GSIZE=2913022398
            FA_NAME="GRCh38.primary_assembly.genome.fa.gz"
            ;;
        mm)
            INDEX_DIR="$mm_dir"
            GSIZE=2654621783
            FA_NAME="GRCm39.primary_assembly.genome.fa.gz"
            ;;
    esac
    
    BOWTIE2_INDEX="$INDEX_DIR/bowtie2/bowtie2"
    SPIKE_INDEX="$spike_dir/bowtie2/bowtie2"
    FASTA="$INDEX_DIR/$FA_NAME"
    
    # Ensure FASTA index exists
    CHROM_SIZES=$(ensure_fasta_index "$FASTA")
    
    # Setup default BED regions for heatmaps
    BED_REGIONS=(
        "$INDEX_DIR/genes.bed"
        "$INDEX_DIR/genes_protein_coding.bed"
    )
    
    # Add custom BED regions if specified in config
    if [[ ${#CUSTOM_BED_REGIONS[@]} -gt 0 ]]; then
        log_info "Adding ${#CUSTOM_BED_REGIONS[@]} custom BED region(s) for heatmap generation"
        for bed_file in "${CUSTOM_BED_REGIONS[@]}"; do
            if [[ -f "$bed_file" ]]; then
                BED_REGIONS+=("$bed_file")
                log_info "  - $(basename "$bed_file")"
            else
                log_warn "Custom BED file not found, skipping: $bed_file"
            fi
        done
    fi
    
    # Log all BED regions that will be used
    log_info "Using ${#BED_REGIONS[@]} BED region(s) for heatmap generation:"
    for bed_file in "${BED_REGIONS[@]}"; do
        if [[ -f "$bed_file" ]]; then
            log_info "  ✓ $(basename "$bed_file")"
        else
            log_warn "  ✗ $(basename "$bed_file") (not found, will be skipped)"
        fi
    done
}

setup_directories() {
    log_info "Creating output directories..."
    
    FASTQC_DIR="$OUTDIR/01_fastqc"
    TRIM_DIR="$OUTDIR/02_fastp"
    ALIGN_DIR="$OUTDIR/03_alignment"
    PEAK_DIR="$OUTDIR/04_peaks"
    QC_DIR="$OUTDIR/05_multiqc"
    
    mkdir -p "$OUTDIR" "$FASTQC_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$PEAK_DIR" "$QC_DIR"
    
    # Initialize statistics CSV
    STATS_CSV="$OUTDIR/statistics_${NORM_METHOD}.csv"
    if [[ ! -f "$STATS_CSV" ]]; then
        log_info "Initializing statistics file: $STATS_CSV"
        cat > "$STATS_CSV" <<EOF
sample,group,control,target,fq1,fq2,total_reads,mapped_reads,duplicate_reads,filtered_reads,spike_reads,norm_method,scale_factor,frip
EOF
    fi
}

# ============================= Sample Processing Functions =============================

run_fastqc() {
    local sample="$1" fq1="$2" fq2="$3"
    
    local qc1="$FASTQC_DIR/$(get_base_name "$fq1")_fastqc.zip"
    local qc2="$FASTQC_DIR/$(get_base_name "$fq2")_fastqc.zip"
    
    if [[ ! -s "$qc1" || ! -s "$qc2" ]]; then
        log_info "[$sample] Running FastQC"
        fastqc -t "$THREADS" -o "$FASTQC_DIR" --quiet "$fq1" "$fq2" 2>/dev/null
    fi
}

run_trimming() {
    local sample="$1" fq1="$2" fq2="$3"
    local clean1="$TRIM_DIR/${sample}_R1.fq.gz"
    local clean2="$TRIM_DIR/${sample}_R2.fq.gz"
    
    if [[ ! -s "$clean1" || ! -s "$clean2" ]]; then
        log_info "[$sample] Trimming reads with fastp"
        fastp -i "$fq1" -I "$fq2" -o "$clean1" -O "$clean2" \
            --thread "$THREADS" --detect_adapter_for_pe \
            --json "$TRIM_DIR/${sample}.json" \
            --html "$TRIM_DIR/${sample}.html" 2>/dev/null
    fi
    
    echo "$clean1 $clean2"
}

run_alignment() {
    local sample="$1" clean1="$2" clean2="$3"
    
    local sorted_bam="$ALIGN_DIR/${sample}.sorted.bam"
    local markdup_bam="$ALIGN_DIR/${sample}.markdup.bam"
    local flagstat_file="$ALIGN_DIR/${sample}.markdup.flagstat"
    
    if [[ ! -s "$markdup_bam" || ! -s "${markdup_bam}.bai" ]]; then
        log_info "[$sample] Aligning to genome"
        bowtie2 -x "$BOWTIE2_INDEX" -1 "$clean1" -2 "$clean2" -p "$THREADS" \
            --very-sensitive --local --no-mixed --no-discordant \
            -X "$MAX_FRAG_LENGTH" \
            --rg-id "$sample" --rg "SM:${sample}\tPL:ILLUMINA" \
            2> "$ALIGN_DIR/${sample}.bowtie2.log" \
            | sambamba view -t "$THREADS" -S -f bam /dev/stdin 2>/dev/null \
            | sambamba sort -t "$THREADS" -m "$SORT_MEM_LIMIT" \
                -o "$sorted_bam" /dev/stdin 2>/dev/null
        
        log_info "[$sample] Marking duplicates"
        sambamba markdup -t "$THREADS" "$sorted_bam" "$markdup_bam" 2>/dev/null
        sambamba index -t "$THREADS" "$markdup_bam" 2>/dev/null
        sambamba flagstat "$markdup_bam" > "$flagstat_file" 2>/dev/null
    fi
    
    # Extract statistics
    local total_reads=$(get_flagstat_value "$flagstat_file" "in total" 1)
    local mapped_reads=$(get_flagstat_value "$flagstat_file" "mapped (.*N/A)" 1)
    local duplicate_reads=$(get_flagstat_value "$flagstat_file" "duplicates" 1)
    
    echo "$total_reads $mapped_reads $duplicate_reads"
}

run_spike_alignment() {
    local sample="$1" clean1="$2" clean2="$3"
    
    local spike_sorted_bam="$ALIGN_DIR/${sample}.spike.sorted.bam"
    local spike_markdup_bam="$ALIGN_DIR/${sample}.spike.markdup.bam"
    local spike_flagstat="$ALIGN_DIR/${sample}.spike.markdup.flagstat"
    
    if [[ ! -s "$spike_markdup_bam" ]]; then
        log_info "[$sample] Aligning to spike-in genome"
        bowtie2 -x "$SPIKE_INDEX" -1 "$clean1" -2 "$clean2" -p "$THREADS" \
            --very-sensitive --local --no-mixed --no-discordant \
            -X "$MAX_FRAG_LENGTH" \
            --rg-id "$sample" --rg "SM:${sample}\tPL:ILLUMINA" \
            2> "$ALIGN_DIR/${sample}.spike.bowtie2.log" \
            | sambamba view -t "$THREADS" -S -f bam /dev/stdin 2>/dev/null \
            | sambamba sort -t "$THREADS" -m "$SORT_MEM_LIMIT" \
                -o "$spike_sorted_bam" /dev/stdin 2>/dev/null
        
        sambamba markdup -t "$THREADS" "$spike_sorted_bam" "$spike_markdup_bam" 2>/dev/null
        sambamba index -t "$THREADS" "$spike_markdup_bam" 2>/dev/null
        sambamba flagstat "$spike_markdup_bam" > "$spike_flagstat" 2>/dev/null
    fi
    
    local spike_reads=$(get_flagstat_value "$spike_flagstat" "mapped (.*N/A)" 1)
    echo "$spike_reads"
}

run_filtering() {
    local sample="$1"
    local markdup_bam="$ALIGN_DIR/${sample}.markdup.bam"
    local filtered_bam="$ALIGN_DIR/${sample}.filtered.bam"
    local flagstat_filtered="$ALIGN_DIR/${sample}.filtered.flagstat"
    
    if [[ ! -s "$flagstat_filtered" ]]; then
        log_info "[$sample] Filtering BAM (MAPQ≥10, proper pairs, no chrM/dups)"
        samtools view -@ "$THREADS" -b -f 2 -q 10 -F 1804 \
            -e 'rname != "chrM"' -o "$filtered_bam" "$markdup_bam" 2>/dev/null
        sambamba index -t "$THREADS" "$filtered_bam" 2>/dev/null
        sambamba flagstat -t "$THREADS" "$filtered_bam" > "$flagstat_filtered" 2>/dev/null
    fi
    
    local filtered_reads=$(get_flagstat_value "$flagstat_filtered" "mapped (.*N/A)" 1)
    echo "$filtered_reads"
}

generate_bigwig() {
    local sample="$1" scale_factor="$2"
    local filtered_bam="$ALIGN_DIR/${sample}.filtered.bam"
    local bigwig="$ALIGN_DIR/${sample}.${NORM_METHOD}.bw"
    
    if [[ ! -s "$bigwig" ]]; then
        log_info "[$sample] Generating BigWig (scale: $scale_factor)"
        bamCoverage -b "$filtered_bam" -o "$bigwig" -of bigwig \
            -p "$THREADS" --binSize 10 --scaleFactor "$scale_factor" 2>/dev/null
    fi
}

process_sample_alignment() {
    local sample="$1" group="$2" control="$3" target="$4" fq1="$5" fq2="$6"
    
    log_info ""
    log_info "========== Processing alignment: $sample =========="
    
    # Skip if already processed
    if grep -q "^${sample}," "$STATS_CSV"; then
        log_info "[$sample] Already processed (found in $STATS_CSV)"
        return 0
    fi
    
    # FastQC
    run_fastqc "$sample" "$fq1" "$fq2"
    
    # Trimming
    read clean1 clean2 <<< $(run_trimming "$sample" "$fq1" "$fq2")
    
    # Alignment and marking duplicates
    read total_reads mapped_reads duplicate_reads <<< $(run_alignment "$sample" "$clean1" "$clean2")
    log_info "[$sample] Stats - Total: $total_reads, Mapped: $mapped_reads, Duplicates: $duplicate_reads"
    
    # Spike-in alignment (optional)
    local spike_reads=0
    if [[ "$RUN_SPIKE" == true ]]; then
        spike_reads=$(run_spike_alignment "$sample" "$clean1" "$clean2")
        log_info "[$sample] Spike-in reads: $spike_reads"
    fi
    
    # Filtering
    local filtered_reads=$(run_filtering "$sample")
    log_info "[$sample] Filtered reads: $filtered_reads"
    
    # Calculate normalization scale factor
    local spike_free_file="$ALIGN_DIR/SpikeFree_SF.txt"
    local scale_factor
    scale_factor=$(calculate_scale_factor "$NORM_METHOD" "$filtered_reads" \
        "$spike_reads" "$sample" "$spike_free_file")
    log_info "[$sample] Normalization: $NORM_METHOD, Scale factor: $scale_factor"
    
    # Generate BigWig
    generate_bigwig "$sample" "$scale_factor"
    
    # Append to statistics (frip will be updated later)
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n' \
        "$sample" "$group" "$control" "$target" "$fq1" "$fq2" \
        "$total_reads" "$mapped_reads" "$duplicate_reads" "$filtered_reads" \
        "$spike_reads" "$NORM_METHOD" "$scale_factor" >> "$STATS_CSV"
    
    log_info "[$sample] Alignment completed"
}

# ============================= Peak Calling =============================

call_peaks_for_sample() {
    local sample="$1" group="$2" control="$3" target="$4"
    
    log_info ""
    log_info "========== Calling peaks: $sample =========="
    
    local filtered_bam="$ALIGN_DIR/${sample}.filtered.bam"
    local peak_file="$PEAK_DIR/${sample}_peaks.${PEAK_TYPE}Peak"
    
    if [[ ! -s "$peak_file" ]]; then
        log_info "[$sample] Running MACS3 ($PEAK_TYPE peaks)"
        
        local macs_args=("-t" "$filtered_bam" -f BAMPE -g "$GSIZE" \
            -n "$PEAK_DIR/${sample}" --nomodel --keep-dup all -q 0.1)
        
        [[ "$PEAK_TYPE" == "broad" ]] && macs_args+=("--broad")
        
        if [[ -n "$control" ]]; then
            local control_bam="$ALIGN_DIR/${control}.filtered.bam"
            if [[ -f "$control_bam" ]]; then
                macs_args+=("-c" "$control_bam")
            else
                log_warn "[$sample] Control BAM not found: $control_bam"
            fi
        fi
        
        macs3 callpeak "${macs_args[@]}" 2>/dev/null
    fi
    
    log_info "[$sample] Peak calling completed"
}

create_consensus_peaks() {
    local target="$1"
    local output_prefix="$2"
    shift 2
    local peak_files=("$@")
    
    if [[ ${#peak_files[@]} -eq 0 ]]; then
        log_warn "No peak files provided for consensus: $target"
        return 1
    fi
    
    local consensus_bed="${PEAK_DIR}/${output_prefix}_consensus.bed"
    
    if [[ ! -s "$consensus_bed" ]]; then
        log_info "Creating consensus peaks for: $target (${#peak_files[@]} samples)"
        merge_peaks_to_consensus \
            -i "${peak_files[*]}" \
            -o "$consensus_bed" \
            -g "$CHROM_SIZES" || return 1
    fi
    
    echo "$consensus_bed"
}

# ============================= QC Functions =============================

collect_samples_by_target() {
    # Collect samples organized by target
    # Populates: TARGET_SAMPLES (associative array)
    declare -gA TARGET_SAMPLES
    declare -ga ALL_SAMPLES ALL_BAMS ALL_BWS ALL_TARGETS
    
    while IFS=, read -r sample group control target fq1 fq2 || [[ -n "${sample:-}" ]]; do
        [[ "$sample" == "sample" ]] && continue
        [[ -z "$sample" ]] && continue
        
        local bam="$ALIGN_DIR/${sample}.filtered.bam"
        local bw="$ALIGN_DIR/${sample}.${NORM_METHOD}.bw"
        
        [[ ! -f "$bam" ]] && continue
        
        ALL_SAMPLES+=("$sample")
        ALL_BAMS+=("$bam")
        ALL_BWS+=("$bw")
        
        if [[ -n "$target" ]]; then
            if [[ -z "${TARGET_SAMPLES[$target]:-}" ]]; then
                TARGET_SAMPLES[$target]="$sample"
                ALL_TARGETS+=("$target")
            else
                TARGET_SAMPLES[$target]="${TARGET_SAMPLES[$target]} $sample"
            fi
        fi
    done < <(tail -n +2 "$SAMPLESHEET" | tr -d '\r')
}

run_bam_qc() {
    local -a bam_array=("$@")
    
    if [[ ${#bam_array[@]} -eq 0 ]]; then
        return 0
    fi
    
    local -a name_array=()
    for bam in "${bam_array[@]}"; do
        name_array+=("$(basename "$bam" .filtered.bam)")
    done
    
    # Fragment length distribution
    if [[ ! -s "$QC_DIR/BAM_fragment_length.pdf" ]]; then
        log_info "[QC] Calculating fragment length distribution"
        bamPEFragmentSize -b "${bam_array[@]}" -p "$THREADS" \
            --maxFragmentLength "$MAX_FRAG_LENGTH" \
            --histogram "$QC_DIR/BAM_fragment_length.pdf" \
            --samplesLabel "${name_array[@]}" \
            --outRawFragmentLengths "$QC_DIR/BAM_fragment_length.txt" >/dev/null
    fi
    
    # Fingerprint
    if [[ ! -s "$QC_DIR/BAM_fingerprint.pdf" ]]; then
        log_info "[QC] Generating fingerprint plot"
        plotFingerprint -b "${bam_array[@]}" -p "$THREADS" \
            --labels "${name_array[@]}" \
            --plotFile "$QC_DIR/BAM_fingerprint.pdf" --skipZeros 2>/dev/null
    fi
}

run_bigwig_qc() {
    local prefix="$1"
    shift
    local -a bw_array=("$@")
    
    if [[ ${#bw_array[@]} -eq 0 ]]; then
        return 0
    fi
    
    local -a name_array=()
    for bw in "${bw_array[@]}"; do
        name_array+=("$(basename "$bw" .${NORM_METHOD}.bw)")
    done
    
    local npz_file="$QC_DIR/BigWig_${prefix}_${NORM_METHOD}.npz"
    
    # Compute summary
    if [[ ! -s "$npz_file" ]]; then
        log_info "[QC] Computing BigWig summary: $prefix"
        multiBigwigSummary bins -b "${bw_array[@]}" -p "$THREADS" \
            --labels "${name_array[@]}" -o "$npz_file" >/dev/null
    fi
    
    # PCA plot
    if [[ ! -s "$QC_DIR/BigWig_${prefix}_${NORM_METHOD}_PCA.pdf" ]]; then
        log_info "[QC] Generating PCA plot: $prefix"
        plotPCA -in "$npz_file" \
            -o "$QC_DIR/BigWig_${prefix}_${NORM_METHOD}_PCA.pdf" 2>/dev/null
    fi
    
    # Correlation heatmap
    if [[ ! -s "$QC_DIR/BigWig_${prefix}_${NORM_METHOD}_correlation.pdf" ]]; then
        log_info "[QC] Generating correlation heatmap: $prefix"
        plotCorrelation -in "$npz_file" --corMethod spearman \
            --whatToPlot heatmap \
            -o "$QC_DIR/BigWig_${prefix}_${NORM_METHOD}_correlation.pdf" 2>/dev/null
    fi
}

run_heatmap_qc() {
    local prefix="$1"
    local bed_file="$2"
    shift 2
    local -a bw_array=("$@")
    
    if [[ ${#bw_array[@]} -eq 0 || ! -f "$bed_file" ]]; then
        return 0
    fi
    
    local bed_name=$(basename "$bed_file" .bed)
    local -a name_array=()
    for bw in "${bw_array[@]}"; do
        name_array+=("$(basename "$bw" .${NORM_METHOD}.bw)")
    done
    
    log_info "[QC] Generating heatmaps for: $prefix - $bed_name"
    
    # TSS heatmap
    local heatmap_tss="$QC_DIR/${bed_name}_${prefix}_tss_${NORM_METHOD}_heatmap.pdf"
    if [[ ! -s "$heatmap_tss" ]]; then
        local mat_tss="$QC_DIR/${bed_name}_${prefix}_tss_${NORM_METHOD}_mat.gz"
        [[ ! -s "$mat_tss" ]] && \
            computeMatrix reference-point -S "${bw_array[@]}" -R "$bed_file" \
                -o "$mat_tss" -p "$THREADS" -b 3000 -a 3000 --skipZeros \
                --samplesLabel "${name_array[@]}" --quiet 2>/dev/null
        plotHeatmap -m "$mat_tss" -out "$heatmap_tss" --colorMap Blues 2>/dev/null
    fi
    
    # Gene body heatmap
    local heatmap_body="$QC_DIR/${bed_name}_${prefix}_body_${NORM_METHOD}_heatmap.pdf"
    if [[ ! -s "$heatmap_body" ]]; then
        local mat_body="$QC_DIR/${bed_name}_${prefix}_body_${NORM_METHOD}_mat.gz"
        [[ ! -s "$mat_body" ]] && \
            computeMatrix scale-regions -S "${bw_array[@]}" -R "$bed_file" \
                -o "$mat_body" -p "$THREADS" -b 3000 -a 3000 --skipZeros \
                --regionBodyLength 5000 --samplesLabel "${name_array[@]}" \
                --quiet 2>/dev/null
        plotHeatmap -m "$mat_body" -out "$heatmap_body" --colorMap Blues 2>/dev/null
    fi
}

calculate_frip() {
    local prefix="$1"
    local consensus_bed="$2"
    shift 2
    local -a bam_array=("$@")
    
    if [[ ${#bam_array[@]} -eq 0 || ! -f "$consensus_bed" ]]; then
        return 0
    fi
    
    local -a name_array=()
    for bam in "${bam_array[@]}"; do
        name_array+=("$(basename "$bam" .filtered.bam)")
    done
    
    local frip_tab="$QC_DIR/${prefix}_consensus_frip.tab"
    local frip_plot="$QC_DIR/${prefix}_consensus_frip.pdf"
    
    if [[ ! -s "$frip_tab" ]]; then
        log_info "[QC] Computing FRiP: $prefix"
        plotEnrichment -p "$THREADS" -b "${bam_array[@]}" \
            --BED "$consensus_bed" --labels "${name_array[@]}" \
            --outRawCounts "$frip_tab" --extendReads \
            --plotFile "$frip_plot" \
            --plotTitle "FRiP - $prefix (Consensus Peaks)" 2>/dev/null
    fi
    
    # Update STATS_CSV with FRiP values
    if [[ -s "$frip_tab" ]]; then
        update_frip_in_stats "$frip_tab"
    fi
}

update_frip_in_stats() {
    local frip_file="$1"
    
    # Parse FRiP file and update stats CSV
    while IFS=$'\t' read -r sample frip || [[ -n "$sample" ]]; do
        [[ "$sample" == "sample" ]] && continue
        [[ -z "$sample" ]] && continue
        
        # Remove .filtered.bam suffix if present
        sample=${sample%.filtered.bam}
        
        # Update the FRiP column in stats CSV
        awk -F',' -v s="$sample" -v f="$frip" 'BEGIN{OFS=","} 
            NR==1 {print; next}
            $1==s {$NF=f; print; next}
            {print}
        ' "$STATS_CSV" > "${STATS_CSV}.tmp" && mv "${STATS_CSV}.tmp" "$STATS_CSV"
    done < <(tail -n +2 "$frip_file")
}

run_global_qc() {
    log_info ""
    log_info "========== Running Global QC =========="
    
    # Collect all samples organized by target
    collect_samples_by_target
    
    if [[ ${#ALL_SAMPLES[@]} -eq 0 ]]; then
        log_warn "No samples found for QC"
        return 0
    fi
    
    log_info "Found ${#ALL_SAMPLES[@]} samples across ${#ALL_TARGETS[@]} targets"
    
    # BAM-level QC (all samples)
    run_bam_qc "${ALL_BAMS[@]}"
    
    # BigWig QC for all samples
    run_bigwig_qc "all_samples" "${ALL_BWS[@]}"
    
    # Heatmaps for all samples over genomic regions
    for bed_file in "${BED_REGIONS[@]}"; do
        [[ -f "$bed_file" ]] && run_heatmap_qc "all_samples" "$bed_file" "${ALL_BWS[@]}"
    done
    
    # Target-specific QC
    for target in "${ALL_TARGETS[@]}"; do
        log_info "Processing target-specific QC: $target"
        
        # Get samples for this target
        local -a target_samples target_bams target_bws
        read -ra target_samples <<< "${TARGET_SAMPLES[$target]}"
        
        for sample in "${target_samples[@]}"; do
            target_bams+=("$ALIGN_DIR/${sample}.filtered.bam")
            target_bws+=("$ALIGN_DIR/${sample}.${NORM_METHOD}.bw")
        done
        
        # BigWig QC for target
        run_bigwig_qc "$target" "${target_bws[@]}"
        
        # Heatmaps for target
        for bed_file in "${BED_REGIONS[@]}"; do
            [[ -f "$bed_file" ]] && run_heatmap_qc "$target" "$bed_file" "${target_bws[@]}"
        done
        
        # Consensus peak heatmap and FRiP for target
        local target_consensus="$PEAK_DIR/${target}_consensus.bed"
        if [[ -f "$target_consensus" ]]; then
            run_heatmap_qc "$target" "$target_consensus" "${target_bws[@]}"
            [[ "$RUN_FRIP" == true ]] && calculate_frip "$target" "$target_consensus" "${target_bams[@]}"
        fi
    done
    
    # All-samples consensus peak QC
    local all_consensus="$PEAK_DIR/all_samples_consensus.bed"
    if [[ -f "$all_consensus" ]]; then
        run_heatmap_qc "all_samples" "$all_consensus" "${ALL_BWS[@]}"
        [[ "$RUN_FRIP" == true ]] && calculate_frip "all_samples" "$all_consensus" "${ALL_BAMS[@]}"
    fi
    
    # MultiQC
    log_info "[QC] Generating MultiQC report"
    multiqc "$OUTDIR" --force -o "$QC_DIR" --title "CUT&Tag Report" 2>/dev/null || true
    
    log_info "Global QC completed"
}

# ============================= Consensus Peak Generation =============================

generate_consensus_peaks() {
    log_info ""
    log_info "========== Generating Consensus Peaks =========="
    
    # Collect samples by target
    declare -A target_peaks
    declare -a all_peak_files
    
    while IFS=, read -r sample group control target fq1 fq2 || [[ -n "${sample:-}" ]]; do
        [[ "$sample" == "sample" ]] && continue
        [[ -z "$sample" ]] && continue
        
        local peak_file="$PEAK_DIR/${sample}_peaks.${PEAK_TYPE}Peak"
        [[ ! -f "$peak_file" ]] && continue
        
        all_peak_files+=("$peak_file")
        
        if [[ -n "$target" ]]; then
            if [[ -z "${target_peaks[$target]:-}" ]]; then
                target_peaks[$target]="$peak_file"
            else
                target_peaks[$target]="${target_peaks[$target]} $peak_file"
            fi
        fi
    done < <(tail -n +2 "$SAMPLESHEET" | tr -d '\r')
    
    # Create consensus for each target
    for target in "${!target_peaks[@]}"; do
        local -a peaks
        read -ra peaks <<< "${target_peaks[$target]}"
        create_consensus_peaks "$target" "$target" "${peaks[@]}"
    done
    
    # Create consensus for all samples
    if [[ ${#all_peak_files[@]} -gt 0 ]]; then
        create_consensus_peaks "all_samples" "all_samples" "${all_peak_files[@]}"
    fi
}

# ============================= Main Pipeline =============================

main() {
    local config_file="${1:-}"
    
    # Setup logging
    local log_file="./pipeline.log"
    exec > >(tee -a "$log_file") 2>&1
    
    log_info "============================================================"
    log_info "CUT&Tag Analysis Pipeline"
    log_info "============================================================"
    
    # Load configuration
    load_config "$config_file"
    
    # Display configuration
    log_info "Configuration:"
    log_info "  Species: $SPECIES"
    log_info "  Samplesheet: $SAMPLESHEET"
    log_info "  Output directory: $OUTDIR"
    log_info "  Threads: $THREADS"
    log_info "  Normalization: $NORM_METHOD"
    log_info "  Peak calling: $CALL_PEAK ($PEAK_TYPE)"
    log_info "  Spike-in: $RUN_SPIKE"
    log_info "  FRiP calculation: $RUN_FRIP"
    log_info "============================================================"
    
    # Setup genome indices
    setup_genome_indices
    
    # Setup directories
    setup_directories
    
    # Check samplesheet
    if [[ ! -s "$SAMPLESHEET" ]]; then
        log_warn "Samplesheet not found: $SAMPLESHEET"
        log_info "Generating template samplesheet..."
        generate_samplesheet "." "$SAMPLESHEET"
        log_info "Please edit $SAMPLESHEET and re-run the pipeline"
        exit 0
    fi
    
    # Phase 1: Alignment
    log_info ""
    log_info "========== Phase 1: Alignment =========="
    safe_read_csv "$SAMPLESHEET" process_sample_alignment
    
    # Phase 2: Peak Calling
    if [[ "$CALL_PEAK" == true ]]; then
        log_info ""
        log_info "========== Phase 2: Peak Calling =========="
        safe_read_csv "$SAMPLESHEET" call_peaks_for_sample
        
        # Phase 3: Consensus Peaks
        generate_consensus_peaks
    fi
    
    # Phase 4: QC
    run_global_qc
    
    # Summary
    log_info ""
    log_info "============================================================"
    log_info "Pipeline completed successfully!"
    log_info "============================================================"
    log_info "Results:"
    log_info "  Statistics → $STATS_CSV"
    log_info "  Alignments → $ALIGN_DIR/"
    log_info "  Peaks → $PEAK_DIR/"
    log_info "  QC Reports → $QC_DIR/"
    if [[ "$CALL_PEAK" == true ]]; then
        log_info "  Consensus peaks:"
        find "$PEAK_DIR" -name "*_consensus.bed" -exec echo "    {}" \;
    fi
    log_info "  Log file → $log_file"
    log_info "============================================================"
}

# ============================= Entry Point =============================

# Run main pipeline
main "$@"