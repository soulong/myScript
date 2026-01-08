#!/usr/bin/env bash
set -euo pipefail

# ============================= Global Configuration =============================
readonly SPECIES="chm13"                # "hs", "chm13" or "mm"
readonly SAMPLESHEET="samplesheet_CutTag.csv"
readonly OUTDIR="result"
readonly THREADS=18
readonly SORT_MEM_LIMIT="16G"        # sambamba sort MEM limit

# Max insert length
readonly MAX_FRAG_LENGTH=1000

# Peak calling config
readonly CALL_PEAK=true             # MACS2
readonly PEAK_TYPE="narrow"           # "broad" or "narrow"
readonly RUN_FRIP=true

# Spike-in
readonly RUN_SPIKE=true

# Normalization
readonly NORM_METHOD="Spike"         # CPM | Spike | SpikeFree


# Index paths (use absolute or relative; avoid ~ in scripts)
# readonly INDEX_ROOTDIR="${HOME}/Index"
readonly INDEX_ROOTDIR="/mnt/f/index"
readonly HS_INDEX_DIR="$INDEX_ROOTDIR/hs/v49"
readonly CHM13_INDEX_DIR="$INDEX_ROOTDIR/hs/chm13"
readonly MM_INDEX_DIR="$INDEX_ROOTDIR/mm/vM38"
readonly SPIKE_DIR="$INDEX_ROOTDIR/Ecoli_novoprotein"

# Output subdirs
readonly FASTQCDIR="$OUTDIR/01_fastqc"
readonly TRIMDIR="$OUTDIR/02_fastp"
readonly ALIGNDIR="$OUTDIR/03_alignment"
readonly PEAKDIR="$OUTDIR/04_peaks"
readonly QCDIR="$OUTDIR/05_multiqc"

# Spike-free scale factor file (must exist if NORM_METHOD=SpikeFree)
readonly SPIKE_FREE_SF_FILE="$ALIGNDIR/SpikeFree_SF.txt"

# Genome setup
case "$SPECIES" in
  chm13)
    readonly INDEX_DIR="$CHM13_INDEX_DIR"
    readonly GSIZE=2913022398
    readonly FA_NAME="chm13v2.0.fa.gz"
    ;;
  hs)
    readonly INDEX_DIR="$HS_INDEX_DIR"
    readonly GSIZE=2913022398
    readonly FA_NAME="GRCh38.primary_assembly.genome.fa.gz"
    ;;
  mm)
    readonly INDEX_DIR="$MM_INDEX_DIR"
    readonly GSIZE=2654621783
    readonly FA_NAME="GRCm39.primary_assembly.genome.fa.gz"
    ;;
  *)
    echo "ERROR: SPECIES must be 'hs' 'chm13' 'mm'"
    exit 1
    ;;
esac

readonly BOWTIE2_INDEX="$INDEX_DIR/bowtie2/bowtie2"
readonly SPIKE_INDEX="$SPIKE_DIR/bowtie2/bowtie2"
readonly FASTA="$INDEX_DIR/$FA_NAME"

readonly CHROM_SIZES="${FASTA%.gz}.fai"
# Create .fai if missing (samtools faidx requires uncompressed FASTA)
if [[ ! -f "$CHROM_SIZES" ]]; then
  echo "Creating chromosome sizes file: $CHROM_SIZES"
  tmp_fa="${FASTA%.gz}"
  gzip -d -k "$FASTA" && samtools faidx "$tmp_fa" && rm -f "$tmp_fa"
fi

# Pre-defined BED regions for heatmaps
BED_REGIONS=(
  "$INDEX_DIR/genes.bed"
  "$INDEX_DIR/genes_protein_coding.bed"
  # "$INDEX_DIR/genes_lncRNA.bed"
)


# ============================= Helper Functions =============================
log()   { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
die()   { echo "ERROR: $*"; exit 1; }
require_file() { [[ -s "$1" ]] || die "Missing/empty required file: $1"; }

# Safely read samplesheet (skip header, handle DOS line endings)
safe_loop() {
  local file="$1"; shift
  require_file "$file"
  while IFS=, read -r sample group control fq1 fq2 || [[ -n "${sample:-}" ]]; do
    [[ "$sample" == "sample" ]] && continue
    [[ -z "$sample" ]] && continue
    "$@" "$sample" "$group" "$control" "$fq1" "$fq2"
  done < <(tail -n +2 "$file" | tr -d '\r')
}

# Extract value from samtools flagstat output
get_flagstat_value() {
  local file="$1" pattern="$2" col="${3:-1}"
  awk -v p="$pattern" -v c="$col" '
    $0 ~ p { print $c; exit }
  ' "$file"
}

merge_consensus_peak() {
    # Merge all individual peak files into one consensus BED (union of all peaks)
    # Options:
    #   -p <glob>   Input peak files pattern (default: *_peaks.*Peak)
    #   -o <bed file> Output file (default: consensus_peaks.bed)
    #   -G <file>   Genome chrom.sizes file (required for proper sorting)

    local peak_glob="*_peaks.*Peak"
    local outfile="consensus_peaks.bed"
    local genome=""

    while getopts "p:o:G:" opt; do
        case $opt in
            p) peak_glob="$OPTARG" ;;
            o) outfile="$OPTARG" ;;
            G) genome="$OPTARG" ;;
            *) echo "Usage: $0 [-p peak_glob] [-o output_prefix] [-G genome_file]" >&2; return 1 ;;
        esac
    done

    [[ -z "$genome" ]] && { echo "Error: -G genome file is required" >&2; return 1; }

    local sample_peaks=($peak_glob)

    if [[ ${#sample_peaks[@]} -eq 0 || "${sample_peaks[0]}" == "$peak_glob" ]]; then
        echo "Error: No peak files found matching '$peak_glob'" >&2
        return 1
    fi

    echo "Merging ${#sample_peaks[@]} peak files into $outfile"

    # Concatenate all files → extract first 3 cols → sort → merge within samples → merge across → sort by genome
    cat "${sample_peaks[@]}" | \
        cut -f1-3 | \
        bedtools sort | \
        bedtools merge | \
        bedtools sort -g "$genome" > "$outfile"

    if [[ -s "$outfile" ]]; then
        echo "Consensus peaks written to: $outfile ($(wc -l < "$outfile") regions)"
    else
        echo "Error: No consensus peaks generated."
        rm -f "$outfile"
        return 1
    fi
}

get_base_fastq_name() {
    local file="$1"
    local name=$(basename "$file")
    # Remove known extensions
    name=${name%.fastq.gz}
    name=${name%.fq.gz}
    echo "$name"
}

# ============================= Process Single Sample (Alignment Phase) =============================
process_alignment() {
  local sample="$1" group="$2" control="$3" fq1="$4" fq2="$5"

  log ""
  log ">>>>>>>>>>>>>>> Processing alignment for sample: $sample <<<<<<<<<<<<<<<<<"

  # Resume check
  if grep -q "^${sample}," "$STATS_CSV"; then
    log "[SKIP] Already processed (found in $STATS_CSV)"
    return 0
  fi

  # Input validation
  # require_file "$fq1"
  # require_file "$fq2"

  # File paths
  local clean1="$TRIMDIR/${sample}_R1.fq.gz"
  local clean2="$TRIMDIR/${sample}_R2.fq.gz"
  local sorted_bam="$ALIGNDIR/${sample}.sorted.bam"
  local markdup_bam="$ALIGNDIR/${sample}.markdup.bam"
  local flagstat_file="$ALIGNDIR/${sample}.markdup.flagstat"
  local filtered_bam="$ALIGNDIR/${sample}.filtered.bam"
  local flagstat_filtered_file="$ALIGNDIR/${sample}.filtered.flagstat"
  local frag_bw="$ALIGNDIR/${sample}.${NORM_METHOD}.bw"

  # FastQC
  local qc1="$FASTQCDIR/$(get_base_fastq_name "$fq1")_fastqc.zip"
  local qc2="$FASTQCDIR/$(get_base_fastq_name "$fq2")_fastqc.zip"
  if [[ ! -s "$qc1" || ! -s "$qc2" ]]; then
    log "[FastQC] $sample"
    fastqc -t "$THREADS" -o "$FASTQCDIR" --quiet "$fq1" "$fq2" 2>/dev/null
  fi

  # Trimming
  if [[ ! -s "$clean1" ]] || [[ ! -s "$clean2" ]]; then
    log "[Trimming] $sample"
    fastp -i "$fq1" -I "$fq2" -o "$clean1" -O "$clean2" \
      --thread "$THREADS" --detect_adapter_for_pe \
      --json "$TRIMDIR/${sample}.json" --html "$TRIMDIR/${sample}.html" 2>/dev/null
  fi

  # Alignment → Sort → Mark Duplicates
  if [[ ! -s "$markdup_bam" ]] || [[ ! -s "${markdup_bam}.bai" ]]; then
    log "[Alignment] $sample -> $INDEX_DIR"
    bowtie2 -x "$BOWTIE2_INDEX" -1 "$clean1" -2 "$clean2" -p "$THREADS" \
      --very-sensitive --local --no-mixed --no-discordant -X "$MAX_FRAG_LENGTH" \
      --rg-id "$sample" --rg "SM:${sample}\tPL:ILLUMINA" \
      2> "$ALIGNDIR/${sample}.bowtie2.log" \
      | sambamba view -t "$THREADS" -S -f bam /dev/stdin 2>/dev/null \
      | sambamba sort -t "$THREADS" -m "$SORT_MEM_LIMIT" -o $sorted_bam /dev/stdin 2>/dev/null

    log "[Mark Duplicates] $sample"
    sambamba markdup -t "$THREADS" "$sorted_bam" "$markdup_bam" 2>/dev/null
    sambamba index -t "$THREADS" "$markdup_bam" 2>/dev/null
    sambamba flagstat "$markdup_bam" > "$ALIGNDIR/${sample}.markdup.flagstat" 2>/dev/null
  fi
  # Statistics from markdup BAM
  local total_reads=$(get_flagstat_value "$flagstat_file" "in total" 1)
  local mapped_reads=$(get_flagstat_value "$flagstat_file" "mapped (.*N/A)" 1)
  local duplicate_reads=$(get_flagstat_value "$flagstat_file" "duplicates" 1)
  log "[Stats] total=$total_reads, mapped=$mapped_reads, dup=$duplicate_reads"

  # Spike-in alignment (optional)
  local spike_reads=0
  if [[ "$RUN_SPIKE" == true ]]; then
    local spike_sorted_bam="$ALIGNDIR/${sample}.spike.sorted.bam"
    local spike_markdup_bam="$ALIGNDIR/${sample}.spike.markdup.bam"
    local spike_flagstat_file="$ALIGNDIR/${sample}.spike.markdup.flagstat"
    if [[ ! -s "$spike_markdup_bam" ]]; then
      log "[Spike-in Alignment] $sample -> $SPIKE_DIR"
      bowtie2 -x "$SPIKE_INDEX" -1 "$clean1" -2 "$clean2" -p "$THREADS" \
        --very-sensitive --local --no-mixed --no-discordant -X "$MAX_FRAG_LENGTH" \
        --rg-id "$sample" --rg "SM:${sample}\tPL:ILLUMINA" \
        2> "$ALIGNDIR/${sample}.spike.bowtie2.log" \
        | sambamba view -t "$THREADS" -S -f bam /dev/stdin 2>/dev/null \
        | sambamba sort -t "$THREADS" -m "$SORT_MEM_LIMIT" -o "$spike_sorted_bam" /dev/stdin 2>/dev/null

      log "[Mark Duplicates (Spike)] $sample"
      sambamba markdup -t "$THREADS" "$spike_sorted_bam" "$spike_markdup_bam" 2>/dev/null
      sambamba index -t "$THREADS" "$spike_markdup_bam" 2>/dev/null
      sambamba flagstat "$spike_markdup_bam" > "$spike_flagstat_file" 2>/dev/null
    fi
    spike_reads=$(get_flagstat_value "$spike_flagstat_file" "mapped (.*N/A)" 1)
    log "[Stats] spike_reads=$spike_reads"
  fi

  # Filter BAM (proper pairs, MAPQ≥10, no chrM, no dups)
  if [[ ! -s "${flagstat_filtered_file}" ]]; then
    log "[Filter] $sample"
    samtools view -@ "$THREADS" -b -f 2 -q 10 -F 1804 -e 'rname != "chrM"' -o "$filtered_bam" "$markdup_bam" 2>/dev/null
    sambamba index -t "$THREADS" "$filtered_bam" 2>/dev/null
    sambamba flagstat -t "$THREADS" "$filtered_bam" >"$flagstat_filtered_file" 2>/dev/null
  fi
  local filtered_reads=$(get_flagstat_value "$flagstat_filtered_file" "mapped (.*N/A)" 1)
  log "[Stats] filtered_reads=$filtered_reads"

  # Normalization scale factor
  local scale_factor=1.0
  if [[ "$NORM_METHOD" == "CPM" && $filtered_reads -gt 0 ]]; then
    scale_factor=$(awk "BEGIN{printf \"%.6f\", 1e6 / $filtered_reads}")
  elif [[ "$NORM_METHOD" == "Spike" && $spike_reads -gt 0 ]]; then
    scale_factor=$(awk "BEGIN{printf \"%.6f\", 1e4 / $spike_reads}")
  elif [[ "$NORM_METHOD" == "SpikeFree" ]]; then
    if [[ ! -s "$SPIKE_FREE_SF_FILE" ]]; then
      die "SpikeFree mode requires file: $SPIKE_FREE_SF_FILE"
    fi
    local sf_raw
    sf_raw=$(awk -F'\t' -v s="$sample" '
      NR>1 {
        gsub(/\.filtered\.bam$/, "", $1)
        if ($1 == s) { print $7; exit }
      }
    ' "$SPIKE_FREE_SF_FILE")
    if [[ -z "$sf_raw" ]]; then
      die "Sample $sample not found in $SPIKE_FREE_SF_FILE"
    fi
    scale_factor=$(awk "BEGIN{printf \"%.6f\", 1e6 / ($filtered_reads * $sf_raw)}")
  fi
  log "[Normalization] method=$NORM_METHOD, scale_factor=$scale_factor"

  # Generate BigWig
  if [[ ! -s "$frag_bw" ]]; then
    log "[BigWig] $sample"
    bamCoverage -b "$filtered_bam" -o "$frag_bw" -of bigwig -p "$THREADS" \
      --binSize 10 --scaleFactor "$scale_factor" 2>/dev/null
  fi

  # Append to stats CSV (with placeholder for frip)
  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%.6f,\n' \
    "$sample" "$group" "$control" "$fq1" "$fq2" \
    "$total_reads" "$mapped_reads" "$duplicate_reads" "$filtered_reads" "$spike_reads" \
    "$NORM_METHOD" "$scale_factor" \
    >> "$STATS_CSV"

  log "[COMPLETED] Alignment for $sample"
}


# ============================= Process Peaks =============================
process_peaks() {
  local sample="$1" group="$2" control="$3" fq1="$4" fq2="$5"

  log ""
  log ">>>>>>>>>>>>>>> Processing peaks for sample: $sample <<<<<<<<<<<<<<<<<"

  local filtered_bam="$ALIGNDIR/${sample}.filtered.bam"
  local peak_file="$PEAKDIR/${sample}_peaks.${PEAK_TYPE}Peak"

  # Peak calling with MACS2
  if [[ "$CALL_PEAK" == true ]]; then
    if [[ ! -s "$peak_file" ]]; then
      log "[Peak Calling] MACS2 ($PEAK_TYPE) on $sample"
      mkdir -p "$PEAKDIR"
      local macs_args=("-t" "$filtered_bam" -f BAMPE -g "$GSIZE" \
        -n "$PEAKDIR/${sample}" --nomodel --keep-dup all -q 0.2)
      if [[ "$PEAK_TYPE" == "broad" ]]; then
        macs_args+=("--broad")
      fi
      if [[ -n "$control" ]]; then
        local control_bam="$ALIGNDIR/${control}.filtered.bam"
        require_file "$control_bam"
        macs_args+=("-c" "$control_bam")
      fi
      # echo macs2 callpeak "${macs_args[@]}"
      macs2 callpeak "${macs_args[@]}" 2>/dev/null
    fi
  fi

  log "[COMPLETED] Peaks for $sample"
}


# ============================= Global QC =============================
run_qc() {
  local samplesheet_path="$1"
  local consensus_peak="$2"

  log ""
  log ">>>>>>>>>>>>>>> Global QC <<<<<<<<<<<<<<<<<"

  # Collect files in order (or auto-detect)
  local -a bam_array=() name_array=() bw_array=()

  if [[ -n "$samplesheet_path" && -s "$samplesheet_path" ]]; then
    while IFS=, read -r sample _ _ _ _ || [[ -n "${sample:-}" ]]; do
      [[ "$sample" == "sample" ]] && continue
      [[ -z "$sample" ]] && continue
      local bam="$ALIGNDIR/${sample}.filtered.bam"
      local bw="$ALIGNDIR/${sample}.${NORM_METHOD}.bw"
      if [[ -s "$bam" ]]; then
        bam_array+=("$bam")
        name_array+=("$sample")
        bw_array+=("$bw")
      fi
    done < <(tail -n +2 "$samplesheet_path" | tr -d '\r')
  else
    mapfile -t bams < <(find "$ALIGNDIR" -name "*.filtered.bam" | sort)
    for bam in "${bams[@]}"; do
      sample=$(basename "$bam" .filtered.bam)
      bw="$ALIGNDIR/${sample}.${NORM_METHOD}.bw"
      bam_array+=("$bam")
      name_array+=("$sample")
      bw_array+=("$bw")
    done
  fi

  # BAM-level QC
  if (( ${#bam_array[@]} > 0 )); then
    log "[QC] Found ${#bam_array[@]} BAM files"

    # Fragment length
    if [[ ! -s "$QCDIR/BAM_fragment_length.pdf" ]]; then
      log "[QC] BAM Fragment length"
      bamPEFragmentSize -b "${bam_array[@]}" -p "$THREADS" --maxFragmentLength "$MAX_FRAG_LENGTH" \
        --histogram "$QCDIR/BAM_fragment_length.pdf" --samplesLabel "${name_array[@]}" \
        --outRawFragmentLengths "$QCDIR/BAM_fragment_length.txt" >/dev/null
    fi
    # Fingerprint
    log "[QC] BAM Fingerprint"
    if [[ ! -s "$QCDIR/BAM_fingerprint.pdf" ]]; then
      plotFingerprint -b "${bam_array[@]}" -p "$THREADS" --labels "${name_array[@]}" \
        --plotFile "$QCDIR/BAM_fingerprint.pdf" --skipZeros
    fi
  fi

  # BigWig-level QC
  if (( ${#bw_array[@]} > 0 )); then
    log "[QC] Found ${#bw_array[@]} BigWig files (method: $NORM_METHOD)"

    # multiBigwigSummary + Correlation/PCA
    if [[ ! -s "$QCDIR/BigWig_${NORM_METHOD}.npz" ]]; then
      log "[QC] BigWig calculation"
      multiBigwigSummary bins -b "${bw_array[@]}" -p "$THREADS" --labels "${name_array[@]}" \
        -o "$QCDIR/BigWig_${NORM_METHOD}.npz" >/dev/null
    fi
    if [[ ! -s "$QCDIR/BigWig_${NORM_METHOD}_PCA.pdf" ]]; then
      log "[QC] BigWig PCA"
      plotPCA -in "$QCDIR/BigWig_${NORM_METHOD}.npz" -o "$QCDIR/BigWig_${NORM_METHOD}_PCA.pdf"
    fi
    if [[ ! -s "$QCDIR/BigWig_${NORM_METHOD}_correlation.pdf" ]]; then
      log "[QC] BigWig correlation"
      plotCorrelation -in "$QCDIR/BigWig_${NORM_METHOD}.npz" --corMethod spearman \
        --whatToPlot heatmap -o "$QCDIR/BigWig_${NORM_METHOD}_correlation.pdf"
    fi

    # Heatmaps over genes
    local heatmap_color='Blues'
    for bed_path in "${BED_REGIONS[@]}"; do
      [[ ! -f "$bed_path" ]] && continue
      bed_name=$(basename "$bed_path" .bed)
      log "[QC] Heatmap over bed: $bed_name"

      # TSS
      heatmap_tss="$QCDIR/${bed_name}_tss_${NORM_METHOD}_heatmap.pdf"
      if [[ ! -s "$heatmap_tss" ]]; then
        log "[QC] TSS"
        mat_tss="$QCDIR/${bed_name}_tss_${NORM_METHOD}_mat.gz"
        [[ ! -s "$mat_tss" ]] && computeMatrix reference-point -S "${bw_array[@]}" -R "$bed_path" -o "$mat_tss" \
          -p "$THREADS" -b 3000 -a 3000 --skipZeros --samplesLabel "${name_array[@]}" --quiet
        plotHeatmap -m "$mat_tss" -out "$heatmap_tss" --colorMap "$heatmap_color"
      fi

      heatmap_body="$QCDIR/${bed_name}_body_${NORM_METHOD}_heatmap.pdf"
      if [[ ! -s "$heatmap_body" ]]; then
        mat_body="$QCDIR/${bed_name}_body_${NORM_METHOD}_mat.gz"
        log "[QC] gene body"
        [[ ! -s "$mat_body" ]] && computeMatrix scale-regions -S "${bw_array[@]}" -R "$bed_path" -o "$mat_body" \
          -p "$THREADS" -b 3000 -a 3000 --skipZeros --samplesLabel "${name_array[@]}" --quiet --regionBodyLength 5000
        plotHeatmap -m "$mat_body" -out "$heatmap_body" --colorMap "$heatmap_color"
      fi
    done

    # Consensus peak heatmaps and FRiP
    if [[ -n "$consensus_peak" && -s "$consensus_peak" ]]; then
      heatmap_peak="$QCDIR/consensus_peak_${NORM_METHOD}_heatmap.pdf"
      if [[ ! -s "$heatmap_peak" ]]; then
        log "[QC] Heatmap over consensus peak $consensus_peak"
        mat_peak="$QCDIR/consensus_peak_${NORM_METHOD}_mat.gz"
        [[ ! -s "$mat_peak" ]] && computeMatrix reference-point --referencePoint center -S "${bw_array[@]}" -R "$consensus_peak" \
          -o "$mat_peak" -p "$THREADS" -b 3000 -a 3000 --skipZeros --samplesLabel "${name_array[@]}" \
          --outFileSortedRegions "$QCDIR/consensus_peak_${NORM_METHOD}_sorted.bed"
        plotHeatmap -m "$mat_peak" -out "$heatmap_peak" --colorMap "$heatmap_color"
      fi

      # FRiP with plotEnrichment using consensus peaks
      local frip_tab="$QCDIR/consensus_frip.tab"
      local frip_plot="$QCDIR/consensus_frip.pdf"
      if [[ ! -s "$frip_tab" ]]; then
        log "[QC] Computing FRiP using plotEnrichment with consensus peaks"
        plotEnrichment -p "$THREADS" -b "${bam_array[@]}" --BED "$consensus_peak" \
          --labels "${name_array[@]}" --outRawCounts "$frip_tab" --extendReads \
          --plotFile "$frip_plot" --plotTitle "FRiP (Consensus Peaks)"
      fi

    fi

  fi

  # MultiQC
  log "[MultiQC] Generating report"
  multiqc "$OUTDIR" --force -o "$QCDIR" --title "CUT&Tag Report"
}


# ============================= Main Execution =============================

# Logging
LOG_FILE="./pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

log "============================================================"
log "Starting CUT&Tag pipeline"
log "Species: $SPECIES | Spike-in: $RUN_SPIKE | Normalization: $NORM_METHOD"
log "Peak calling: $CALL_PEAK ($PEAK_TYPE) | FRiP: $RUN_FRIP"
log "Log file: $LOG_FILE"
log "============================================================"

# Initialize directories and stats CSV
mkdir -p "$OUTDIR" "$FASTQCDIR" "$TRIMDIR" "$ALIGNDIR" "$PEAKDIR" "$QCDIR"
readonly STATS_CSV="$OUTDIR/statistics_${NORM_METHOD}.csv"
if [[ ! -f "$STATS_CSV" ]]; then
  log "Initializing statistics CSV: $STATS_CSV"
  cat > "$STATS_CSV" <<EOF
sample,group,control,fq1,fq2,total_reads,mapped_reads,duplicate_reads,filtered_reads,spike_reads,norm_method,scale_${NORM_METHOD}
EOF
fi

# Alignment phase
safe_loop "$SAMPLESHEET" process_alignment

# Peak calling phase (individual sample peaks)
if [[ "$CALL_PEAK" == true ]]; then
  safe_loop "$SAMPLESHEET" process_peaks
fi

# Consensus peak calling
consensus_peakfile="$PEAKDIR/consensus.bed"
if [[ ! -s "$consensus_peakfile" ]]; then
  merge_consensus_peak -p "$PEAKDIR/*_peaks.${PEAK_TYPE}Peak" -o "$consensus_peakfile" -G "$CHROM_SIZES"
fi

# Global QC (including consensus FRiP)
run_qc "$SAMPLESHEET" "$consensus_peakfile"

log "============================================================"
log "Pipeline finished successfully!"
log "Statistics → $STATS_CSV"
log "Alignments → $ALIGNDIR/"
log "QC Reports → $QCDIR/"
log "Consensus peaks → $consensus_peakfile"
log "============================================================"