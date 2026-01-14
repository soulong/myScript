# CUT&Tag Analysis Pipeline

A modular and robust pipeline for analyzing CUT&Tag sequencing data with target-specific consensus peak calling and comprehensive QC.

## Features

- ✅ Automated FastQC, trimming, alignment, and peak calling
- ✅ Multiple normalization methods (CPM, Spike-in, SpikeFree)
- ✅ Target-specific and global consensus peak generation
- ✅ Comprehensive QC reports with deepTools and MultiQC
- ✅ FRiP (Fraction of Reads in Peaks) calculation
- ✅ Modular design with reusable helper functions
- ✅ Configurable via YAML or command-line defaults

## Requirements

### Software Dependencies

- **Alignment**: bowtie2, samtools, sambamba
- **Quality Control**: fastqc, fastp, multiqc
- **Peak Calling**: macs3
- **Analysis**: deeptools, bedtools
- **System**: bash ≥4.0, awk, grep, find

### Installation

All tools should be available in your PATH. We recommend using conda:

```bash
conda create -n cuttag \
    bowtie2 samtools sambamba fastqc fastp multiqc \
    macs3 deeptools bedtools -c bioconda -c conda-forge
conda activate cuttag
```

## Quick Start

### 1. Prepare Your Data

Organize your FASTQ files in a directory:

```
project/
├── sample1_R1.fq.gz
├── sample1_R2.fq.gz
├── sample2_R1.fq.gz
├── sample2_R2.fq.gz
└── ...
```

### 2. Generate Samplesheet

Run the pipeline once to auto-generate a template:

```bash
./run_cuttag.sh
```

This creates `samplesheet.csv`. Edit it to add sample information:

```csv
sample,group,control,target,fq1,fq2
H3K4me3_rep1,treatment,,H3K4me3,./sample1_R1.fq.gz,./sample1_R2.fq.gz
H3K4me3_rep2,treatment,,H3K4me3,./sample2_R1.fq.gz,./sample2_R2.fq.gz
H3K27me3_rep1,treatment,,H3K27me3,./sample3_R1.fq.gz,./sample3_R2.fq.gz
IgG_rep1,control,,IgG,./sample4_R1.fq.gz,./sample4_R2.fq.gz
```

**Columns:**
- `sample`: Unique sample identifier
- `group`: Experimental group (e.g., treatment, control)
- `control`: Control sample name for peak calling (optional)
- `target`: Antibody target (e.g., H3K4me3, H3K27ac)
- `fq1`, `fq2`: Paths to R1 and R2 FASTQ files

### 3. Create Configuration (Optional)

Create `config.yml` in your working directory:

```yaml
species: chm13
threads: 24
peak_type: narrow
norm_method: CPM
index_rootdir: /path/to/your/indices
```

### 4. Run Pipeline

```bash
# With default configuration
./run_cuttag.sh

# With custom configuration
./run_cuttag.sh config.yml
```

The pipeline can be located anywhere. It will work in your current directory.

## Pipeline Workflow

### Phase 1: Alignment

1. **FastQC**: Quality control of raw reads
2. **Trimming**: Adapter trimming with fastp
3. **Alignment**: Bowtie2 alignment to reference genome
4. **Filtering**: Remove duplicates, low-quality reads, mitochondrial reads
5. **Normalization**: Generate normalized BigWig tracks
6. **Spike-in** (optional): Align to spike-in genome for normalization

### Phase 2: Peak Calling

1. **Individual Peaks**: MACS3 peak calling per sample
2. **Target Consensus**: Merge peaks within each target group
3. **Global Consensus**: Merge all peaks across all samples

### Phase 3: Quality Control

1. **BAM QC**: Fragment size distribution, fingerprint analysis
2. **BigWig QC**: PCA, correlation heatmaps (target-specific and global)
3. **Genomic Feature Heatmaps**: Signal at TSS and gene bodies
4. **Peak Heatmaps**: Signal at consensus peaks (target-specific and global)
5. **FRiP Calculation**: Enrichment in consensus peaks
6. **MultiQC**: Aggregate report

## Output Structure

```
result/
├── 01_fastqc/              # FastQC reports
├── 02_fastp/               # Trimming reports and clean reads
├── 03_alignment/           # BAM files, BigWig tracks
│   ├── sample.filtered.bam
│   ├── sample.CPM.bw
│   └── ...
├── 04_peaks/               # Peak files
│   ├── sample_peaks.narrowPeak
│   ├── H3K4me3_consensus.bed         # Target-specific consensus
│   ├── H3K27me3_consensus.bed
│   └── all_samples_consensus.bed     # Global consensus
├── 05_multiqc/             # QC reports
│   ├── multiqc_report.html
│   ├── BigWig_H3K4me3_CPM_PCA.pdf    # Target-specific QC
│   ├── BigWig_all_samples_CPM_PCA.pdf
│   └── ...
├── statistics_CPM.csv      # Summary statistics
└── pipeline.log            # Execution log
```

## Configuration Options

### Species Options

- `hs`: Human GRCh38
- `chm13`: Human T2T-CHM13 (default)
- `mm`: Mouse GRCm39

### Peak Types

- `narrow`: Sharp peaks (transcription factors, H3K4me3)
- `broad`: Broad peaks (histone marks like H3K27me3, H3K36me3)

### Normalization Methods

- `CPM`: Counts per million (default)
- `Spike`: Spike-in normalization (requires `run_spike: true`)
- `SpikeFree`: Use pre-computed scale factors from file

## Target-Specific Analysis

The refactored pipeline now performs analysis at two levels:

### 1. Target-Specific

For each unique target in the `target` column:
- Separate consensus peaks
- Target-specific BigWig QC (PCA, correlation)
- Target-specific heatmaps
- Target-specific FRiP calculation

### 2. Global (All Samples)

- Combined consensus peaks from all samples
- Global PCA and correlation analysis
- Global heatmaps and FRiP

This allows you to assess biological replicates within each target and compare across different targets.

## Helper Functions

The `helper.sh` script provides reusable functions:

### Logging
```bash
log_info "Processing sample"
log_warn "Missing file"
log_error "Critical error"
die "Fatal error, exiting"
```

### File Validation
```bash
require_file "/path/to/file"
require_dir "/path/to/dir"
require_command "samtools"
```

### Statistics Extraction
```bash
reads=$(get_flagstat_value "sample.flagstat" "mapped" 1)
```

### Peak Merging
```bash
merge_peaks_to_consensus \
    -i "peak1.bed peak2.bed peak3.bed" \
    -o consensus.bed \
    -g genome.sizes
```

### Normalization
```bash
scale=$(calculate_scale_factor "CPM" "$filtered_reads")
```

## Troubleshooting

### Issue: "Helper script not found"

Ensure `helper.sh` is in the same directory as `run_cuttag.sh`, or the script directory is in your PATH.

### Issue: "Index files not found"

Update `index_rootdir` in your config file to point to your genome index directory.

### Issue: Samplesheet parsing errors

- Ensure no spaces around commas
- Check for DOS line endings (use `dos2unix` if needed)
- Verify all file paths are correct

### Issue: Peak calling fails

- Ensure you have enough filtered reads (>1M recommended)
- Check MACS3 log files in `result/04_peaks/`
- Verify control samples are correctly specified

## Advanced Usage

### Resume Interrupted Pipeline

The pipeline automatically resumes from the last completed step. Simply re-run:

```bash
./run_cuttag.sh
```

Completed samples (present in `statistics_CPM.csv`) are skipped.

### Custom Genome Regions

Add custom BED files to your `config.yml` for additional heatmaps:

```yaml
# Add one or more custom BED regions
bed_region: /path/to/enhancers.bed
bed_region: /path/to/promoters.bed
bed_region: /path/to/custom_peaks.bed
```

Each BED file should be tab-delimited with at least 3 columns (chr, start, end):

```
chr1    1000    2000    enhancer1
chr1    5000    6000    enhancer2
chr2    3000    4000    enhancer3
```

The pipeline will automatically generate heatmaps for:
- **Default regions**: genes.bed, genes_protein_coding.bed (from genome index)
- **Custom regions**: All BED files specified in config
- **Analysis levels**: Both target-specific and all-samples

**Example workflow:**

1. Create custom BED files:
```bash
# Extract promoter regions (TSS ± 2kb)
awk 'BEGIN{OFS="\t"} {
    if($6=="+") print $1, ($2-2000<0?0:$2-2000), $2+2000, $4, ".", $6
    else print $1, ($3-2000<0?0:$3-2000), $3+2000, $4, ".", $6
}' genes.bed > promoters.bed
```

2. Add to config.yml:
```yaml
bed_region: ./promoters.bed
bed_region: ./enhancers.bed
```

3. Run pipeline - it will generate heatmaps for all specified regions!

### Spike-in Normalization

```yaml
run_spike: true
norm_method: Spike
```

Ensure spike-in genome index is available at `$index_rootdir/Ecoli_novoprotein/`.

## Citation

If you use this pipeline, please cite the original tools:

- **Bowtie2**: Langmead & Salzberg, Nat Methods 2012
- **MACS3**: Zhang et al., Genome Biol 2008; updated in 2023
- **deepTools**: Ramírez et al., Nucleic Acids Res 2016
- **MultiQC**: Ewels et al., Bioinformatics 2016

## License

This pipeline is provided as-is for research use.

## Support

For issues or questions:
1. Check the log file: `pipeline.log`
2. Verify all dependencies are installed
3. Ensure input files and indices are correctly formatted