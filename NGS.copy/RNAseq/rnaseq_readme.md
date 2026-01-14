# RNA-seq Analysis Pipeline

A modular, robust, and user-friendly pipeline for RNA-seq data analysis with automatic paired-end detection, quality control, pseudo-alignment, and optional full alignment with quantification.

## Features

- ✅ Automatic paired-end FASTQ file detection
- ✅ Comprehensive quality control (FastQC)
- ✅ Adapter trimming and quality filtering (fastp)
- ✅ Fast transcript quantification (Salmon)
- ✅ Optional genome alignment (STAR)
- ✅ Optional gene-level quantification (featureCounts)
- ✅ Integrated quality reports (MultiQC)
- ✅ Resume capability (skip completed samples)
- ✅ Detailed logging and error handling
- ✅ Configurable via YAML or defaults
- ✅ Standalone execution from any location

## Quick Start

### 1. Install Dependencies

All tools should be available in your PATH:

```bash
conda create -n rnaseq \
    fastqc fastp salmon star samtools subread multiqc \
    -c bioconda -c conda-forge
conda activate rnaseq
```

### 2. Prepare Your Data

Organize FASTQ files in a directory:

```
project/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
├── sample2_R2.fastq.gz
└── ...
```

**Supported naming conventions:**
- `_R1.fastq.gz` / `_R2.fastq.gz`
- `_r1.fq.gz` / `_r2.fq.gz`
- `_1.fastq.gz` / `_2.fastq.gz`

### 3. Create Configuration (Optional)

Create `config.yml` in your working directory:

```yaml
species: hs
fastq_dir: 01.RawData
threads: 24
run_star: false
run_featurecounts: false
library_type: IU
index_rootdir: /path/to/your/indices
```

### 4. Run Pipeline

```bash
# With default configuration
./run_rnaseq.sh

# With custom configuration
./run_rnaseq.sh config.yml
```

The pipeline can be located anywhere and will work in your current directory.

## Pipeline Workflow

### Phase 1: Quality Control & Trimming

1. **FASTQ Detection**: Automatically finds paired-end files
2. **FastQC**: Quality assessment of raw reads
3. **fastp**: Adapter trimming, quality filtering, duplication removal

### Phase 2: Quantification

4. **Salmon**: Fast pseudo-alignment and transcript quantification
   - Automatic library type detection
   - Sequence and GC bias correction
   - Selective alignment validation

### Phase 3: Alignment (Optional)

5. **STAR**: Genome alignment if `run_star: true`
   - Produces sorted BAM files
   - Generates comprehensive mapping statistics

### Phase 4: Gene Quantification (Optional)

6. **featureCounts**: Gene-level read counting if `run_featurecounts: true`
   - Requires STAR BAM files
   - Handles various library types

### Phase 5: Reporting

7. **MultiQC**: Aggregated quality control report
8. **Statistics Summary**: CSV file with all metrics

## Output Structure

```
result/
├── 01_fastqc/                    # FastQC reports
│   ├── sample1_R1_fastqc.html
│   └── ...
├── 02_fastp/                     # Trimming reports and clean reads
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample1.json
│   └── fastp.log
├── 03_salmon/                    # Salmon quantification
│   ├── sample1/
│   │   ├── quant.sf             # Transcript abundance
│   │   ├── quant.genes.sf       # Gene abundance
│   │   └── logs/
│   └── salmon.log
├── 04_star/                      # STAR alignments (if enabled)
│   ├── sample1_Aligned.sortedByCoord.out.bam
│   ├── sample1_Log.final.out
│   └── star.log
├── 05_counts/                    # featureCounts (if enabled)
│   ├── gene_counts.txt
│   └── featurecounts.log
├── 06_multiqc/                   # MultiQC report
│   └── multiqc_report.html
├── statistics.csv                # Summary statistics
└── rnaseq_pipeline.log          # Pipeline execution log
```

## Configuration Options

### Species Options

- `hs`: Human (Homo sapiens)
- `mm`: Mouse (Mus musculus)
- `dm`: Fruit fly (Drosophila melanogaster)
- `ce`: Worm (Caenorhabditis elegans)

### Library Types

For featureCounts strandedness:

- `IU`: Unstranded (default) - Used for standard RNA-seq
- `ISF`: Forward stranded - Used for some Illumina protocols
- `ISR`: Reverse stranded - Most common for modern libraries (dUTP-based)

**How to determine your library type:**

```bash
# Use RSeQC infer_experiment.py or check your library prep protocol
# Common protocols:
# - TruSeq Stranded mRNA: ISR
# - SMARTer Stranded: ISF
# - Standard TruSeq: IU
```

### Analysis Modes

#### Mode 1: Salmon Only (Fast, Recommended)
```yaml
run_star: false
run_featurecounts: false
```
- Fastest option
- Provides transcript and gene-level counts
- Sufficient for most differential expression analysis
- **Use this for standard RNA-seq analysis**

#### Mode 2: Salmon + STAR (With BAM files)
```yaml
run_star: true
run_featurecounts: false
```
- Generates BAM files for visualization (IGV, UCSC)
- Use if you need aligned reads for manual inspection
- Salmon still provides quantification

#### Mode 3: Full Pipeline (STAR + featureCounts)
```yaml
run_star: true
run_featurecounts: true
```
- Traditional RNA-seq workflow
- Gene counts from featureCounts
- Use if you specifically need featureCounts output
- Slower than Salmon alone

## Understanding the Outputs

### Salmon Quantification Files

**`quant.sf`** - Transcript-level abundance:
```
Name                Length  EffectiveLength  TPM      NumReads
ENST00000456328.2   1657    1515.000        0.000000 0.000
ENST00000450305.2   632     490.000         0.358097 1.000
```

**`quant.genes.sf`** - Gene-level abundance (aggregated):
```
Name               Length   EffectiveLength  TPM      NumReads
ENSG00000223972.5  1735.67  1593.67         0.000000 0.000
ENSG00000227232.5  1351.00  1209.00         0.297564 2.000
```

**Key columns:**
- **TPM**: Transcripts Per Million (normalized, comparable across samples)
- **NumReads**: Estimated fragment count

### Statistics CSV

Comprehensive per-sample metrics:

- `raw_r1/r2`: Number of reads in raw FASTQ files
- `clean_r1/r2`: Number of reads after trimming
- `salmon_total`: Total fragments processed by Salmon
- `salmon_mapped`: Number of mapped fragments
- `salmon_rate`: Mapping rate percentage
- `star_total`: Total reads input to STAR (if run)
- `star_unique`: Uniquely mapped reads
- `star_multi`: Multi-mapped reads
- `star_unmapped`: Unmapped reads

### MultiQC Report

Interactive HTML report with:
- FastQC summary (per-base quality, adapter content, duplication)
- fastp trimming statistics
- Salmon mapping rates and bias metrics
- STAR alignment metrics (if run)
- Sample correlation and PCA (if multiple samples)

## Genome Indices Setup

### Directory Structure

```
/path/to/indices/
├── hs/
│   ├── salmon/
│   ├── star/
│   └── gencode.v46.annotation.gtf
├── mm/
│   ├── salmon/
│   ├── star/
│   └── gencode.vM35.annotation.gtf
└── ...
```

### Building Salmon Index

```bash
# Human transcriptome
salmon index \
    -t gencode.v46.transcripts.fa.gz \
    -i salmon_index \
    -p 8 \
    --gencode

# Mouse transcriptome
salmon index \
    -t gencode.vM35.transcripts.fa.gz \
    -i salmon_index \
    -p 8 \
    --gencode
```

**Download transcriptomes:**
- Human: https://www.gencodegenes.org/human/
- Mouse: https://www.gencodegenes.org/mouse/

### Building STAR Index

```bash
# Human genome
STAR --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile gencode.v46.annotation.gtf \
    --sjdbOverhang 100 \
    --runThreadN 8

# Mouse genome
STAR --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles GRCm39.primary_assembly.genome.fa \
    --sjdbGTFfile gencode.vM35.annotation.gtf \
    --sjdbOverhang 100 \
    --runThreadN 8
```

**Download genomes:**
- Human: https://www.gencodegenes.org/human/
- Mouse: https://www.gencodegenes.org/mouse/

## Downstream Analysis

### Importing Salmon Results to R

```R
# Install tximport and DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("tximport", "DESeq2", "GenomicFeatures"))

# Load libraries
library(tximport)
library(DESeq2)

# Create sample metadata
samples <- data.frame(
    sample = c("sample1", "sample2", "sample3", "sample4"),
    condition = c("control", "control", "treatment", "treatment")
)

# Import Salmon quantification
files <- file.path("result/03_salmon", samples$sample, "quant.sf")
names(files) <- samples$sample

# Load transcript-to-gene mapping
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("gencode.v46.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Import with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, 
                                 colData = samples,
                                 design = ~ condition)

# Run differential expression
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "treatment", "control"))
```

### Using featureCounts Output

```R
# Read featureCounts output
counts <- read.table("result/05_counts/gene_counts.txt", 
                     header = TRUE, row.names = 1, skip = 1)

# Remove annotation columns and keep only counts
counts <- counts[, 6:ncol(counts)]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = samples,
                               design = ~ condition)
```

## Troubleshooting

### Issue: "No gzipped FASTQ files found"

**Solution**: Ensure FASTQ files are gzipped and in the correct directory

```bash
# Check files
ls 01.RawData/*.gz

# Compress if needed
gzip 01.RawData/*.fastq
```

### Issue: Salmon index validation failed

**Solution**: Check index directory structure

```bash
# Verify index files exist
ls /path/to/salmon_index/
# Should contain: hash.bin, pos.bin, seq.bin, info.json

# Rebuild if corrupted
salmon index -t transcripts.fa -i salmon_index
```

### Issue: Low mapping rates

**Possible causes:**
1. **Wrong species**: Check if using correct genome
2. **Low quality reads**: Inspect FastQC reports
3. **Contamination**: Check for adapter sequences
4. **Library type mismatch**: For featureCounts, verify strandedness

**Debug:**
```bash
# Check Salmon log
cat result/03_salmon/salmon.log

# Inspect a sample's detailed log
cat result/03_salmon/sample1/logs/salmon_quant.log
```

### Issue: Pipeline interrupted

**Solution**: Simply re-run - the pipeline will resume

```bash
# Resume from where it stopped
./run_rnaseq.sh config.yml
```

Completed samples (in `statistics.csv`) are automatically skipped.

### Issue: Out of memory

**Solutions:**
1. Reduce number of threads
2. Increase system memory
3. Process fewer samples at once

```yaml
# Reduce threads in config
threads: 8
```

## Performance Considerations

### Computational Requirements

**Per sample (30M reads, 2x100bp):**
- **FastQC**: ~2-5 minutes
- **fastp**: ~5-10 minutes
- **Salmon**: ~5-10 minutes (very fast!)
- **STAR**: ~20-30 minutes (slower)
- **featureCounts**: ~5 minutes

**Memory:**
- **Salmon**: ~8-12 GB
- **STAR**: ~30-35 GB (human genome)
- **fastp/FastQC**: ~2-4 GB

**Storage (per sample):**
- Raw FASTQ: ~3-5 GB
- Trimmed FASTQ: ~2-4 GB
- BAM file (if STAR): ~4-6 GB
- Salmon output: ~100-200 MB

### Optimization Tips

1. **Use Salmon only** for standard analysis (fastest)
2. **Run STAR only if needed** for BAM files
3. **Increase threads** for faster processing
4. **Process on compute cluster** for large cohorts

## Best Practices

1. **Always check MultiQC report** after pipeline completion
2. **Inspect low mapping rate samples** - may indicate quality issues
3. **Use Salmon for differential expression** - it's accurate and fast
4. **Keep raw FASTQ files** - disk is cheap, re-sequencing is expensive
5. **Document your library type** - critical for featureCounts
6. **Version control your config** - reproducibility is key

## Citation

If you use this pipeline, please cite the tools:

- **FastQC**: Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data
- **fastp**: Chen et al., Bioinformatics 2018
- **Salmon**: Patro et al., Nat Methods 2017
- **STAR**: Dobin et al., Bioinformatics 2013
- **featureCounts**: Liao et al., Bioinformatics 2014
- **MultiQC**: Ewels et al., Bioinformatics 2016

## Support

For issues:
1. Check `rnaseq_pipeline.log` for detailed errors
2. Verify all dependencies are installed: `check_dependencies`
3. Ensure genome indices are properly built
4. Check that FASTQ file naming is consistent

## License

This pipeline is provided as-is for research use.