conda activate ngs
# conda install nextflow

# download geo data
cd /media/hao/Data1/NGS/2025-04-16_H358_RNAseq_ATACseq/RNAseq

fastq-dump --gzip --split-3 --outdir rawdata GSM7439318
fastq-dump --gzip --split-3 --outdir rawdata GSM7439319
fastq-dump --gzip --split-3 --outdir rawdata GSM7439320
fastq-dump --gzip --split-3 --outdir rawdata GSM7439321
fastq-dump --gzip --split-3 --outdir rawdata GSM7439322
fastq-dump --gzip --split-3 --outdir rawdata GSM7439323
fastq-dump --gzip --split-3 --outdir rawdata GSM7439324
fastq-dump --gzip --split-3 --outdir rawdata GSM7439325
fastq-dump --gzip --split-3 --outdir rawdata GSM7439326


# run pipeline
nextflow run nf-core/rnaseq \
	-profile docker \
    --input samplesheet.csv \
    --outdir results \
    --pseudo_aligner salmon \
    --skip_alignment \
	--extra_salmon_quant_args '--seqBias --gcBias' \
	--salmon_index /media/hao/Data1/index/hs/salmon \
    --gtf /media/hao/Data1/index/hs/GRCh38.gtf \
    --fasta /media/hao/Data1/index/hs/GRCh38.fa
# -resume



nextflow run nf-core/atacseq \
	-profile docker \
	--input samplesheet.csv \
	--macs_gsize 2913022398 \
	--outdir results \
	--mito_name chrM \
	--aligner bowtie2 \
	--bowtie2_index /media/hao/Data1/index/hs/bowtie2 \
	--gtf /media/hao/Data1/index/hs/GRCh38.gtf \
	--fasta /media/hao/Data1/index/hs/GRCh38.fa
