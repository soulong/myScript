## under wsl2 ubuntu24.04-LTS

conda activate ngs
# install nextflow under ngs
# conda install nextflow

# enter analysis directory
cd /mnt/c/Users/haohe/Desktop/CUTTAG/

# run pipeline
nextflow run nf-core/cutandrun \
	-profile docker \
	--input samplesheet_CutTag.csv \
	--outdir ./results \
	--skip_trimming true \
	--remove_mitochondrial_reads true \
	--mito_name chrM \
	--normalisation_mode CPM \
	--normalisation_binsize 1 \
	--peakcaller seacr \
	--seacr_norm norm \
	--seacr_stringent relaxed \
	--macs_gsize 2700000000 \
	--macs2_narrow_peak true \
	--extend_fragments true \
	--use_control false \
	--consensus_peak_mode group \
	--replicate_threshold 1 \
	--bowtie2 /mnt/d/Index/hs/bowtie2 \
	--gtf /mnt/d/Index/hs/GRCh38.gtf \
	--fasta /mnt/d/Index/hs/GRCh38.fa \
	--max_cpus 8 \
	--max_memory '42GB'


# chech session id
nextflow -log
# resume a pre-run pipeline or specific session
nextflow run nf-core/cutandrun -resume [session_ID]

	--skip_trimming true \



