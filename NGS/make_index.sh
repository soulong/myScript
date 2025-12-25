#! /bin/bash
# ubuntu_24.04_LTS
# last edited at 2025-11-13 by Hao He

#################################  setup  #################################
conda activate ngs
# conda env export --no-builds | grep -v "^prefix: " > ngs.yml
# conda env create -f ngs.yml

# set directory
root_dir=/mnt/d/Index/hs/v49 \
	&& mkdir -p $root_dir \
	&& cd $root_dir \
	&& echo ">>>>>> start to processing <<<<<<<"

# hs
fasta_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
cdna_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
gtf_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz

# # mm
# fasta_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
# cdna_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.transcripts.fa.gz
# gtf_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.basic.annotation.gtf.gz



#################################  process data  #################################
# fasta
fasta_gz=$(basename "$fasta_link")
[ ! -f "$fasta_gz" ] && echo "downloading fasta" && wget $fasta_link -O $fasta_gz
#[ ! -f "${fasta_gz%.gz}.fai" ] && echo "index fasta" && gunzip -k $fasta_gz && samtools faidx ${fasta_gz%.gz}
# cdna
cdna_gz=$(basename "$cdna_link")
[ ! -f "$cdna_gz" ] && echo "downloading cdna" && wget $cdna_link -O $cdna_gz
# gtf
gtf_gz=$(basename "$gtf_link")
[ ! -f "$gtf_gz" ] && echo "downloading gtf" && wget $gtf_link -O $gtf_gz


#################################  salmon index  #################################
# salmon use ultra fast quasi mapping on transcript instead of genome mapping
mkdir -p salmon
# generate decoy
zgrep "^>" $fasta_gz | cut -d " " -f 1 > decoys.txt && sed -i -e 's/>//g' decoys.txt
# add custome genome, usually like FPs, Cas9, resistance genes
if [[ -f "custom_genome.fa" && -f "custom_genome.gtf" ]]; then
	echo "add custome genome ...."
	gffread -w custom_cdna.fa -g custom_genome.fa custom_genome.gtf && gzip custom_cdna.fa
	gzip -k custom_genome.fa && \
		zgrep "^>" custom_genome.fa.gz | cut -d " " -f 1 > custom_decoys.txt && \
		sed -i -e 's/>//g' custom_decoys.txt && cat custom_decoys.txt >> decoys.txt
	# generate gentrome
	zcat custom_cdna.fa.gz $cdna_gz custom_genome.fa.gz $fasta_gz > gentrome.fa.gz
else
	zcat $cdna_gz $fasta_gz > gentrome.fa.gz
fi
# index
salmon index -p 8 -t gentrome.fa.gz -d decoys.txt -i salmon --gencode
# remove unused files
if [ $? -eq 0 ]; then rm custom_decoys.txt custom_genome.fa.gz custom_cdna.fa.gz gentrome.fa.gz custom_genome.fa.fai; fi


#################################  bowtie2 index  #################################
# note that bowtie2 is not a junction aware aligner
mkdir -p bowtie2
bowtie2-build --threads 8 $fasta_gz bowtie2/bowtie2


#################################  star index  #################################
# star is specifically designed for RNA mapping on genome
# mkdir -p star
# STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star \
# 	--genomeFastaFiles $fasta_gz --sjdbGTFfile $gtf_gz \
# 	--sjdbOverhang 149 --limitGenomeGenerateRAM 40000000000 # ~ x/1024/1000/1000 GB
