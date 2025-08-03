#! /bin/bash

thread=28
fastq_dir="01.RawData"
bam=false
salmon_index="/media/hao/Data1/index/hs/salmon"
bowtie2_index="/media/hao/Data1/index/hs/bowtie2/bowtie2"
star_index="/media/hao/Data1/index/hs/star"



echo "======================================================="
echo "================== parameter prasing =================="
echo "======================================================="

# get all gzipped fastq files
fastq_files=($(find ${fastq_dir} -type f -name "*.gz" | sort -u))
# echo "${fastq_files[@]}"
# printf '%s\n' "${fastq_files[@]}"

# get fastq file suffix
suffix_read1=$( echo ${fastq_files[@]} | egrep -o "_([rR])?1.f(ast)?q.gz" | sort -u )
suffix_read2=$( echo ${fastq_files[@]} | egrep -o "_([rR])?2.f(ast)?q.gz" | sort -u )
echo -e "\nextract ${suffix_read1} as read-1 suffix"
echo -e "extract ${suffix_read2} as read-2 suffix\n"

# get fastq file read 1/2
fastq_read1=($( printf '%s\n' "${fastq_files[@]}" | grep ${suffix_read1} | sort -u))
fastq_read2=($( printf '%s\n' "${fastq_files[@]}" | grep ${suffix_read2} | sort -u))
# echo "${fastq_read1[@]}"
# echo "${fastq_read2[@]}"

# get fastq file basenames
basenames=()
for i in ${fastq_read1[@]}
	do
		name=$(basename $i)
		basename=${name%${suffix_read1}}
		# echo $basename
		basenames+=($basename)
	done
echo -e "\n------------- pharse paried samples ------------"
printf '%s\n' "${basenames[@]}"
echo -e "-------------------------------------------------\n"



echo -e "\n\n"
echo "======================================================="
echo "=========== step 1: run fastqc on raw data ============"
echo "======================================================="

fastqc_output_dir="fastqc"

if [ ! -d "${fastqc_output_dir}" ]; then

	mkdir -p ${fastqc_output_dir}
	echo `date` >> ${fastqc_output_dir}/fastqc.log

	SECONDS=0

	fastqc -t ${thread} -o ${fastqc_output_dir} ${fastq_files[@]} \
		2>> ${fastqc_output_dir}/fastqc.log

	duration=${SECONDS}
	echo ">> $((${duration} / 60)) minutes and $((${duration} % 60)) seconds elapsed"

else
  echo "existing fastqc directory found, skip fastqc step"
fi



echo -e "\n\n"
echo "======================================================="
echo "================= step 2: run fastp ==================="
echo "======================================================="

fastq_clean_dir="fastq_clean"
mkdir -p ${fastq_clean_dir} 
echo `date` >> ${fastq_clean_dir}/fastp.log

for i in ${!basenames[@]}
  do
	    if [ -f "${fastq_clean_dir}/${basenames[$i]}_R1.fastq.gz" ]; then
	       echo "${basenames[$i]} exists, jump to next one"
	      continue
	    fi

	    trap 'echo stopped, removing created files; \
	      rm ${fastq_clean_dir}/${basenames[$i]}*; \
	      exit 1' SIGINT SIGTERM

	    echo -e "\n>>>>>>> trim on ${basenames[$i]} ..."

	    SECONDS=0

	    fastp -w $thread \
	    	-i ${fastq_read1[$i]} \
	      -I ${fastq_read2[$i]}\
	      -o ${fastq_clean_dir}/${basenames[$i]}_R1.fastq.gz \
	      -O ${fastq_clean_dir}/${basenames[$i]}_R2.fastq.gz \
	      -j ${fastq_clean_dir}/${basenames[$i]}.json \
	      -h ${fastq_clean_dir}/${basenames[$i]}.html \
	      -R "${basenames[$i]}" \
	      2>> ${fastq_clean_dir}/fastp.log
	          
	    duration=$SECONDS
	    echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
	    sleep 1
  done

trap '' SIGINT SIGTERM



echo -e "\n\n"
echo "=============================================================="
echo "=========== step 3: pseudo alignment with salmon ============="
echo "=============================================================="

# # make salmon index
# salmon index -p 8 -i "/mnt/d/Index/hs/salmon" -t "/mnt/d/Index/hs/cdna.fa"

salmon_dir="salmon"
salmon_input_dir=${fastq_clean_dir}
mkdir -p ${salmon_dir}
echo `date` >> ${salmon_dir}/salmon.log

for i in ${!basenames[@]}
	do
		if [ -d "${salmon_dir}/${basenames[$i]}" ]; then
	        echo "${basenames[$i]} exists, jump to next one"
	        continue
	    fi

	    trap 'echo stopped, removing created files; \
	        rm -r ${salmon_dir}/${basenames[$i]}*; \
	        exit 1' SIGINT SIGTERM

		echo -e "\n>>>>>>> salmon quant on ${basenames[$i]} ..."

		SECONDS=0

		salmon quant -i ${salmon_index} -l A -p $thread \
	        -1 ${salmon_input_dir}/${basenames[$i]}_R1.fastq.gz \
	        -2 ${salmon_input_dir}/${basenames[$i]}_R2.fastq.gz \
	        --seqBias --gcBias --validateMappings \
	        -o ${salmon_dir}/${basenames[$i]} \
			2>> ${salmon_dir}/salmon.log

		duration=$SECONDS
	    echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
	    sleep 1
	done

trap '' SIGINT SIGTERM





echo -e "\n\n"
echo "=============================================================="
echo "============== step 3: alignment with star ================"
echo "=============================================================="

if [ "$bam" == "true" ]; then

	star_dir="bam"
	star_input_dir=${fastq_clean_dir}
	mkdir -p ${star_dir}
	echo `date` >> ${star_dir}/star.log

	for i in ${!basenames[@]}
		do
			if [ -d "${star_dir}/${basenames[$i]}" ]; then
		        echo "${basenames[$i]} exists, jump to next one"
		        continue
		    fi

		    trap 'echo stopped, removing created files; \
		        rm -r ${star_dir}/${basenames[$i]}*; \
		        exit 1' SIGINT SIGTERM

			mkdir -p ${star_dir}/${basenames[$i]}

			echo -e "\n>>>>>>> star align on ${basenames[$i]} ..."

			SECONDS=0

			STAR --runThreadN $thread --genomeDir ${star_index} --readFilesCommand zcat \
			--readFilesIn ${star_input_dir}/${basenames[$i]}_R1.fastq.gz ${star_input_dir}/${basenames[$i]}_R2.fastq.gz \
			--outFileNamePrefix ${star_dir}/${basenames[$i]} \
			--outSAMtype BAM SortedByCoordinate

			# echo "samtools sort ..."
			# samtools view -@ $thread -Sbh ${star_dir}/$basename/Aligned.out.sam | \
			# samtools sort -@ $thread -m 2G -o ${star_dir}/$basename/$basename.bam

			# echo "samtools index ..."
			# samtools index ${star_dir}/$basename/$basename.bam

			# rm ${star_dir}/$basename/*.sam
			
			duration=$SECONDS
		    echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
		    sleep 1
		done

	trap '' SIGINT SIGTERM

fi



# echo -e "\n\n"
# echo "=============================================================="
# echo "============ step 3: featurecounts on gene level ============="
# echo "=============================================================="

# if [ "$bam" == "true" ]; then

# 	featurecounts_dir="featurecounts"
# 	featurecounts_input_dir=${star_dir}
# 	# strandness: IU = 0, ISF = 1, ISR = 2
# 	featurecounts_strandness=0
# 	mkdir -p ${featurecounts_dir}
# 	echo `date` >> ${featurecounts_dir}/featurecounts.log

# 	if [ ! -f ${featurecounts_dir}/gene_featurecounts.txt ]; then
		
# 		SECONDS=0

# 		featureCounts -T $thread -s ${featurecounts_strandness} \
# 			-t exon -g gene_id -Q 20 -p -C --donotsort \
# 			-a ${gtf_file} \
# 			-o ${featurecounts_dir}/gene_featurecounts.txt `ls $featurecounts_input_dir/*/*.bam` \
# 			2>> ${featurecounts_dir}/featurecounts.log

# 		duration=$SECONDS
# 		echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
# 		sleep 1
# 	else
# 		echo "featureCounts result has found, skip"
# 	fi

# fi



echo -e "\n\n"
echo "=============================================================="
echo "================= step 4: multiQC evaluation ================="
echo "=============================================================="

multiqc_dir="multiqc"

if [ ! -d ${multiqc_dir} ]; then

	mkdir -p ${multiqc_dir}

	SECONDS=0

	multiqc -f . -o ${multiqc_dir}

	duration=$SECONDS
	echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
	sleep 1
	
else
	echo "multiqc directory has found, skip"
fi




