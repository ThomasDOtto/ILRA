#!/bin/bash
### Modified by José Luis Ruiz to improve performance and to be used as part of the ILRA pipeline in 2022. Check out the changelog to get a broad overview of the updates


### Checking and setting arguments by default if not provided
readRoot=$1
fragmentSize=$2
referenceORG=$3
start=$4
end=$5
parallel_block_size=$6; export parallel_block_size
cores=$7
java_memory=$8; export java_memory # export to be used also in other scripts
seq_parts=$9 # iCORN2 may be too slow if there are hundreds or thousands of contigs in the sequences to correct, because by default it divides the sequences and process the contigs in parallel and ReorderSam maps and gets the bam file for each one of them. When a number is provided in this parameter, iCORN2 divides the sequences in "seq_parts" parts, and then it process them in parallel
low_mem_mode=${10}
export _JAVA_OPTIONS="-Xms10g -Xmx$java_memory"

echo -e "We are automatically setting ICORN2_HOME to:"
ICORN2_HOME=$(dirname "$0"); export ICORN2_HOME
echo $ICORN2_HOME

if [ -z "$readRoot" ] || [ -z "$fragmentSize" ] || [ -z "$referenceORG" ] || [ -z "$start" ] || [ -z "$end" ] || [ -z "$parallel_block_size" ] || [ ! -f "$referenceORG" ] ; then
	echo "Oh $(whoami), you have to call me and provide Illumina short reads and a correct reference
### Usage:
icorn2.serial_bowtie2.sh <Illumina reads root name> <Illumina reads fragment size> <Reference genome> <Number iteration to start> <Number iteration to stop> <Number of parts to process simultaneously> <Cores> <Low memory mode>

Multithreading would be used in several steps, such as the mapping stage and the GATK variant calls. For more debugging information, please set the environmental variable ICORN2_VERBOSE to an integer, for example 'ICORN2_VERBOSE=1; export ICORN2_VERBOSE'

iCORN2 expects Illumina short reads whose root name is provided with absolute paths (the files should be prepared to follow the name convention <Illumina reads root name>_1.fastq.gz and <Illumina reads root name>_2.fastq.gz)
To continue a correction use <Number iteration to start> with one unit more than the iteration that you want to resume, for example to continue after three iterations just put 4
To activate the low memory mode, provide 'yes' as the argument. Otherwise, use 'no' or leave empty ('no' is used as default)
iCORN2 is dividing the input sequences in the constituent parts, which will be by default simultaneously processed. To this end, the number of cores will be automatically distributed, but the memory usage may increase. If the number of parts to process simultaneoursly is larger than the number of cores, or larger than the number of input sequences, or if you provide 'yes' to the low memory argument, the parts will be sequentially processed instead, using half the cores. If still running into memory issues, try and decrease the cores even more or reduce the number of parts that are simultaneously processed...

"
	exit 1;

fi

if [ -z "$cores" ] ; then
	cores=1
fi
echo -e "Current cores: $cores cores\n"

if [ -z "$low_mem_mode" ] ; then
	low_mem_mode="no"
fi

if [ -z "$java_memory" ]; then
	java_memory=100g
	echo "Number of threads: "$cores
	export _JAVA_OPTIONS="-Xms10g -Xmx$java_memory"
fi

if [ -z "$seq_parts" ]; then
	seq_parts=0	
fi
export seq_parts

if [ "$parallel_block_size" -gt $(grep -c ">" $referenceORG) ]; then
	parallel_block_size=$(grep -c ">" $referenceORG)
fi

cores_split=$((cores / parallel_block_size)) # subset of cores for the simultaneous SNP calling threads

if [ $cores_split -eq 0 ]; then
	low_mem_mode="yes"
fi

export low_mem_mode
mkdir -p tmp_dir # make local temporal dir to store gatk and picard

### Setting up the reference sequences...
refRoot=$(echo $referenceORG | perl -nle '@ar=split(/\//); print $ar[($ar[scalar(@ar)]-1)]') # For getting a copy of the reference in the folder, without | in the name
if [ ! -f "ICORN2.$refRoot.$start" ] ; then
	cat $referenceORG | perl -e 'while(<STDIN>) {if (/>(\S+)/){ chomp; s/\|/_/g; s/\./_/g; print "$_\n";} else {print}}' > ICORN2.$refRoot.$start
fi


### Decompressing the Illumina short reads...
# This may be time-consuming, but unfortunately it is mandatory, since SNP-o-matic in the iCORN2's correction step does not allow gzipped .fastq.
# Since it's supossed to only use a core, decompression of both sets of pairs can be done simultaneously...
parallel --verbose -j 2 "pigz -dfc -p $cores $readRoot\_{}.fastq.gz > $PWD/"${readRoot##*/}"\_{}.fastq" ::: {1..2} &> processing_decompressing_log_out.txt
readRoot_uncompressed=$PWD/"${readRoot##*/}"


### Executing iCORN2...
for ((i=$start;$i<=$end;i++)); do
	echo -e "\n\n\n#### ITERATION ++++ $i"
	time1=`date +%s`
### Call the mapper
	echo -e "\nCalling the mapper...\n"
	$ICORN2_HOME/icorn2.mapper.sh ICORN2.$refRoot.$i 13 3 $readRoot_uncompressed ICORN2_$i 1200 $cores &>> icorn2.mapper_log_out.txt
	echo -e "\nMapper DONE\n"
### Call SNP caller and correction in splitted sequences:
	cd ICORN2_$i
	ln -sf ../ICORN2.$refRoot.$i ref.fa
# Separating the sequences and subsetting the alignment (with updated java, make sure not to use the java 1.7 that's required by iCORN2's GenomeAnalysisTK.jar, but a more recent one, for example the 1.8 that GATK requires))
	echo -e "\nSeparating the sequences and analyzing them with parallel...\n"
	echo "Make sure not to use the java 1.7 that's required by iCORN2's GenomeAnalysisTK.jar, but a more recent one, for example the 1.8 that GATK requires..."
	# Separating the sequences and creating and array
	# iCORN2 may be too slow if there are hundreds or thousands of contigs in the sequences to correct, because by default it divides the sequences and process the contigs in parallel and ReorderSam maps and gets the bam file for each one of them.
	# If the option is provided by the user, iCORN2 divides the sequences in "seq_parts" parts, and then it process them in parallel. Even if the fasta files being processed are then not a single contigs, overall is way faster if there were hundreds or thousands of contigs in the original sequence
	# iCORN2 is likely going to fail if any of the sequences being processed is larger than or around 250-300Mb. This will try to avoid it, but it'll never cut away individual contigs or chromosomes...
	if [ $seq_parts -eq 0 ]; then
		cat ref.fa | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) "_splitter-part.fa")} print $0 > filename }'		
	else
		fasta-splitter --n-parts $seq_parts --nopad --measure seq --line-length 0 --out-dir $PWD --part-num-prefix .splitter-part- ref.fa
	fi
	maxsize=260000000
	for f in $(ls -l | grep "splitter-part" | sed 's,  *, ,g' | cut -d " " -f 5); do
		if (( f > maxsize)); then echo $f >> larger_contigs.txt; fi
	done
	if [ -s larger_contigs.txt ]; then rm $(ls | grep "splitter-part") && fasta-splitter --part-size $maxsize --nopad --measure seq --line-length 0 --out-dir $PWD --part-num-prefix .splitter-part- ref.fa; fi
	arr=($(find . -name "*splitter-part*" -exec basename {} \; | grep -v ref.fa)); length_arr=${#arr[@]}
	echo -e "\nPlease note iCORN2 is likely going to fail if any of the sequences being processed is larger than or around 300Mb...Keep this in mind if there are segmentation errors later on in SNP-o-matic\n"
	
	# Parallel processing of the sequences...
	parallel --verbose -j $cores samtools faidx {} ::: ${arr[@]} &> processing_split_parts_faidx_log_out.txt
	# Manual parallelization if GNU's parallel is not available... 
	#count1=1; block=0
	#while [ $count1 -le $length_arr ]; do
	#	count2=1
	#	(
	#	while [ $count2 -le $parallel_block_size ] && [ $count1 -le $length_arr ]; do
	#		samtools faidx ${arr[$count1 - 1]} &
	#		echo "sequence=${arr[$((count1 - 1))]}"
	#		(( count1++ ))
	#		(( count2++ ))
	#		echo "count1=$count1"; echo "count2=$count2"
	#	done
	#	(( block++ ))
	#	echo "block=$block"; echo -e "\n"
	#	) 2>&1 | cat -u >> processing_split_parts_faidx_log_out.txt
	#	eval $(grep "count1=" processing_split_parts_faidx_log_out.txt | tail -1)
	#	eval $(grep "block=" processing_split_parts_faidx_log_out.txt | tail -1)
	#done

	parallel --verbose -j $cores java -XX:-UseParallelGC -XX:ParallelGCThreads=$cores -jar $ICORN2_HOME/picard.jar CreateSequenceDictionary -R {}.fa -O {}.dict -TMP_DIR ../tmp_dir ::: $(echo ${arr[@]} | sed 's,.fa,,g') &> processing_split_parts_createseqdictionary_log_out.txt
	# Manual parallelization if GNU's parallel is not available...
	#count1=1; block=0
	#while [ $count1 -le $length_arr ]; do
	#	count2=1
	#	(
	#	while [ $count2 -le $parallel_block_size ] && [ $count1 -le $length_arr ]; do
	#		element=${arr[$count1 - 1]}; element2=${element%.*}
	#		java -XX:-UseParallelGC -XX:ParallelGCThreads=$((cores / 2)) -jar $ICORN2_HOME/picard.jar CreateSequenceDictionary R=$element O=$element2.dict TMP_DIR=../tmp_dir &
	#		echo "sequence=$element"
	#		(( count1++ ))
	#		(( count2++ ))
	#		echo "count1=$count1"; echo "count2=$count2"
	#	done
	#	(( block++ ))
	#	echo "block=$block"; echo -e "\n"
	#	) 2>&1 | cat -u >> processing_split_parts_createseqdictionary_log_out.txt
	#	eval $(grep "count1=" processing_split_parts_createseqdictionary_log_out.txt | tail -1)
	#	eval $(grep "block=" processing_split_parts_createseqdictionary_log_out.txt | tail -1)
	#done

	if [ $low_mem_mode == "no" ]; then
# Executing in parallel each subset of the sequences
		echo -e "Cores to be used simultaneously to process each split sequence: $cores_split"
		echo -e "\nProcessing simultaneously splitted sequences in blocks of at most $parallel_block_size elements, please change the variable 'blocks_size (-b)' in the ILRA.sh main script if required, for example because less cores available or running into memory issues...\n"
		echo -e "\nPlease be aware that temp files adding up to hundreds of GBs can be created during the parallel execution of iCORN2...\n"

		#java_memory_2=${java_memory::-1} # Excessive use of RAM by the parallel execution of ReorderSam, so try and limit within the user specified range
		#java_memory_reorder_sam="$((java_memory_2 / parallel_block_size))"g
		#export _JAVA_OPTIONS="-Xms10g -Xmx$java_memory_reorder_sam"

		#parallel --verbose -j $parallel_block_size java -XX:-UseParallelGC -XX:ParallelGCThreads=$cores_split -jar $ICORN2_HOME/picard.jar ReorderSam -INPUT out.sorted.markdup.bam -OUTPUT {}.fa.bam -SEQUENCE_DICTIONARY {}.dict -REFERENCE_SEQUENCE {}.fa -S true -VERBOSITY WARNING -COMPRESSION_LEVEL 1 -CREATE_INDEX true -TMP_DIR ../tmp_dir ::: $(echo ${arr[@]} | sed 's,.fa,,g') &> processing_split_parts_reordersam_log_out.txt
		#export _JAVA_OPTIONS="-Xms10g -Xmx$java_memory"
		# Manual parallelization if GNU's parallel is not available...
		#count1=1; block=0
		#while [ $count1 -le $length_arr ]; do
		#	count2=1
		#	(
		#	while [ $count2 -le $parallel_block_size ] && [ $count1 -le $length_arr ]; do
		#		element=${arr[$count1 - 1]}; element2=${element%.*}
		#		java -XX:-UseParallelGC -XX:ParallelGCThreads=$cores_split -jar $ICORN2_HOME/picard.jar ReorderSam INPUT=out.sorted.markdup.bam OUTPUT=$element.bam SEQUENCE_DICTIONARY=$element2.dict REFERENCE_SEQUENCE=$element S=true VERBOSITY=WARNING COMPRESSION_LEVEL=1 CREATE_INDEX=true TMP_DIR=../tmp_dir &
		#		echo "sequence=$element"
		#		(( count1++ ))
		#		(( count2++ ))
		#		echo "count1=$count1"; echo "count2=$count2"
		#	done
		#	(( block++ ))
		#	echo "block=$block"; echo -e "\n"
		#	) 2>&1 | cat -u >> processing_split_parts_reordersam_log_out.txt
		#	eval $(grep "count1=" processing_split_parts_reordersam_log_out.txt | tail -1)
		#	eval $(grep "block=" processing_split_parts_reordersam_log_out.txt | tail -1)
		#done
		
		# parallel --verbose -j $((parallel_block_size * 3)) icorn2.snpcall.correction.sh {} {}.bam $cores_split $readRoot_uncompressed $fragmentSize {}.$(($i+1)) $i ::: ${arr[@]} &> processing_split_parts_correction_icorn2_log_out.txt
		# 3X the size of the blocks because it's not very intensive and snp-o-matic only uses one core, so it can accumulate
		# Manual parallelization if GNU's parallel is not available...
		#count1=1; block=0
		#while [ $count1 -le $length_arr ]; do
		#	count2=1
		#	(
		#	while [ $count2 -le $parallel_block_size ] && [ $count1 -le $length_arr ]; do
		#		icorn2.snpcall.correction.sh ${arr[$count1 - 1]} ${arr[$count1 - 1]}.bam $cores_split $readRoot_uncompressed $fragmentSize ${arr[$count1 - 1]}.$(($i+1)) $i &
		#		echo "sequence=${arr[$count1 - 1]}"
		#		(( count1++ ))
		#		(( count2++ ))
		#		echo "count1=$count1"; echo "count2=$count2"
		#	done
		#	(( block++ ))
		#	echo "block=$block"; echo -e "\n"
		#	) 2>&1 | cat -u >> processing_split_parts_correction_icorn2_log_out.txt
		#	eval $(grep "count1=" processing_split_parts_correction_icorn2_log_out.txt | tail -1)
		#	eval $(grep "block=" processing_split_parts_correction_icorn2_log_out.txt | tail -1)
		#done
		# This block was working with a previous version, which performed all the ReorderSam before the correction. However, this accumulated a lot of bam files, requesting too much space somethimes. So I incorporated the ReorderSam execution within the icorn2.snpcall.correction.sh
		
		parallel --verbose -j $((parallel_block_size * 2)) icorn2.snpcall.correction.sh {} $cores_split $readRoot_uncompressed $fragmentSize {}.$(($i+1)) $i ::: ${arr[@]} &> processing_split_parts_correction_icorn2_log_out.txt
		# 2X the size of the blocks because it's not very intensive and snp-o-matic only uses one core, so it can accumulate
		#
	elif [ $low_mem_mode == "yes" ]; then
# Executing sequentially each subset of the sequences
		echo -e "Cores to be used to sequentially process each split sequence: $((cores / 2))"
		for part in $(find . -name "*.fa" -exec basename {} \; | grep -v ref.fa); do
			java -XX:-UseParallelGC -XX:ParallelGCThreads=$((cores / 2)) -jar $ICORN2_HOME/picard.jar ReorderSam INPUT=out.sorted.markdup.bam OUTPUT=$part.bam SEQUENCE_DICTIONARY=${part%.*}.dict REFERENCE_SEQUENCE=$part S=true VERBOSITY=WARNING COMPRESSION_LEVEL=1 CREATE_INDEX=true TMP_DIR=../tmp_dir
			icorn2.snpcall.correction.sh $part $part.bam $((cores / 2)) $readRoot_uncompressed $fragmentSize $part.$(($i+1)) $i
		done
	fi
	cd ../
# Merge the different splits:
	cat $(find $PWD/ICORN2_$i -name "*.fa.$(($i+1))") > ICORN2.$refRoot.$(($i+1)); rm $(find $PWD/ICORN2_$i -name "*.fa.$(($i+1))")
	echo "SNP" > ICORN2_$i/iter.$i.General.stats; awk '/^SNP$/,/^INS$/' $(find $PWD/ICORN2_$i -name "*General.stats" | grep -v iter) | sed "/SNP/d" | sed "/INS/d" | sort >> ICORN2_$i/iter.$i.General.stats
	echo "INS" >> ICORN2_$i/iter.$i.General.stats; awk '/^INS$/,/^DEL$/' $(find $PWD/ICORN2_$i -name "*General.stats" | grep -v iter) | sed "/INS/d" | sed "/DEL/d" | sort >> ICORN2_$i/iter.$i.General.stats
	echo "DEL" >> ICORN2_$i/iter.$i.General.stats; awk '/^DEL$/,/^HETERO$/' $(find $PWD/ICORN2_$i -name "*General.stats" | grep -v iter) | sed "/DEL/d" | sed "/HETERO/d" | sort >> ICORN2_$i/iter.$i.General.stats
	echo "HETERO" >> ICORN2_$i/iter.$i.General.stats; awk '/^HETERO$/,/^Rej.SNP$/' $(find $PWD/ICORN2_$i -name "*General.stats" | grep -v iter) | sed "/HETERO/d" | sed "/Rej.SNP/d" | sort >> ICORN2_$i/iter.$i.General.stats
	echo "Rej.SNP" >> ICORN2_$i/iter.$i.General.stats; awk '/^Rej.SNP$/,/^Rej.INS$/' $(find $PWD/ICORN2_$i -name "*General.stats" | grep -v iter) | sed "/Rej.SNP/d" | sed "/Rej.INS/d" | sort >> ICORN2_$i/iter.$i.General.stats
	echo "Rej.INS" >> ICORN2_$i/iter.$i.General.stats; awk '/^Rej.INS$/,/^Rej.DEL$/' $(find $PWD/ICORN2_$i -name "*General.stats" | grep -v iter) | sed "/Rej.INS/d" | sed "/Rej.DEL/d" | sort >> ICORN2_$i/iter.$i.General.stats
	echo "Rej.DEL" >> ICORN2_$i/iter.$i.General.stats; awk -v RS='(^|\n)Rej.DEL\n' 'END{printf "%s", $0}' $(find $PWD/ICORN2_$i -name "*General.stats" | grep -v iter) | sort >> ICORN2_$i/iter.$i.General.stats
# Collecting results:
	echo "Current corrections:"; perl $ICORN2_HOME/icorn2.collectResults.pl .
	time2=`date +%s`
	echo -e "\nIteration $i DONE"
	echo -e "Elapsed time (secs): $((time2-time1))"; echo -e "Elapsed time (hours): $(echo "scale=2; $((time2-time1))/3600" | bc -l)"
done
rm "$readRoot_uncompressed"_1.fastq "$readRoot_uncompressed"_2.fastq; rm -rf tmp_dir

echo -e "\n\n\nTo look in into a correction, open the file ending with .1, .2, .3... in artemis and load the gff file onto it..."
