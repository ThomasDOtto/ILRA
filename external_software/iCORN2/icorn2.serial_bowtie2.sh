#!/bin/bash
### Modified by José Luis Ruiz to improve performance and to be used as part of the ILRA pipeline in 2022. Check out the changelog to get a broad overview of the updates


### Checking and setting arguments by default if not provided
readRoot=$1
fragmentSize=$2
referenceORG=$3
start=$4
end=$5
splitter_parts=$6
low_mem_mode=$7; export low_mem_mode # export to be used also in other scripts

echo -e "We are automatically setting ICORN2_HOME to $(dirname "$0")\n"
ICORN2_HOME=$(dirname "$0"); export ICORN2_HOME

if [ -z "$readRoot" ] || [ -z "$fragmentSize" ] || [ -z "$referenceORG" ] || [ -z "$start" ] || [ -z "$end" ] || [ -z "$splitter_parts" ] || [ ! -f "$referenceORG" ] ; then
	echo "Oh $(whoami), you have to call me and provide Illumina short reads and a correct reference
### Usage:
icorn2.serial_bowtie2.sh <Illumina reads root name> <Illumina reads fragment size> <Reference genome> <Number iteration to start> <Number iteration to stop> <Number of parts to split the sequences> <Low memory mode>

Multithreading would be used in several steps, such as the mapping stage and the GATK variant calls. Please set the environmental variable ICORN2_THREADS, for example for using 8 cores 'ICORN2_THREADS=8; export ICORN2_THREADS
For more debugging information, please set the environmental variable ICORN2_VERBOSE to an integer, for example 'ICORN2_VERBOSE=1; export ICORN2_VERBOSE'

iCORN2 expects Illumina short reads whose root name is provided with absolute paths (the files should be prepared to follow the name convention <Illumina reads root name>_1.fastq.gz and <Illumina reads root name>_2.fastq.gz)
To continue a correction after three iteration just put 4 in the <Number iteration to start>, or one unit more than the iteration that you want to resume
To activate the low memory mode, provide 'yes' as the argument. Otherwise, use 'no' or leave empty ('no' is used as default)
iCORN2 is dividing the input sequences in N parts, which will be by default simultaneously processed. To this end, the number of cores will be automatically distributed, but the memory usage may increase. If the number of parts to split the input sequences is larger than the number of ICORN2_THREADS, or larger than the number of input sequences, or if you provide 'yes' to the low memory argument, the parts will be sequentially processed instead, using half the ICORN2_THREADS. If still running into memory issues, try and decrease the cores even more or reduce the number of parts the input sequences are divided into...

"
	exit 1;

fi

if [ -z "$ICORN2_THREADS" ] ; then
	ICORN2_THREADS=1; export ICORN2_THREADS
fi
echo -e "Current ICORN2_THREADS: $ICORN2_THREADS cores\n"
echo -e "Please change and export the environmental variable ICORN2_THREADS if necessary"

if [ -z "$low_mem_mode" ] ; then
	low_mem_mode="no"; export low_mem_mode # export to be used also in other scripts
fi

if [ "$splitter_parts" -gt $(grep -c ">" $referenceORG) ]; then
	splitter_parts=$(grep -c ">" $referenceORG)
fi

cores_split=$((ICORN2_THREADS / splitter_parts)) # subset of cores for the simultaneous SNP calling threads

if [ $cores_split -eq 0 ]; then
	low_mem_mode="yes"; export low_mem_mode # export to be used also in other scripts
fi

mkdir -p tmp_dir # make local temporal dir to store gatk and picard

### Setting up the reference sequences...
refRoot=$(echo $referenceORG | perl -nle '@ar=split(/\//); print $ar[($ar[scalar(@ar)]-1)]') # For getting a copy of the reference in the folder, without | in the name
if [ ! -f "ICORN2.$refRoot.$start" ] ; then
	cat $referenceORG | perl -e 'while(<STDIN>) {if (/>(\S+)/){ chomp; s/\|/_/g; s/\./_/g; print "$_\n";} else {print}}' > ICORN2.$refRoot.$start
fi


### Decompressing the Illumina short reads...
# This may be time-consuming, but unfortunately it is mandatory, since SNP-o-matic in the iCORN2's correction step does not allow gzipped .fastq.
# Since it's supossed to only use a core, decompression of both sets of pairs is done simultaneously...
(
for i in {1..2}; do
	pigz -dfc -p $ICORN2_THREADS $readRoot\_$i.fastq.gz > $PWD/"${readRoot##*/}"\_$i.fastq &
done
) 2>&1 | cat -u >> processing_decompressing_log_out.txt
readRoot_uncompressed=$PWD/"${readRoot##*/}"


### Executing iCORN2...
for ((i=$start;$i<=$end;i++)); do
	echo -e "\n\n\n#### ITERATION ++++ $i"
	time1=`date +%s`
### Call the mapper
	echo -e "\nCalling the mapper...\n"
	$ICORN2_HOME/icorn2.mapper.sh ICORN2.$refRoot.$i 13 3 $readRoot_uncompressed ICORN2_$i 1200 $ICORN2_THREADS &> icorn2.mapper_log_out.txt
	echo -e "\nMapper DONE\n"
### Call SNP caller and correction in splitted sequences:
	cd ICORN2_$i
	ln -sf ../ICORN2.$refRoot.$i ref.fa
# Separating the sequences and subsetting the alignment (with updated java, make sure not to use the java 1.7 that's required by iCORN2's GenomeAnalysisTK.jar, but a more recent one, for example the 1.8 that GATK requires))
	echo -e "\nSeparating the sequences and batch analyzing them...\n"
	fasta-splitter.pl --n-parts $splitter_parts --nopad --measure seq --line-length 0 --out-dir $PWD ref.fa	
	arr=($(find . -name "*part*.fa" -exec basename {} \;)); length_arr=${#arr[@]}
	count1=1; block=0
	while [ $count1 -le $length_arr ]; do
		count2=1
		(
		while [ $count2 -le $splitter_parts ] && [ $count1 -le $length_arr ]; do
			samtools faidx ${arr[$count1 - 1]} &
			echo "sequence=${arr[$((count1 - 1))]}"
			(( count1++ ))
			(( count2++ ))
			echo "count1=$count1"; echo "count2=$count2"
		done
		(( block++ ))
		echo "block=$block"; echo -e "\n"
		) 2>&1 | cat -u >> processing_split_parts_faidx_log_out.txt
		eval $(grep "count1=" processing_split_parts_faidx_log_out.txt | tail -1)
		eval $(grep "block=" processing_split_parts_faidx_log_out.txt | tail -1)
	done
		
	count1=1; block=0
	while [ $count1 -le $length_arr ]; do
		count2=1
		(
		while [ $count2 -le $splitter_parts ] && [ $count1 -le $length_arr ]; do
			element=${arr[$count1 - 1]}; element2=${element%.*}
			java -XX:-UseParallelGC -XX:ParallelGCThreads=$((ICORN2_THREADS / 2)) -jar $ICORN2_HOME/picard.jar CreateSequenceDictionary R=$element O=$element2.dict TMP_DIR=../tmp_dir &				
			echo "sequence=$element"
			(( count1++ ))
			(( count2++ ))
			echo "count1=$count1"; echo "count2=$count2"
		done
		(( block++ ))
		echo "block=$block"; echo -e "\n"
		) 2>&1 | cat -u >> processing_split_parts_createseqdictionary_log_out.txt
		eval $(grep "count1=" processing_split_parts_createseqdictionary_log_out.txt | tail -1)
		eval $(grep "block=" processing_split_parts_createseqdictionary_log_out.txt | tail -1)
	done
	if [ $low_mem_mode == "no" ]; then
# Executing in parallel each subset of the sequences
		echo -e "Cores to be used simultaneously to process each split sequence: $cores_split"
		echo -e "\nProcessing simultaneously splitted sequences in blocks of at most $splitter_parts elements, please manually change the variable 'iCORN2_blocks_size' in the ILRA.sh main script if required, for example because less cores available or running into memory issues...\n)"		
		count1=1; block=0
		while [ $count1 -le $length_arr ]; do
			count2=1
			(
			while [ $count2 -le $splitter_parts ] && [ $count1 -le $length_arr ]; do
				element=${arr[$count1 - 1]}; element2=${element%.*}
				java -XX:-UseParallelGC -XX:ParallelGCThreads=$cores_split -jar $ICORN2_HOME/picard.jar ReorderSam INPUT=out.sorted.markdup.bam OUTPUT=$element.bam SEQUENCE_DICTIONARY=$element2.dict REFERENCE_SEQUENCE=$element S=true VERBOSITY=WARNING COMPRESSION_LEVEL=1 CREATE_INDEX=true TMP_DIR=../tmp_dir &
				echo "sequence=$element"
				(( count1++ ))
				(( count2++ ))
				echo "count1=$count1"; echo "count2=$count2"
			done
			(( block++ ))
			echo "block=$block"; echo -e "\n"
			) 2>&1 | cat -u >> processing_split_parts_reordersam_log_out.txt
			eval $(grep "count1=" processing_split_parts_reordersam_log_out.txt | tail -1)
			eval $(grep "block=" processing_split_parts_reordersam_log_out.txt | tail -1)
		done

		count1=1; block=0
		while [ $count1 -le $length_arr ]; do
			count2=1
			(
			while [ $count2 -le $splitter_parts ] && [ $count1 -le $length_arr ]; do				
				icorn2.snpcall.correction.sh ${arr[$count1 - 1]} ${arr[$count1 - 1]}.bam $cores_split $readRoot_uncompressed $fragmentSize ${arr[$count1 - 1]}.$(($i+1)) $i &				
				echo "sequence=${arr[$count1 - 1]}"
				(( count1++ ))
				(( count2++ ))
				echo "count1=$count1"; echo "count2=$count2"
			done
			(( block++ ))
			echo "block=$block"; echo -e "\n"
			) 2>&1 | cat -u >> processing_split_parts_correction_icorn2_log_out.txt
			eval $(grep "count1=" processing_split_parts_correction_icorn2_log_out.txt | tail -1)
			eval $(grep "block=" processing_split_parts_correction_icorn2_log_out.txt | tail -1)
		done		
	elif [ $low_mem_mode == "yes" ]; then
# Executing sequentially each subset of the sequences
		echo -e "Cores to be used to sequentially process each split sequence: $((ICORN2_THREADS / 2))"
		for part in $(ls | grep -e ".part.*.fa$"); do
			java -XX:-UseParallelGC -XX:ParallelGCThreads=$((ICORN2_THREADS / 2)) -jar $ICORN2_HOME/picard.jar ReorderSam INPUT=out.sorted.markdup.bam OUTPUT=$part.bam SEQUENCE_DICTIONARY=${part%.*}.dict REFERENCE_SEQUENCE=$part S=true VERBOSITY=WARNING COMPRESSION_LEVEL=1 CREATE_INDEX=true TMP_DIR=../tmp_dir
			icorn2.snpcall.correction.sh $part $part.bam $((ICORN2_THREADS / 2)) $readRoot_uncompressed $fragmentSize $part.$(($i+1)) $i
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
echo -e "\n\n\nTo look in into a correction, open the file ending with .1 in artemis and load the gff file onto it..."
