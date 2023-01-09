#!/bin/bash
### Modified by José Luis Ruiz to be used as part of the ILRA pipeline in 2022


### Check arguments
genome=$1
cores=$2
readRoot=$3
fragmentSize=$4
outname=$5
iternum=$6

### Create the bam file:
java_memory_2=${java_memory::-1} # Excessive use of RAM by the parallel execution of ReorderSam, so try and limit within the user specified range
# 1.25X the size of the blocks in the parallel execution of this script to make sure RAM is not filled and jobs can accumulate
if [ $seq_parts -eq 0 ]; then
	#java_memory_reorder_sam="$((java_memory_2 / (( parallel_block_size * 5/4 )) ))"g
	java_memory_reorder_sam="$((java_memory_2 / (( parallel_block_size )) ))"g
else
	java_memory_reorder_sam="$((java_memory_2 / seq_parts ))"g
fi

if [ ! -z "$length_arr" ] ; then
	java_memory_reorder_sam="$((java_memory_2 / length_arr ))"g
fi

export _JAVA_OPTIONS="-Xms5g -Xmx$java_memory_reorder_sam"
java -XX:-UseParallelGC -XX:ParallelGCThreads=$cores -Djava.io.tmpdir=../tmp_dir -jar $ICORN2_HOME/picard.jar ReorderSam -INPUT out.sorted.markdup.bam -OUTPUT $genome.bam -SEQUENCE_DICTIONARY ${genome%.*}.dict -REFERENCE_SEQUENCE $genome -S true -VERBOSITY WARNING -COMPRESSION_LEVEL 1 -CREATE_INDEX true -TMP_DIR ../tmp_dir

export _JAVA_OPTIONS="-Xms5g -Xmx$java_memory"

aln=$genome.bam
cores=$((cores * 2)) # it's not very intensive and snp-o-matic/Pilon mostly use only one core, so it can accumulate...

genome_size=$(stat --printf="%s" $genome)

if (( genome_size < maxsize)); then
	echo -e "\n\nCorrecting $genome with iCORN2...\n\n"
	if [ ! -z "$ICORN2_VERBOSE" ] ; then
	  echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -I out.sorted.markdup.bam -o gatk_variants.realign.intervals -R ref.fa -nt $cores"
	  echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T IndelRealigner -I out.sorted.markdup.bam -o gatk.realigned.bam -targetIntervals gatk_variants.realign.intervals -R ref.fa"
	  echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper -I gatk.realigned.bam -pnrm POOL -o gatk_variants.variants.vcf -ploidy 1 -glm POOLBOTH -R ref.fa -nt $cores"
	  echo "To achieve performance improvements, an approach based on dividing the sequences and computing simultaneously has been implemented"
	fi

	mkdir -p $genome.out
	### Call the SNP caller
	echo -e "\nCalling the SNP caller... $genome\n"
	$ICORN2_HOME/jdk1.7.0_80/bin/java -XX:-UseParallelGC -XX:ParallelGCThreads=$cores -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $aln -o $genome.out/gatk_variants.realign.${genome%.*}.intervals -R $genome -nt $cores &> $genome.out/gatk.1.${genome%.*}.log_out.txt
	$ICORN2_HOME/jdk1.7.0_80/bin/java -XX:-UseParallelGC -XX:ParallelGCThreads=$cores -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T IndelRealigner -I $aln -o $genome.out/gatk.realigned.${genome%.*}.bam -targetIntervals $genome.out/gatk_variants.realign.${genome%.*}.intervals -R $genome &> $genome.out/gatk.2.${genome%.*}.log_out.txt
	$ICORN2_HOME/jdk1.7.0_80/bin/java -XX:-UseParallelGC -XX:ParallelGCThreads=$cores -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper -I $genome.out/gatk.realigned.${genome%.*}.bam -pnrm POOL -o $genome.out/gatk_variants.variants.${genome%.*}.vcf -ploidy 1 -glm POOLBOTH -R $genome -nt $cores &> $genome.out/gatk.3.${genome%.*}.log_out.txt
	rm $genome.out/gatk.realigned.${genome%.*}.bam $genome.out/gatk.realigned.${genome%.*}.bai $aln ${aln%.*}.bai
	echo -e "\nSNP caller DONE... $genome\n"

	### Call iCORN2 starter
	echo -e "\nCalling icorn2.correct.sh... $genome\n"
	$ICORN2_HOME/icorn2.correct.sh $genome $genome.out $readRoot $fragmentSize $outname &> ../out.ICORN2.CORRECT.${genome%.*}.$iternum.o
	echo -e "\nicorn2.correct.sh DONE... $genome\n"
	return=$?
	if [ "$return" != "0" ] ; then
	  echo "See errors of the iCORN2 correction script $ICORN2_HOME/icorn2.correct.sh"
	  exit 1;
	fi
elif (( genome_size > maxsize)); then
	echo -e "\n\nCorrecting $genome with Pilon...\n\n"
	if [ ! -z "$length_arr" ] ; then
		java_memory_Pilon="$((java_memory_2 / length_arr ))"g
		export _JAVA_OPTIONS="-Xms5g -Xmx$java_memory_Pilon"		
	fi
	pilon --genome $genome --frags $aln --output $genome"_pilon.fa.$(( iternum + 1)).tmp" --fix bases --outdir $PWD &>> $genome"_pilon_log_out.txt"
	# pilon --genome $genome --frags $aln --output $genome"_pilon.fa.$(( iternum + 1)).tmp" --fix bases --changes --outdir $PWD &>> $genome"_pilon_log_out.txt" # To obtain a list of the changes
	# mv $genome"_pilon.fa.$(( iternum + 1)).tmp.changes" $genome"_pilon.changes"
fi
export _JAVA_OPTIONS="-Xms5g -Xmx$java_memory"

echo -e "\nITERATION $iternum - FRAGMENT $genome: DONE"