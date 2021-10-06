#!/bin/bash
### Modified by José Luis Ruiz to be used as part of the ILRA pipeline in 2022


### Check arguments
genome=$1
kmer=$2
stepsize=$3
readRoot=$4
export readRoot
resultname=$5
insertSize=$6
cores=$7
tmp=$$

if [ -z "$insertSize" ] || [ -z "$genome" ] || [ -z "$kmer" ] || [ -z "$stepsize" ] || [ -z "$readRoot" ] || [ -z "$resultname" ] ; then
	echo "Usage: icorn2.bowtie2.sh <Genome> <k-mer> <StepSize> <Illumina reads root name> <ResultsName> <Insertsize Range> <cores>"
	exit;
fi

if [ -z "$cores" ] ; then
	cores=$ICORN2_THREADS
	if [ -z "$ICORN2_THREADS" ] ; then
		cores=1
	fi
fi

samtools faidx $genome
mkdir $resultname


### Executing bowtie2...
bowtie2-build --threads $cores $genome $genome &> output.bowtiebuild.txt
bowtie2 -x $genome -p $ICORN2_THREADS -X $insertSize --very-sensitive -N 1 -L 31 --rdg 5,2 -1 "$readRoot"_1.fastq -2 "$readRoot"_2.fastq | awk -va=$readRoot '{if ($1~/^@/) {print} else {print $0"\tRG:Z:"a}}' | samtools view -@ $cores -f 0x2 -u -t $genome.fai - | samtools reheader -c 'sed "$ a\@RG\tID:$readRoot\tSM:1"' - | samtools sort -@ $cores -u -m 2G -o $resultname/out.bam -
echo -e "\nbowtie2 DONE"

### Executing MarkDuplicates with updated java (make sure not to use the java 1.7 that's needed for the modified GenomeAnalysisTK.jar, but a more recent one)
$ICORN2_HOME/jdk-16.0.2/bin/java -XX:ParallelGCThreads=$cores -XX:-UseParallelGC -XX:+UseSerialGC -jar $ICORN2_HOME/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT M=/dev/null ASSUME_SORTED=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=$resultname/out.bam O=$resultname/out.sorted.markdup.bam
return=$?
if [ "$return" != "0" ] ; then
	echo "Sorry, MarkDuplicates failed..."
	exit 1;
else 
	samtools index -@ $cores $resultname/out.sorted.markdup.bam
	rm $resultname/out.bam; echo -e "\nMarkDuplicates DONE"
fi

