#!/bin/bash
### Modified by José Luis Ruiz to be used as part of the ILRA pipeline in 2022


### Check arguments
genome=$1
kmer=$2
stepsize=$3
readRoot=$4; export readRoot # Export to later add the read group flag to the sam header
resultname=$5
insertSize=$6
cores=$7
tmp=$$

echo -e "\nCurrent ICORN2_HOME directory is $ICORN2_HOME\n"

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
echo -e "Current ICORN2_THREADS: $cores cores\n"

samtools faidx $genome; mkdir -p $resultname

### Executing bowtie2...
echo -e "\nCalling bowtie2..."
bowtie2-build --threads $cores $genome $genome &> $resultname/output.bowtiebuild.txt
bowtie2 -t -x $genome -p $cores -X $insertSize --very-sensitive -N 1 -L 31 --rdg 5,2 -1 "$readRoot"_1.fastq -2 "$readRoot"_2.fastq | awk -va=$readRoot '{if ($1~/^@/) {print} else {print $0"\tRG:Z:"a}}' | samtools view -@ $cores -f 0x2 -u -t $genome.fai - | samtools reheader -c 'sed "$ a\@RG\tID:$readRoot\tSM:1"' - | samtools sort -@ $cores -n -u -m 2G -o $resultname/out.bam -
rm $genome.* # Remove the old bowtie2 index, not required anymore
echo -e "\nbowtie2 DONE"

### Executing MarkDuplicates
echo -e "\nCalling MarkDuplicates (GATK)..."
$ICORN2_HOME/gatk MarkDuplicatesSpark -I $resultname/out.bam -O $resultname/out.sorted.markdup.bam -VS SILENT -OBI true --create-output-bam-splitting-index false --spark-master local[$cores] &> $resultname/MarkDuplicatesSpark_log.out.txt
# Add the argument -M marked_dup_metrics.txt if you want to gather the stats, but by default omitted as it slows down
return=$?
if [ "$return" != "0" ] ; then
	echo "Sorry, MarkDuplicatesSpark failed... Maybe due to memory? Maybe due to incorrect java version? (OpenJDK v8 required by GATK). Please check the log..."
	exit 1;
else
	rm $resultname/out.bam; echo -e "\nMarkDuplicates DONE"
fi