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
	cores=1
fi
echo -e "Current cores: $cores\n"

if [ $low_spa_mode == "yes" ]; then
	output_folder=/dev/shm
	echo -e "\n You have selected low_spa mode for iCORN2. Accordingly, some steps are going to be performed in /dev/shm, which will increase memory usage or maybe not work in some systems..."
else
	output_folder=$PWD
fi

if [ $low_mem_mode == "yes" ]; then
	cores=$((cores / 2))
	echo -e "\n You have selected low_memory mode for iCORN2. Accordingly, Bowtie2 and MarkDuplicatesSpark will be executed using half the cores to try and limit memory usage. If still running into memory issues, try manually decreasing even more..."
fi

samtools faidx $genome 
mkdir -p $output_folder/$resultname; mkdir -p $resultname

### Executing bowtie2...
echo -e "\nCalling Bowtie2..."
bowtie2-build --threads $cores $genome $output_folder/$genome &> $resultname/bowtiebuild.log_out.txt
bowtie2 -t -x $output_folder/$genome -p $cores -X $insertSize --very-sensitive -N 1 -L 31 --rdg 5,2 -1 "$readRoot"_1.fastq.gz -2 "$readRoot"_2.fastq.gz | awk -va=$readRoot '{if ($1~/^@/) {print} else {print $0"\tRG:Z:"a}}' | samtools view -@ $cores -f 0x2 -u -t $genome.fai - | samtools reheader -c 'sed "$ a\@RG\tID:$readRoot\tSM:1"' - | samtools sort -@ $((cores / 2)) -n -l 9 -m 10G -o $output_folder/$resultname/out.bam -
rm $output_folder/$genome.* # Remove the bowtie2 index
echo -e "\nBowtie2 DONE"

### Executing MarkDuplicates
echo -e "\nCalling MarkDuplicatesSpark (GATK)..."
# Add the argument -M marked_dup_metrics.txt if you want to gather the stats, but by default omitted as it slows down
echo -e "\nThe use of OpenJDK v8 is going to try and be forced, by adding to the PATH once... but custom paths may be required here...\n"
PATH=$ICORN2_HOME/jdk8u302-b08/bin:$PATH
$ICORN2_HOME/gatk MarkDuplicatesSpark -I $output_folder/$resultname/out.bam -O $output_folder/$resultname/out.sorted.markdup.bam -VS SILENT -OBI true --create-output-bam-splitting-index false --tmp-dir $TMPDIR --spark-master local[$((cores / 2))] &> $resultname/MarkDuplicatesSpark_log.out.txt

### Ending:
return_markdup=$?
export return_markdup
if [ "$return_markdup" != "0" ] ; then
	echo -e "\nSorry, MarkDuplicatesSpark failed... Maybe due to excessive memory required? (consider using less cores, setting up Xmx in _JAVA_OPTIONS variable?). Maybe due to excessive open files (see ulimit -a)?. Maybe due to incorrect JAVA version (OpenJDK v8 is strictly required by GATK, and has to be in the PATH, double check the PATH)? Please check the log MarkDuplicatesSpark_log.out.txt to find more information...\n"
	exit 1;
else
	rm $output_folder/$resultname/out.bam
	rm -rf $TMPDIR/*
fi
echo -e "\nMarkDuplicates DONE"
