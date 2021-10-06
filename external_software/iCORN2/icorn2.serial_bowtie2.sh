#!/bin/bash
### Modified by José Luis Ruiz to be used as part of the ILRA pipeline in 2022


### Check and set arguments by default if not provided
readRoot=$1
fragmentSize=$2
referenceORG=$3
start=$4
end=$5

if [ -z "$ICORN2_HOME" ] ; then
	echo "We are automatically setting ICORN2_HOME to $(dirname "$0")"
	ICORN2_HOME=$(dirname "$0"); export ICORN2_HOME
fi
if [ -z "$end" ] || [ -z "$readRoot" ] || [ ! -f "$referenceORG" ] ; then
	echo "Oh $(whoami), you have to call me and provide Illumina short reads and a correct reference
### Usage:
icorn2.serial_bowtie2.sh <Illumina reads root name> <Illumina reads fragment size> <Reference genome> <Number iteration to start> <Number iteration to stop>

iCORN2 expects Illumina short reads whose root name is provided with absolute paths (the files should be prepared to follow the name convention <Illumina reads root name>_1.fastq and <Illumina reads root name>_2.fastq)
To continue a correction after three iteration just put 4 in the <Number iteration to start>, or one unit more than the iteration that you want to resume

### Further parameters via environmental variables:
Multithreading would be used in several steps, such as the mapping stage and the GATK variant calls. Please set the variable ICORN2_THREADS, for example for using 8 cores 'ICORN2_THREADS=8; export ICORN2_THREADS
For debugging information, please set the variable ICORN2_VERBOSE to an integer, for example 'ICORN2_VERBOSE=1; export ICORN2_VERBOSE'
"
	exit 1;

fi

if [ -z "$ICORN2_THREADS" ] ; then
	ICORN2_THREADS=1; export ICORN2_THREADS
fi

refRoot=$(echo $referenceORG | perl -nle '@ar=split(/\//); print $ar[($ar[scalar(@ar)]-1)]') # For getting a copy of the reference in the folder, without | in the name
if [ ! -f "ICORN2.$refRoot.$start" ] ; then
	cat $referenceORG | perl -e 'while(<STDIN>) {if (/>(\S+)/){ chomp; s/\|/_/g; s/\./_/g; print "$_\n";} else {print}}' > ICORN2.$refRoot.$start
fi
reference=ICORN2.$refRoot


### Executing iCORN2...
for ((i=$start;$i<=$end;i++)); do
	directory=ICORN2_$i
	echo -e "\n\n\nIteration ++++ $i"

	### Call the mapper
	$ICORN2_HOME/icorn2.bowtie2.sh ICORN2.$refRoot.$i 13 3 $readRoot ICORN2_$i 1200 $ICORN2_THREADS

	### Call the SNP caller
	cd ICORN2_$i
	ln -s ../$reference.$i ref.fa; samtools faidx ref.fa
	if [ ! -z "$ICORN2_VERBOSE" ] ; then
		echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -I out.sorted.markdup.bam -o gatk_variants.realign.intervals -R ref.fa -nt $ICORN2_THREADS"
		echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T IndelRealigner -I out.sorted.markdup.bam -o gatk.realigned.bam -targetIntervals gatk_variants.realign.intervals -R ref.fa"
		echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper -I gatk.realigned.bam -pnrm POOL -o gatk_variants.variants.vcf -ploidy 1 -glm POOLBOTH -R ref.fa -nt $ICORN2_THREADS"
	fi
	$ICORN2_HOME/jdk-16.0.2/bin/java -XX:ParallelGCThreads=$ICORN2_THREADS -jar $ICORN2_HOME/picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict
	$ICORN2_HOME/jdk1.7.0_80/bin/java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -I out.sorted.markdup.bam -o gatk_variants.realign.intervals -R ref.fa -nt $ICORN2_THREADS &> out.gatk.1.txt
	$ICORN2_HOME/jdk1.7.0_80/bin/java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T IndelRealigner -I out.sorted.markdup.bam -o gatk.realigned.bam -targetIntervals gatk_variants.realign.intervals -R  ref.fa &> out.gatk.2.txt
	$ICORN2_HOME/jdk1.7.0_80/bin/java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper -I gatk.realigned.bam -pnrm POOL -o gatk_variants.variants.vcf -ploidy 1 -glm POOLBOTH -R ref.fa -nt $ICORN2_THREADS &> out.gatk.3.txt
	echo -e "\nSNP caller DONE"; cd ..

	### Call iCORN2 starter
	echo "Calling iCORN2 starter"
	$ICORN2_HOME/icorn2.correct.sh $directory $readRoot $fragmentSize $reference.$(($i+1)) > out.ICORN2.CORRECT.$i.o
	echo "icorn2.correct.sh finished"
	return=$?
	if [ "$return" != "0" ] ; then
		echo "See errors in the iCORN2 correction script: $ICORN2_HOME/icorn2.correct.sh $directory $readRoot $fragmentSize $reference.$(($i+1))"
		exit 1;
	fi
	echo -e "\nIteration $i DONE"
	echo "Current corrections:"; perl $ICORN2_HOME/icorn2.collectResults.pl .
done

echo -e "\n\n\nTo look in into a correction, open the file ending with .1 in artemis and load the gff file onto it..."
