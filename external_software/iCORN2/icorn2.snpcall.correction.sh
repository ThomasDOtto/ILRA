#!/bin/bash
### Modified by José Luis Ruiz to be used as part of the ILRA pipeline in 2022


### Check arguments
genome=$1
aln=$2
cores=$3
readRoot=$4
fragmentSize=$5
outname=$6
iternum=$7

if [ ! -z "$ICORN2_VERBOSE" ] ; then
  echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -I out.sorted.markdup.bam -o gatk_variants.realign.intervals -R ref.fa -nt $ICORN2_THREADS"
  echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T IndelRealigner -I out.sorted.markdup.bam -o gatk.realigned.bam -targetIntervals gatk_variants.realign.intervals -R ref.fa"
  echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper -I gatk.realigned.bam -pnrm POOL -o gatk_variants.variants.vcf -ploidy 1 -glm POOLBOTH -R ref.fa -nt $ICORN2_THREADS"
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
echo -e "\nITERATION $iternum - FRAGMENT $genome: DONE\n"
