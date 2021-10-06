#!/bin/bash
### Modified by José Luis Ruiz to be used as part of the ILRA pipeline in 2022


### Check and set arguments, variables...
dir=$1
readRoot=$2
fragmentSize=$3
resultname=$4
tmp=$$
SNPOMATIC_HOME=$ICORN2_HOME; export SNPOMATIC_HOME    
PERL5LIB=$ICORN2_HOME:$PERL5LIB; export PERL5LIB 


### Executing correction
cd $dir; echo "Calling icorn2.Correct.pl"
echo "perl $ICORN2_HOME/icorn2.Correct.pl ref.fa gatk_variants.variants.vcf ../$readRoot $fragmentSize $resultname"
perl $ICORN2_HOME/icorn2.Correct.pl ref.fa gatk_variants.variants.vcf $readRoot $fragmentSize $resultname

### Get it into correct sequence length
perl $ICORN2_HOME/fasta2singleLine.pl $resultname $resultname.$tmp
perl $ICORN2_HOME/fasta2multipleLine.pl $resultname.$tmp $resultname 80

cp $resultname ../; cd ..
