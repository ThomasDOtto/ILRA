#!/bin/bash
# Copyright (c) 2011-2015 Genome Research Ltd.
# Author: Thomas D. Otto <tdo@sanger.ac.uk>
#
# This file is part of ABACAS2.
#
# ABACAS2 is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

### Modified by José Luis Ruiz to be used as part of the ILRA pipeline in 2022

reference=$1
contig=$2
cores=$3
ABA_MIN_LENGTH=$4
ABA_MIN_IDENTITY=$5


### Make sure that perl finds the libraries
PERL5LIB=$(dirname "$0"):$PERL5LIB; export PERL5LIB
echo -e "We are adding to PERL5LIB the directory $(dirname "$0")..."


### Help
if [ -z "$contig" ] || [ -z "$reference" ] || [ -z "$ABA_MIN_LENGTH" ] ; then
	echo "

*** Abacas II ***       For any disturbation with this program, please don't blame Sammy!

### Usage
abacas2.nonparallel.sh <Reference> <Contigs to order> <Cores to nucmer>. Optional: <min aligment length> <Identity cutoff>

reference:           Fasta (or multi-fasta) against which the contigs should be ordered
contig:              Contigs or query that should be ordered
Min aligment length: Alignment length significance threshold (default 200)
Identity cutoff:     Threshold for identity to place contigs (default 95)

### Further parameters via envinromental variables:
ABA_CHECK_OVERLAP=1; export ABA_CHECK_OVERLAP # This will try to overlap contigs
ABA_splitContigs=1; export ABA_splitContigs # This will split contigs. This is good to split the orign, and to find rearrangement. A split contigs has the suffix _i (i the part)
ABA_WORD_SIZE=20; export ABA_WORD_SIZE # This sets the word size. This is critical for speed issues in nucmer. Default is 20
ABA_COMPARISON=promer; export ABA_COMPARISON # This sets promer as the comparison method to use (default). It can be changed to 'nucmer'
Check the script 'abacas2.showACT.sh' if you want to automatically check the comparisons in the Artemis Comparison Tool (ACT)

### Results
Results are all files beginning with 'Res.*'. This includes GFF files for contig placements, coverage plots to be loaded into Artemis, FASTA sequences for each contiguated pseudochromosome and the 'Res.abacasBin.fna' file for all input sequences that could not be mapped to the reference (aka the 'bin').

### License
ABACAS2 is released under GPL3.

### Citation
PMID: 19497936
"
	exit;
fi


### Check and set arguments by default if not provided
if [ -z "$ABA_MIN_LENGTH" ] ; then
	ABA_MIN_LENGTH=200;
fi
if [ -z "$ABA_MIN_IDENTITY" ] ; then
	ABA_MIN_IDENTITY=95;
fi
if [ -z "$ABA_CHECK_OVERLAP" ] ; then
        ABA_CHECK_OVERLAP=0;
	export ABA_CHECK_OVERLAP
fi
if [ "$ABA_MIN_IDENTITY" -gt "99" ] ; then
	echo "Your identity might be too high $ABA_MIN_IDENTITY > 99 "
	exit;
fi

tmp=$$
sed 's/|/_/g' $reference > Ref.$tmp
reference=Ref.$tmp
ln -s $contig Contigs.$tmp
contig=Contigs.$tmp
export ABA_MIN_LENGTH ABA_MIN_IDENTITY contig reference


### ABACAS2 runComparison
abacas2.runComparison.sh $reference $contig $cores $ABA_splitContigs


### ABACAS2 TillingGraph
for x in `grep '>' $reference | perl -nle '/>(\S+)/;print $1' ` ; do
	abacas2.doTilingGraph.pl $x.coords $contig Res
done


### ABACAS2 binning
# abacas2.bin.sh $contig Res.abacasBin.fna && grep -v '>'  Res.abacasBin.fna | awk 'BEGIN {print ">Bin.union"} {print}' > Res.abacasBin.oneSeq.fna_ && cat Res*fna > Genome.abacas.fasta && bam.correctLineLength.sh Genome.abacas.fasta  &> /dev/null && mv  Res.abacasBin.fna_ Res.abacasBin.fna
# JLRuiz update 2022. The line above is the original one in Abacas. However, there was always an error regarding the last statement: "mv Res.abacasBin.fna_ Res.abacasBin.fna"
# It seems "Res.abacasBin.fna_" had never been defined, nor "Res.abacasBin.fna" seemed to be used after that attempted mv command. I'm not sure on the aim, but I have commented the original line and added the same without the last mv statement so we avoid the errors
# In any case, the contigs that are identified in the bin by abacas2 are included in the output sequences, just not merged into a single one. If this is something desirable, here is where the code has to be modified

abacas2.bin.sh $contig Res.abacasBin.fna
grep -v '>' Res.abacasBin.fna | awk 'BEGIN {print ">Bin.union"} {print}' > Res.abacasBin.oneSeq.fna_
cat Res*fna > Genome.abacas.fasta
bam.correctLineLength.sh Genome.abacas.fasta &> /dev/null


### Blasting:
ref=$reference
pre=Res
if [ -z "$pre" ] ; then
   echo "usage <reference> Res"

fi

tmp=$$
ln -s $ref REF.$tmp

mkdir Reference; cd Reference
SeparateSequences.pl ../REF.$tmp
cd ..; mkdir comp

count=0
for nameRes in `grep '>' $ref | perl -nle 's/\|/_/g;/>(\S+)/; print $1'` ; do
	let count++;
	if [ $count -gt 200 ] ; then
		echo "too many contigs to bsub!\n";
		exit
	fi
	formatdb -p F -i Reference/$nameRes
	megablast -F T -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -a $cores -o comp/comp.$nameRes.blast
done


