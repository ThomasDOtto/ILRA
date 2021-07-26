#!/bin/bash

#
# File: IPA.runSMALT.pl
# Time-stamp: old
# $Id: $
#
# (c) Copyright: Genome Research Ltd. 2017
#
# Author: Thomas Dan Otto
# Licencse: GNU General Public License v3.0
#
# Description:  little script to align with maq a read pair Illumina/Solexa lane onto a genome to produce a bam file for Artemis.
# samtools and smalt must be in the PATH

genome=$1
kmer=$2
stepsize=$3
fastqF=$4
fastqR=$5
resultname=$6
insertSize=$7
threads=$8

#randome stamp for temporary file name
tmp=$$

if [ -z "$insertSize" ] ;
		then
		echo "Usage: <genome> <k-mer> <stepsize> <Illumina forward reads> <Illumina reverse reads; 0 if single Reads> <resultname> <Insertsize Range:i.e. 350, 0 if single> optional: for more specific parameter, set the SMALT_PARAMETER enrovimental variable"
		echo
		exit;
	fi;


smalt index -k $kmer -s $stepsize Genome.index.$tmp $genome > out.$tmp.txt

echo "Extra parameter (SMALT_PARAMETER): $SMALT_PARAMETER";
if [ "$fastqR" = "0" ] ; then
  smalt map $SMALT_PARAMETER -n 6  -f samsoft  -o $resultname.$tmp.sam Genome.index.$tmp $fastqF >> out.$tmp.txt


else

# Check if the insert size is in the correct format:


# Generating sam with SSAHA2 (>= version 2.5, must be in path)

	echo "	smalt map $SMALT_PARAMETER  -f samsoft -i $insertSize -o $resultname.$tmp.sam Genome.index.$tmp $fastqF $fastqR >> out.$resultname.$tmp.txt"

	smalt map $SMALT_PARAMETER    -f samsoft -i $insertSize -o $resultname.$tmp.sam Genome.index.$tmp $fastqF $fastqR >> out.$resultname.$tmp.txt
	smaltOK=$?
	if [ $smaltOK -ne 0 ] ; then
		echo "Problems with SMALT: $smaltOK"
		exit $smaltOK
	fi

fi
# rm  Genome.index.$tmp.*
if [ -z "$tmp.sam" ] ; then
	echo "Problem generating SAM file.";
	exit;
fi

# Modified by JlRuiz 2021

# Convert reference to fai for bam file
# samtools faidx $genome


# SAM to BAM
# samtools import $genome.fai $resultname.$tmp.sam $resultname.$tmp.bam
#samtools view -S -b -h $resultname.$tmp.sam > $resultname.test.$tmp.bam
#rm $resultname.$tmp.sam
#order the bam file
samtools sort -@ $threads -o $resultname.bam $resultname.$tmp.sam
rm $resultname.$tmp.sam
#index the bam file, to get the bai file.
samtools index -@ $threads $resultname.bam

# rm out.$resultname.$tmp.txt  $resultname.$tmp.sam
