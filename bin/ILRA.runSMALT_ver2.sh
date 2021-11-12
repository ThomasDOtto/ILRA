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

### Modified by José Luis Ruiz to be used as part of the ILRA pipeline in 2022

### Check and set arguments
genome=$1
kmer=$2
stepsize=$3
fastqF=$4
fastqR=$5
resultname=$6
insertSize=$7
threads=$8
tmp=$$

if [ -z "$genome" ] || [ -z "$kmer" ] || [ -z "$stepsize" ] || [ -z "$fastqF" ] || [ -z "$insertSize" ] ; then
	echo "Usage: <genome> <k-mer> <stepsize> <Illumina forward reads> <Illumina reverse reads; 0 if single Reads> <resultname> <Insertsize Range:i.e. 350, 0 if single> optional: for more specific parameter, set the SMALT_PARAMETER enrovimental variable"
	exit;
fi


### Executing smalt...
smalt index -k $kmer -s $stepsize Genome.index.$tmp $genome > out.$tmp.txt
echo -e "\n\nExtra parameter (SMALT_PARAMETER): $SMALT_PARAMETER"
if [ "$fastqR" = "0" ] ; then
	smalt map $SMALT_PARAMETER -f samsoft Genome.index.$tmp $fastqF | samtools view -@ threads -f 0x2 -u - | samtools sort -@ $threads -u -m 2G --write-index -o $resultname.bam##idx##$resultname.bam.bai - >> out.$tmp.txt
else
	smalt map $SMALT_PARAMETER -f samsoft -i $insertSize Genome.index.$tmp $fastqF $fastqR | samtools view -@ threads -f 0x2 -u - | samtools sort -@ $threads -u -m 2G --write-index -o $resultname.bam##idx##$resultname.bam.bai - >> out.$resultname.$tmp.txt
	smaltOK=$?
	if [ $smaltOK -ne 0 ] ; then
		echo "Problems with SMALT: $smaltOK"
		exit $smaltOK
	fi
fi

