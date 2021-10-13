#!/bin/bash
# Modified by José Luis Ruiz from https://github.com/wudustan/fastq-info to be used as part of the ILRA pipeline in 2022

# Version 1.2
# Original author: Raymond Kiu Raymond.Kiu@quadram.ac.uk

fastq1=$1  #input 1
fastq2=$2  #input 2
fasta=$3   # input 3

#To calculate read lenghts 
readlength1=$(awk 'NR % 4 == 2 { s += length($1); t++} END {print s/t}' <(zcat $fastq1)) #calculate read length for fastq file 1
readlength2=$(awk 'NR % 4 == 2 { s += length($1); t++} END {print s/t}' <(zcat $fastq2)) #calculate read length for fastq file 2
totalreadlength=$(echo "$readlength1+$readlength2"|bc) #calculate total read length
readlength=$(echo "$totalreadlength/2"|bc) #calculate average read length by parsing into bc command - which truncates the decimal numbers by default in multiplication/division

#To calculate read counts
readcount1=$(zcat $fastq1| echo $((`wc -l`/4))) #read count in file 1
readcount2=$(zcat $fastq2| echo $((`wc -l`/4))) #read count in file 2
totalreadcount=$(echo "$readcount1+$readcount2"|bc) #total read count of file 1 and file 2
readcount=$(echo "$totalreadcount/2"|bc) #average read count - parsed into bc command to be rounded

#To calculate genome size
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $fasta > "$fasta"-single
cat "$fasta"-single|awk 'NR%2==0'| awk '{print length($1)}' > "$fasta"-oneline
genomesize=$(awk '{ sum += $1 } END { print sum }' "$fasta"-oneline)

#sequencing depth= LN/G = read length * 2(paired-end) * average reads / genome size
LN=$(echo "$readlength*2*$readcount")
G=$genomesize
coverage=$(echo "$LN/$G"|bc)

#print outputs
echo -n -e "sample_name\t"
echo -n -e "\tread_length\t"
echo -n -e "reads\t"
echo -n -e "genome\t"
echo -e "coverage"
echo -n -e "$fasta\t"
echo -n -e "$readlength\t"
echo -n -e "$readcount\t"
echo -n -e "$genomesize\t"
echo -e "$coverage"

#remove intermediary files
rm "$fasta"-single;
rm "$fasta"-oneline;

