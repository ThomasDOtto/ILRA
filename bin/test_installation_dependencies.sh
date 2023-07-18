#!/bin/bash

if [[ $debug == "all" || $debug == "step1" ]]; then
	type ILRA.removesmalls.pl >/dev/null 2>&1 || { echo >&2 "I require the ILRA folder, all files within, the ILRA/bin subfolder, the ILRA/external_software subfolder... to be in the PATH. Aborting..."; exit 1; }
	type assembly-stats >/dev/null 2>&1 || { echo >&2 "I require assembly-stats but it's not installed or available in the PATH. Aborting..."; exit 1; }	
elif [[ $debug == "all" || $debug == "step2" ]]; then
	type formatdb >/dev/null 2>&1 || { echo >&2 "I require formatdb but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type megablast >/dev/null 2>&1 || { echo >&2 "I require megablast but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type parallel >/dev/null 2>&1 || { echo >&2 "I require GNU's parallel but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type smalt >/dev/null 2>&1 || { echo >&2 "I require smalt but it's not installed or available in the PATH. Aborting..."; exit 1; }
 	type time >/dev/null 2>&1 || { echo >&2 "I require time but it's not installed or available in the PATH. Aborting..."; exit 1; }
elif [[ $debug == "all" || $debug == "step3" ]]; then
	type abacas2.nonparallel.sh >/dev/null 2>&1 || { echo >&2 "I require ABACAS2 and its dependencies but something is not installed or available in the PATH. Aborting..."; exit 1; }
	type fastaq >/dev/null 2>&1 || { echo >&2 "I require Fastaq but it's not installed or available in the PATH. Aborting..."; exit 1; }
elif [[ $debug == "all" || $debug == "step4" ]]; then
	type icorn2.serial_bowtie2.sh >/dev/null 2>&1 || { echo >&2 "I require iCORN2 and its dependencies but something is not installed or available in the PATH. Aborting..."; exit 1; }
	type bowtie2 >/dev/null 2>&1 || { echo >&2 "I require Bowtie2 but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type java >/dev/null 2>&1 || { echo >&2 "I require java, particularly v1.7.0 for iCORN2, but it's not installed or available in the PATH. Aborting..."; exit 1; }
elif [[ $debug == "all" || $debug == "step5" ]]; then
	type winnowmap >/dev/null 2>&1 || { echo >&2 "I require winnowmap but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type meryl >/dev/null 2>&1 || { echo >&2 "I require meryl but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type circlator >/dev/null 2>&1 || { echo >&2 "I require Circlator but it's not installed or available in the PATH. Aborting..."; exit 1; }
elif [[ $debug == "all" || $debug == "step6" ]]; then
	type kraken2 >/dev/null 2>&1 || { echo >&2 "I require kraken2 but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type rextract >/dev/null 2>&1 || { echo >&2 "I require recentrifuge but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type fasta_to_fastq.pl >/dev/null 2>&1 || { echo >&2 "I require fasta_to_fastq.pl but it's not installed or available in the PATH. Aborting..."; exit 1; }
        type vecscreen >/dev/null 2>&1 || { echo >&2 "I require vecscreen but it's not installed or available in the PATH. Aborting..."; exit 1; }
        type VSlistTo1HitPerLine.sh >/dev/null 2>&1 || { echo >&2 "I require VSlistTo1HitPerLine.sh but it's not installed or available in the PATH. Aborting..."; exit 1; }
        type blastdbcmd >/dev/null 2>&1 || { echo >&2 "I require blastdbcmd but it's not installed or available in the PATH. Aborting..."; exit 1; }
        type blastn >/dev/null 2>&1 || { echo >&2 "I require blastn but it's not installed or available in the PATH. Aborting..."; exit 1; }
        type bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools but it's not installed or available in the PATH. Aborting..."; exit 1; }
elif [[ $debug == "all" || $debug == "step7" || $quality_step == "yes" ]]; then
	type gtool.py >/dev/null 2>&1 || { echo >&2 "I require gtool.py but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type fastq-info.sh >/dev/null 2>&1 || { echo >&2 "I require fastq-info.sh but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type pigz >/dev/null 2>&1 || { echo >&2 "I require pigz but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type fastqc >/dev/null 2>&1 || { echo >&2 "I require fastqc but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type quast.py >/dev/null 2>&1 || { echo >&2 "I require quast but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type busco >/dev/null 2>&1 || { echo >&2 "I require busco but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type meryl >/dev/null 2>&1 || { echo >&2 "I require meryl but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type syri >/dev/null 2>&1 || { echo >&2 "I require syri but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type plotsr >/dev/null 2>&1 || { echo >&2 "I require plotsr but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type asm2stats.pl >/dev/null 2>&1 || { echo >&2 "I require asm2stats.pl but it's not installed or available in the PATH. Aborting..."; exit 1; }
        type asm2stats.minmaxgc.pl >/dev/null 2>&1 || { echo >&2 "I require asm2stats.minmaxgc.pl but it's not installed or available in the PATH. Aborting..."; exit 1; }
fi
