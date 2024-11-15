#!/bin/bash
# Authors: José Luis Ruiz.pk & Thomas Dan Otto
# License: GNU General Public License v3.0
# https://github.com/ThomasDOtto/ILRA

########## ILRA
base64 -d <<<"H4sIAAAAAAAAA7VVyUoEMRC95yvq5ChoR+YbPIkgiHgK1HjTg8zB8SL18aa2LE2PjJ0xTaeTl6pXL1t1AEAusLZ0nsJznGsx0imh0WoMQBKC/qSxFoLGMzJPh8xseT5pmSCmBY/OhuASr4J3kblk9tIqfV0x+VYTKF3S+QqgbRfQcImfjCBecYUOIZp5BJkIFhFtPHE0eh914QQbhF22TdKSGiP38xNBW25Su9mOYytA0Cx0rFzc3vnyQJQRhWyAraMGkqBIfTyXREDujBQ0mkW0h9dFnrwnmeCFF4Qpq0nbRdu9spO8phVRHm+xPnFWSGtRzKcncaXHqI9XJfGAmASJkyO5/vxOeZC/TJL9k7ydSdNNWCzFx/iwIGpzDeqcLcWZIa8h5gjgQxFVWxsPTZJF5hAx/H7RzlRkARdCrc4GR8uRSZ1SpptZmUawVRpoTjQBjWArRZQ8pTmfmJxapGJ6TLHD1LnYjYio/wwnb/8ihpFl0jQTYalzWERtVPJ6SwrmWbyzszT4nyKwF8ErETs7z2hnEIH+W6tnQqYdZWv8TEidp40Fw6iaBBw+E5JqoT0TVqg5E9LV9K8YuYV8h7ejAUawlSI2/V3fMPl6bF3CnPEw0wAW7h+evt6/4QKe3/Yfr59wN8Hj4bCfYHu73YYfi5OSru8KAAA=" | gunzip
echo -e "\n\n\n\nCurrent time: $(date)\n"
# Pipeline for finishing and the Improvement of Long Reads Assemblies
start=`date +%s`

##### Arguments / Variables
# A string with command options and an array with arguments
options=$@
arguments=($options)
# Get arguments by looping through an index
index=0
for argument in $options; do
# Incrementing index
	index=`expr $index + 1`
# Gather the parameters
	case $argument in
		-h*) echo "ILRA v1.5.3. usage: ILRA.sh [options]
		-h | -help # Type this to get help
		-a | -assembly # Name of the long reads assembly to correct (FASTA format, can be gzipped)
		-f | -filter_contig_size # Size threshold (bp) to filter the contigs (by default, 5000)
		-I | -Illumina_reads # Root name of the paired-end Illumina reads (FASTQ format, must be gzipped and name root_1.fastq.gz and root_2.fastq.gz prior execution)
		-C | -Correction_Illumina_reads # Whether illumina reads are provided and all steps of correction should be performed ('no' /'yes' by default)
		-c | -corrected_reads # Corrected long reads for circularization of the contigs containing the strings by -s and -S (FASTQ or fasta format, can be gzipped, and Circlator is going to be skipped if not provided)
		-R | -Range_insert_size # Insert size range (bp) of the Illumina reads to use in mapping steps (by default, 800)
		-i | -iterations_polishing # Number of iterations to perform in polishing by iCORN2 or Pilon (by default, 3)
		-r | -reference # Reference file (full pathway, FASTA format, if not provided ABACAS2 and some quality control steps are going to be skipped)
		-g | -gff_file # Reference annotation file (full pathway, GFF format, if not provided some quality control steps are going to be skipped)		
		-s | -seq_circularize # Regex pattern to find in the contig names and circularize (by default 'MT|M76611|API')
		-L | -Long_reads_technology # Technology of long reads sequencing ('ont' or 'pb', by default)
		-T | -Taxon_ID # NCBI taxon ID or a scientific name (Plasmodium by default, 5820)
		-t | -threads # Number of cores to use in multithreaded steps (by default, 40)
		-e | -ending_telomere_seq_1 # Telomere-associated sequence to search in the first strand (by default, 'CCCTAAACCCTAAACCCTAAA')
		-E | -Ending_telomere_seq_2 # Telomere-associated sequence to search in the second strand (by default, 'TTTAGGGTTTAGGGTTTAGGG')
		-o | -output # Output folder (absolute pathway)
		-n | -name # Base name of the output file		
		-d | -debug_step # For debug, step to remove the content of the corresponding folder and resume a failed run ('step1', 'step2a', 'step2b', 'step3', 'step4', 'step4i', 'step5', 'step6', 'step7' or 'all' to execute everything, by default)
		-D | -databases # Folder for storing the databases for the decontamination step (by default, the 'databases' folder under ILRA main folder)
		-K | -Kraken2_fast_mode # Kraken2 fast mode, consisting on copying the Kraken2 database to /dev/shm (RAM) so execution is faster ('yes' / 'no' by default)
		-k | -Kraken2_databases # Folder within the folder databases (-D) containing the database used by Kraken2 (by default, 'standard_eupathdb_48_kraken2_db')
		-b | -block_size # Block size for parallel processing (by default, 5)
		-Bi | -blast_block_size # Block size for the initial parallel processing in the first megaBLAST step (by default the argument -b, block_size, but may be necessary to change since megaBLAST can be less computationally expensive)
		-Bl | -blast_block_size_large # Block size for parallel processing in the first megaBLAST step for the larger contigs whose processing failed due to RAM requirements (by default 2, 1 would mean sequential processing of the large contigs and should not fail)
		-As | -abacas2_split # Number of parts to split and process in parallel in ABACAS2 (by default the argument -b, block_size, but may be necessary to decrease due to memory issues)
		-Ab | -abacas2_blast # Whether to do blast within ABACAS2 to compare with the reference and display in ACT (1 by default, which means blasting, or 0)
		-Ol | -overlap_length # Threshold in the length of overlap to consider a contig contained into other (by default, 2000 bp)
		-Oi | -overlap_identity # Threshold in the percentage of identity to consider a contig contained into other (by default, 99)
		-Of | -overlap_fraction # Threshold in the fraction (percentage) of a contig contained into other to consider overlapping (by default, 90)
		-F | -filter_using_Illumina # Do you have even Illumina reads coverage? If yes contigs in the long reads assembly not covered by Illumina short reads by a certain threshold specified in other parameters will be filtered out ('no' /'yes' by default)		
		-Mc | -merging_coverage_threshold # Threshold in the fraction of mean genome coverage (percentage) of Illumina short reads at an overlap between contigs to consider their merging (by default, 0.5)
		-Md | -merging_coverage_deviation # Positive deviation to sum to the fraction of mean genome coverage (percentage) of Illumina short reads at an overlap between contigs to consider their merging (by default, 0.1)
		-p | -pilon # Whether to use pilon instead of iCORN2 ('yes'/'no' by default)
		-P | -parts_icorn2_split # Number of parts to split the input sequences of iCORN2 before processing them (0 by default, which means no splitting)
		-q | -quality_assesment # Whether to execute the final step 7, which may be slow, for assessing the quality of the corrected assembly, gathering sequences, analyzing telomeres... ('no'/'yes' by default)
		-Q | -BUSCO database for quality assessment # The name of the BUSCO database to be used in the quality assessment step. Automatic lineage selected by default if user does not input one of the datasets in 'busco --list-datasets' here (e.g. bacteria_odb10)
		-Mj | -java_memory # Max Java memory (heap space) to be used ('XXg', by default 240g=240GB used)
		-lm | -low_memory # Activate low memory mode for iCORN2 ('yes'/'no' by default)
		-ls | -low_space # Activate low space mode for iCORN2 and some steps (e.g., uncompressing fastq reads, iterative mapping, MarkDuplicatesSpark, etc) are processed in /dev/shm and occupy RAM ('yes'/'no' by default)
		-m | -mode # Add 'taxon' to execute decontamination based on taxonomic classification by kraken2, add 'blast' to execute decontamination based on BLAST against databases as requested by the DDBJ/ENA/Genbank submission, add 'both' to execute both approaches, and add 'light' to execute ILRA in light mode and skip these steps (this is the default option)" && exit 1;;
		-a*) assembly=${arguments[index]} ;;
		-o*) dir=${arguments[index]} ;;
		-c*) correctedReads=${arguments[index]} ;;
		-n*) name=${arguments[index]} ;;
		-r*) reference=${arguments[index]} ;;
		-I*) illuminaReads=${arguments[index]} ;;
		-C*) perform_correction=${arguments[index]} ;;
		-s*) seqs_circl=${arguments[index]} ;;
		-i*) number_iterations_icorn=${arguments[index]} ;;
		-t*) cores=${arguments[index]} ;;
		-f*) contigs_threshold_size=${arguments[index]} ;;
		-F*) contigs_illumina_filter=${arguments[index]} ;;
		-R*) InsertsizeRange=${arguments[index]} ;;
		-T*) taxonid=${arguments[index]} ;;
		-g*) gff_file=${arguments[index]} ;;
		-e*) telomere_seq_1=${arguments[index]} ;;
		-E*) telomere_seq_2=${arguments[index]} ;;
		-L*) seq_technology=${arguments[index]} ;;
		-m*) mode=${arguments[index]} ;;
		-Mj*) java_memory=${arguments[index]} ;;
		-Ol*) overlap_length=${arguments[index]} ;;
		-Oi*) overlap_identity=${arguments[index]} ;;
		-Of*) overlap_fraction=${arguments[index]} ;;
		-Mc*) merging_coverage_threshold=${arguments[index]} ;;
		-Md*) merging_coverage_deviation=${arguments[index]} ;;
		-lm*) low_mem=${arguments[index]} ;;
  		-ls*) low_spa=${arguments[index]} ;;
		-d*) debug=${arguments[index]} ;;
		-D*) databases=${arguments[index]} ;;
		-K*) kraken2_fast=${arguments[index]} ;;
		-k*) kraken2_databases=${arguments[index]} ;;
		-b*) blocks_size=${arguments[index]} ;;
		-Ab*) abacas2_blast=${arguments[index]} ;;
		-Bi*) blast_blocks_size=${arguments[index]} ;;
		-Bl*) blast_blocks_size_large=${arguments[index]} ;;
		-p*) pilon=${arguments[index]} ;;
		-P*) parts_icorn2_split=${arguments[index]} ;;
		-As*) abacas2_split=${arguments[index]} ;;
		-q*) quality_step=${arguments[index]} ;;
		-Q*) quality_step_busco_db=${arguments[index]} ;;
	esac
done
export name; export telomere_seq_1; export telomere_seq_2

echo -e "Check help with parameter '-h' or read the GitHub README for full example usage"
echo -e "In case you want to test ILRA with a smaller subset of example reads, check the subfolder /path/to/ILRA/test_data and the 'Quick Start' in the GitHub README"
# echo -e "Current date/time: $(date)\n"

##### Checking Arguments / Variables:
echo -e "\nI'm now quickly checking and showing the arguments that are going to be used in the ILRA run...\n"
if [[ -d "$dir/1.Filtering" ]]; then
	echo -e "PLEASE keep in mind the output directory already exists (found the folder '$dir/1.Filtering', so this is potentially a rerunning. ILRA is going to clean the existing folders for each step and rerun. If you want to run ILRA only from a particular step, please use the argument '-d '...\n"
 	echo -e "If the process is not aborted, ILRA is proceeding in..."
	secs=$((1 * 30))
	while [ $secs -gt 0 ]; do
	   echo -ne "$secs\033[0K\r"
	   sleep 1
	   : $((secs--))
	done
fi

if [ -z "$cores" ]; then
	cores=40
	echo "Number of threads: "$cores
fi

if [ -z "$java_memory" ]; then
	java_memory=240g
	echo "Number of threads: "$cores
fi
export _JAVA_OPTIONS="-Xms5g -Xmx$java_memory"

if [ -z "$low_mem" ]; then
	low_mem="no"
fi

if [ -z "$low_spa" ]; then
	low_spa="no"
fi

if [ -z "$parts_icorn2_split" ]; then
	parts_icorn2_split=0
fi

if [[ $assembly == /* ]]; then
	echo -e "Assembly provided correctly"
elif [ -z "$assembly" ]; then
	echo -e "ASSEMBLY NOT PROVIDED, exiting... Check out ILRA.sh -h"
	exit 1
else
	echo -e "Assuming the assembly is in the current pathway... If errors please provide absolute pathways as arguments"
	assembly=$PWD/$assembly
fi

if [[ $dir == /* ]]; then
	echo -e "Setting up the directory to work..."
else
	dir=$PWD/$dir
fi

if [[ $correctedReads == /* ]]; then
	echo -e "Corrected long reads provided correctly"
elif [ -z "$correctedReads" ]; then
	echo -e "CORRECTED READS NOT PROVIDED, Circlator is going to be skipped"
	echo -e "If errors, please provide absolute pathways as arguments"
else
	correctedReads=$PWD/$correctedReads
	echo -e "Assuming corrected reads are in the current pathway..."
	echo -e "If errors, please provide absolute pathways as arguments"
fi

if [ -z "$perform_correction" ]; then
	perform_correction="yes"
elif [ $perform_correction == "no" ]; then
	echo -e "You have deactivated the filtering out of contigs in the assembly not covered by short Illumina reads..."
	perform_correction=""
fi

if [ -z "$contigs_illumina_filter" ]; then
	contigs_illumina_filter="yes"
elif [ $contigs_illumina_filter == "no" ]; then
	echo -e "You have deactivated correction based on using short reads to find and remove overlapping contigs..."	
fi

if [ -z "$blocks_size" ]; then
	blocks_size=5
fi

if [ -z "$abacas2_split" ]; then
	abacas2_split=$blocks_size
fi
if [ -z "$blast_blocks_size" ]; then
	blast_blocks_size=$blocks_size
fi
if [ -z "$blast_blocks_size_large" ]; then
	blast_blocks_size_large=2
fi
if [ -z "$abacas2_blast" ]; then
	abacas2_blast=1
fi

if [[ -z "$illuminaReads" ]] && [[ $perform_correction == "yes" ]]; then
	echo -e "YOU HAVE NOT PROVIDED ILLUMINA READS. Please rerun the pipeline with the argument "-C no" if you don't want to use Illumina reads or perform correction. ILRA will then skip some steps accordingly..."
	exit 1
elif [[ $illuminaReads == /* ]]; then
	echo -e "Checking Illumina short reads..."
elif [[ $illuminaReads != /* ]]; then
	echo -e "Not sure where Illumina reads are... Assuming they are in the current pathway... If errors please provide absolute pathways as arguments"
	illuminaReads=$PWD/$illuminaReads
fi

if [ -f $illuminaReads\_1.fastq.gz ]; then
	echo "Good, ILRA is detecting the naming required for the Illumina reads: _1.fastq.gz and _2.fastq.gz..."
else
	if [ "$perform_correction" == "yes" ]; then
		echo -e "PLEASE be aware of the naming required for the Illumina reads: _1.fastq.gz and _2.fastq.gz"
		echo -e "ILRA is not detecting the Illumina reads files as you provided them. If you want to use Illumina reads for polishing, please check naming, paths and that the files do exist. Otherwise rerun the pipeline with the argument "-C no" if you don't want to use Illumina reads or perform correction. ILRA will skip some steps accordingly..."
		exit 1
	fi
fi

if [ -z "$reference" ]; then
	echo "NO REFERENCE GENOME PROVIDED. ABACAS2 will be skipped ..."
	doAbacas2=0
else
	doAbacas2=1
	if [[ $reference == /* ]]; then
		echo -e "Reference provided correctly"
	else
		echo -e "Not sure where reference is... Assuming it is in the current pathway... If errors please provide absolute pathways as arguments"
		reference=$PWD/$reference	
	fi
fi

if [ ! -z "$gff_file" ]; then
	if [[ $gff_file == /* ]]; then
		echo -e "GFF file with genes provided correctly"
	else
		echo -e "Not sure where gff file is... Assuming it is in the current pathway... If errors please provide absolute pathways as arguments"
		gff_file=$PWD/$gff_file	
	fi
fi

if [ -z "$mode" ]; then
	mode="light"
	echo "LIGHT MODE activated, steps for decontamination and preparation for online databases steps will be skipped ..."
else
	echo "ILRA execution mode (-m) is: "$mode
fi

if [ -z "$kraken2_fast" ]; then
	kraken2_fast="no"
fi

if [ -z "$kraken2_databases" ]; then
	kraken2_databases="standard_eupathdb_48_kraken2_db"
fi

if [ -z "$pilon" ]; then
	pilon="no"
	echo "You are using iCORN2 for the correction via Illumina short reads. If you want to use pilon, please provide the argument '-p yes'..."
elif [[ $pilon == "yes" ]]; then
	echo "You are using pilon for the correction via Illumina short reads ..."
else
	pilon="no"
	echo "You are using iCORN2 for the correction via Illumina short reads. If you want to use pilon, please provide the argument '-p yes'..."
fi

if [ -z "$quality_step" ] || [ $quality_step == "yes" ]; then
	quality_step="yes"
	echo "You are requesting the last step for assessing the quality of the corrected assembly, gathering sequences, analyzing telomeres... etc"
elif [ $quality_step == "no" ]; then
	echo "You are skipping the last step for assessing the quality of the corrected assembly, gathering sequences, analyzing telomeres... etc"
fi
export quality_step

echo -e "\nFinal arguments used:"
if [ -z "$number_iterations_icorn" ]; then
	number_iterations_icorn=3
	echo "Number of iCORN2 iterations: "$number_iterations_icorn
fi

if [ -z "$seqs_circl" ]; then
	seqs_circl="MT|M76611|API"
	echo "Seqs to circularize: "$seqs_circl
fi

if [ -z "$contigs_threshold_size" ]; then
	contigs_threshold_size=5000
	echo "Length threshold to discard contigs (bp): "$contigs_threshold_size
fi

if [ -z "$InsertsizeRange" ]; then
	InsertsizeRange=800
	echo "Insert size range Illumina short reads (bp): "$InsertsizeRange
fi

if [ -z "$taxonid" ]; then
	taxonid=5820
	echo "NCBI taxon id: "$taxonid
fi

if [ -z "$databases" ]; then
	databases=$(dirname $0)/databases
fi

if [ $mode == "light" ]; then
	echo "Execution mode is $mode. Databases are/should be located in the folder $databases"
elif [[ $mode == "blast" || $mode == "taxon" || $mode == "both" ]]; then
	mkdir -p $databases
	echo "Execution mode is $mode. Databases are/should be located in the folder $databases"
fi

if [ $seq_technology == "pb" ]; then
	long_reads_technology="pacbio"
	echo "Long reads sequencing technology: PacBio"
	echo "If Oxford Nanopore were used, please provide the argument 'ont'"
elif [ $seq_technology == "ont" ]; then
	long_reads_technology="ont2d"
	echo "Long reads sequencing technology: ONT"
	echo "If PacBio were used, please provide the argument 'pb'"
else
	long_reads_technology="pacbio"
	echo "Using by default long reads sequencing technology: PacBio"
	echo "If this needs to be change please provide 'ont' as the last argument"
fi

if [ -z "$telomere_seq_1" ]; then
	telomere_seq_1="CCCTAAACCCTAAACCCTAAA"
fi

if [ -z "$telomere_seq_2" ]; then
	telomere_seq_2="TTTAGGGTTTAGGGTTTAGGG"
fi

if [ -z "$debug" ]; then
	debug="all"
fi
export debug

if [[ $assembly = *.gz ]]; then
	mkdir -p $dir; pigz -p $cores -dc $assembly > $dir/$(basename ${assembly%.*})
	assembly=$dir/$(basename ${assembly%.*})
fi

if [[ $blocks_size -gt $(grep -c ">" $assembly) ]]; then
	blocks_size=$(grep -c ">" $assembly)
fi

if [ -z "$overlap_length" ]; then
	overlap_length=2000
fi
if [ -z "$overlap_identity" ]; then
	overlap_identity=99
fi
if [ -z "$overlap_fraction" ]; then
	overlap_fraction=90
fi
if [ -z "$merging_coverage_threshold" ]; then
	merging_coverage_threshold=0.5
fi
if [ -z "$merging_coverage_deviation" ]; then
	merging_coverage_deviation=0.1
fi

echo -e "assembly="$assembly
echo -e "dir="$dir
if [[ $correctedReads == /* ]]; then
	echo -e "correctedReads="$correctedReads
fi
echo -e "name="$name
echo -e "reference="$reference
echo -e "illuminaReads="$illuminaReads
echo -e "cores="$cores
echo -e "cores_SMALT="$cores
echo -e "cores_iCORN2="$cores
echo -e "seqs_circl="$seqs_circl
echo -e "number_iterations_icorn="$number_iterations_icorn
echo -e "contigs_threshold_size="$contigs_threshold_size
echo -e "InsertsizeRange="$InsertsizeRange
echo -e "The telomere sequences used are:\nLeft:\t" $telomere_seq_1"\nRight:\t" $telomere_seq_2
echo -e "The block size for parallel processing is:" $blocks_size

echo -e "\nPATH with paths to all tools must be properly set: export PATH=...:\$PATH either in the bashrc or bash_profile files in HOME directory, activating a conda environment... etc"
echo -e "\nPlease be aware that many scripts (bash, perl...) within the ILRA folder are used, and you may need to manually change the interpreter in the corresponding shebang (first line #!) statements so everything works in your system. Many software and dependencies are also required and are being automatically checked at the beginning of each step execution. You can resume a step for which you did not have the software ready with -d argument. The pipeline will exit if any required software is not found in the variable PATH of your system, which you likely need to change accordingly... You may also need to make scripts executable (chmod), to make source, export... etc, so that the PATH variable and others are available for all scripts"
echo -e "If the pipeline keeps running, congrats! All the required software is available"

##### Checking if there's an interrupted ILRA run that should be resumed
if [[ $debug == "all" ]] && [[ -d "$dir" ]]; then
	step_to_resume=$(ls $dir | egrep '^[0-9].' | sort | tail -1 | sed 's,\..*,,g')
	echo -e "\n The directory $dir already exists, and ILRA is going to resume the step$step_to_resume in few seconds... If this is not correct, please stop ILRA with Ctrl+C and remove the directory, or specify another output dir with '-o'\n"
	secs=$((1 * 30))
	while [ $secs -gt 0 ]; do
		echo -ne "$secs\033[0K\r"
		sleep 1
		: $((secs--))
	done
	export debug="step$step_to_resume"
fi

##### Checking the installed databases:
if [[ $mode == "blast" || $mode == "taxon" || $mode == "both" ]]; then
	echo -e "\nPLEASE note some local databases for the decontamination step (kraken2 and blast according to DDBJ/ENA/Genbank requirements) are needed"
	echo -e "These databases are large, particularly the kraken2 database, so please be aware that the RAM memory usage at the step 6 of ILRA may reach hundreds of GBs. Rerun with '-m light' or '-m blast' to skip if not acceptable"
	echo -e "ILRA is now going to check or give you instructions so the databases are downloaded and placed in the corresponding folders"
	if [[ -f $databases/taxdump/names.dmp ]] && [[ -f $databases/taxdump/nodes.dmp ]] && [ "$(find -L $databases -name hash.k2d | wc -l)" -gt 0 ]; then
		echo -e "\nGood, ILRA is detecting all of the required databases for decontamination based on taxonomic classification in "$databases"\n"
	else
		echo -e "\nILRA will exit until these steps are performed. Please note 'wget' may not complete the download and then uncompressing would give errors. You would need to remove any incomplete file and restart download"
		echo -e "The databases names.dmp, nodes.dmp, merged.dmp and delnodes.dmp have to be downloaded by the user executing: mkdir -p $databases/taxdump && cd $databases/taxdump && wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip && unzip taxdmp.zip && rm taxdmp.zip"
		echo -e "The kraken2 database has to be built by the user and placed within the $databases folder, under a folder with your preferred name. Please execute the following and read these comments if you want to build a kraken2 v2.1.2 database containing the standard kraken2 recommended sequences and the EuPathDB v46 preprocessed and cleaned sequences (http://ccb.jhu.edu/data/eupathDB/):"
		echo -e "cd $databases && cores=64"
		echo -e "kraken2-build --download-taxonomy --threads \$cores --db standard_eupathdb_46_kraken2_db & #### The download commands can be sent to the background so they are simultaneously processed"
		echo -e "kraken2-build --download-library archaea --threads \$cores --db standard_eupathdb_46_kraken2_db &"
		echo -e "kraken2-build --download-library viral --threads \$cores --db standard_eupathdb_46_kraken2_db &"
		echo -e "kraken2-build --download-library plasmid --threads \$cores --db standard_eupathdb_46_kraken2_db &"
		echo -e "kraken2-build --download-library bacteria --threads \$cores --db standard_eupathdb_46_kraken2_db & #### This is around a download of 150GB"
		echo -e "kraken2-build --download-library fungi --threads \$cores --db standard_eupathdb_46_kraken2_db &"
		echo -e "kraken2-build --download-library UniVec_Core --threads \$cores --db standard_eupathdb_46_kraken2_db &"
		echo -e "kraken2-build --download-library protozoa --threads \$cores --db standard_eupathdb_46_kraken2_db &"
		echo -e "# kraken2-build --download-library nt --threads \$cores --db standard_eupathdb_46_kraken2_db & #### With the nt database results would more precise, but it's huge. It works if you are patient, but it would involve a download of ~450GB, and then building would require to allocate at least 1TB. Then the final database would be of ~450GB of size, and require that much RAM to build and run"
		echo -e "# kraken2-build --add-to-library chr1.fa --threads \$cores --db standard_eupathdb_46_kraken2_db #### If you want to add custom sequences or a particular genome you suspect contamination for... check kraken2 manual to handle taxonomy"
		echo -e "# To add EuPathDB:"
		echo -e "mkdir -p standard_eupathdb_46_kraken2_db/library/added && cd standard_eupathdb_46_kraken2_db/library/added"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/*DB46.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/seqid2taxid.map &"
		echo -e "for f in *.tgz; do tar -xvzf \$f; done && rm *.tgz \$(ls | grep log)"
		echo "awk '{ print "'"TAXID""\t"$1"\t"$2 }'"' seqid2taxid.map > prelim_map_1.txt && rm seqid2taxid.map"
		echo -e "# To build the kraken2 database:"
  		echo -e "cd $databases"
		echo -e "kraken2-build --build --threads \$cores --db standard_eupathdb_46_kraken2_db #### This would take around 5 hours and require around 75GB of RAM. If you included nt database, around 24 hours and 500GB of RAM. For the larger database including nt, you may need to include the flag --fast-build to avoid stalling or too much time building"
		echo -e "kraken2-inspect --db standard_eupathdb_46_kraken2_db --threads \$cores > standard_eupathdb_46_kraken2_db_k2_inspect.txt"
		echo -e "kraken2-build --clean --threads \$cores --db standard_eupathdb_46_kraken2_db #### This would remove the downloaded sequences and keep only the kraken2 database. This must be done to save space, but first double check the included sequences and the output of kraken2 inspect. If you execute this command and then you want to make any change to the database, you would have to download everything again"

		echo -e "\n\nIt is recommended to freshly prepare updated databases, but alternatively, you can download the already built older database: cd $databases && wget --no-check-certificate" '"bit.ly/kraken2_dbs" -O standard_eupathdb_48_kraken2_db.tar.xz && tar -xf standard_eupathdb_48_kraken2_db.tar.xz && rm standard_eupathdb_48_kraken2_db.tar.xz'
		echo -e "Please note that if you followed ILRA installation instructions of kraken2 v2.1.2, the kraken2-build script has been modified so database building is faster (improved masking, https://github.com/DerrickWood/kraken2/pull/39) and a bug in sequences download from ncbi has been corrected (https://github.com/DerrickWood/kraken2/issues/571)"
		echo -e "Please note that if you have allocated enough RAM, copying the kraken2 database to the faster RAM, something like /dev/shm/, and then pointing to the library there in kraken2 execution with the flag --memory-mapping would greatly improve speed, particularly if multiple runs (https://github.com/DerrickWood/kraken2/issues/451). ILRA includes this mode with the flag -K"
		echo -e "The execution of kraken2 classification with the suggested database (archaea + viral + plasmid + bacteria + fungi + UniVec_Core + EuPathDB) will require ~70GB of RAM."
		exit 1
	fi
	if [[ $mode == "blast" || $mode == "both" ]]; then
		echo -e "Several databases for conforming to DDBJ/ENA/Genbank requirements are needed, checking..."
	elif [[ -f $databases/contam_in_euks.fa ]] && [[ -f $databases/contam_in_prok.fa ]] && [[ -f $databases/adaptors_for_screening_euks.fa ]] && [[ -f $databases/adaptors_for_screening_proks.fa ]] && [[ -f $databases/mito.ndb ]] && [[ -f $databases/taxdb.btd ]] && [[ -f $databases/rrna ]]; then
		echo -e "\nGood, ILRA is detecting all of the required databases for decontamination based on blasting in "$databases"\n"
	else
		echo -e "Several databases for conforming to DDBJ/ENA/Genbank requirements are needed. Please download them and note that 'wget' may not complete the download and then uncompressing would give errors. You would need to remove any incomplete file and restart download. Please execute:"
		echo -e "cd $databases/"
		echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz && pigz -df -p 2 contam_in_euks.fa.gz && makeblastdb -in contam_in_euks.fa -dbtype nucl -out contam_in_euks.fa -title contam_in_euks.fa"
		echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_prok.fa && makeblastdb -in contam_in_prok.fa -dbtype nucl"
		echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_euks.fa && formatdb -p F -i adaptors_for_screening_euks.fa"
		echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_proks.fa && formatdb -p F -i adaptors_for_screening_proks.fa"
		echo -e "wget https://ftp.ncbi.nlm.nih.gov/blast/db/mito.tar.gz && tar -xvzf mito.tar.gz"
		echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/rrna.gz && pigz -df -p 2 rrna.gz && makeblastdb -in rrna -dbtype nucl -out rrna -title rrna"
		echo -e "Alternatively, to download all databases please execute: cd $databases && wget --no-check-certificate" '"bit.ly/ILRA_dbs" -O databases.tar.xz && tar -xf databases.tar.xz && rm databases.tar.xz'
		echo -e "This last option would download all databases at once, but you may prefer to fresh download them from NCBI, as they may be outdated in the compilation..."
		echo -e "ILRA is not detecting the required databases to decontaminate and you are not in the light mode, so the pipeline is exiting. Please double check the instructions just printed above and that the databases marked as not detected are in "$databases
		if [[ -f $databases/contam_in_euks.fa ]]; then echo "contam_in_euks.fa detected"; else echo "contam_in_euks.fa not detected"; fi
		if [[ -f $databases/contam_in_prok.fa ]]; then echo "contam_in_prok.fa detected"; else echo "contam_in_prok.fa not detected"; fi
		if [[ -f $databases/adaptors_for_screening_euks.fa ]]; then echo "adaptors_for_screening_euks.fa detected"; else echo "adaptors_for_screening_euks.fa not detected"; fi
		if [[ -f $databases/adaptors_for_screening_proks.fa ]]; then echo "adaptors_for_screening_proks.fa detected"; else echo "adaptors_for_screening_proks.fa not detected"; fi
		if [[ -f $databases/mito.ndb ]]; then echo "mito.ndb detected"; else echo "mito.ndb not detected"; fi
		if [[ -f $databases/taxdb.btd ]]; then echo "taxdb.btd detected"; else echo "taxdb.btd not detected"; fi
		if [[ -f $databases/rrna ]]; then echo "rrna detected"; else echo "rrna not detected"; fi
		ls -Rslh $databases
		exit 1
	fi
fi


##### ILRA execution STEPS:
# 1.  Discard smaller contigs
# 2.  megablast -F F
# 2a. Delete contained contigs
# 2b. Find overlaps
# 3.  ABACAS2 without overlap checking for reordering (resolve trivial gaps, reorder contigs, get the chromosome names, rename the sequences...)
# 4.  Error correction (iCORN2 or pilon)
# 5.  Circularization or organelles or any sequence reqired
# 6.  Decontamination. Taxonomic classification via kraken2/recentrifuge and final filtering for decontamination and databases upload (blast). Rename sequences if step 3 skipped
# 7.  Evaluate the assemblies, get telomere sequences counts, GC stats, sequencing depth, converting files...

#### 1. Discard contigs smaller than a threshold, convert to single line, and remove special characters and punctuation but ">" and "_" from the fasta headers, which is important for tools such as iCORN2:
if [[ $debug == "all" || $debug == "step1" ]]; then
	test_installation_dependencies.sh
	time1=`date +%s`
	echo -e "\n\nSTEP 1: Size filtering starting..."; echo -e "Current date/time: $(date)\n"
	mkdir -p $dir/1.Filtering; cd $dir/1.Filtering; rm -rf *
	echo -e "### Excluded contigs and length based on length threshold of $contigs_threshold_size: (ILRA.removesmalls.pl)" > ../Excluded.contigs.fofn
	ILRA.fasta2singleLine.pl $assembly | ILRA.removesmalls.pl $contigs_threshold_size - | sed 's/|/_/g' | ILRA.fasta2singleLine.pl - | awk -F ' ' '{ if ($0 ~ /^>/) { print $1;} else { print $0}}' | sed '/[[:punct:]]*/{s/[^[:alnum:][:space:]>_]/_/g}' > 01.assembly.fa
	short_contigs=$(grep -v "###" ../Excluded.contigs.fofn | cut -f1 | tr '\n' '|'); ILRA.removesmalls.R $PWD/01.assembly.fa $assembly $short_contigs $cores
	formatdb -p F -i $dir/1.Filtering/01.assembly.fa
	echo -e "\nPlease check the file 'Contigs_length_stats.txt' and the provided plots to assess whether the chosen length threshold is the most appropriate\n"
	echo -e "\nThe sequence of the $(echo $short_contigs | sed 's,|$,,g' | tr '|' '\n' | wc -l) discarded contigs aligned to the filtered assembly with: $(grep 'overall alignment rate' Contigs_length_stats.txt)"
 	echo "Before this step:"; assembly-stats $assembly | head -n 2
	echo -e "\nAfter this step:"; assembly-stats 01.assembly.fa | head -n 2
	echo -e "\nSTEP 1: DONE"; echo -e "Current date/time: $(date)"
	time2=`date +%s`; echo -e "STEP 1 time (secs): $((time2-time1))"; echo -e "STEP 1 time (hours): $(echo "scale=2; $((time2-time1))/3600" | bc -l)"
debug="all"
fi

#### 2. MegaBLAST
if [[ $debug == "all" || $debug == "step2a" ]]; then
	test_installation_dependencies.sh
	if [[ ! -d "$dir/1.Filtering" ]]; then
		mkdir -p $dir/1.Filtering; ln -sf $assembly $dir/1.Filtering/01.assembly.fa
	fi
	time1=`date +%s`
	echo -e "\n\nSTEP 2: MegaBLAST starting..."; echo -e "Current date/time: $(date)\n"
	mkdir -p $dir/2.MegaBLAST; cd $dir/2.MegaBLAST; rm -rf *
	if [ $cores -ge $blast_blocks_size ]; then
		echo -e "\nProcessing in the MegaBLAST simultaneously the individual contigs in blocks of at most $blast_blocks_size elements, please manually change the variable '-B' in the ILRA.sh main script if required, for example because less cores available or running into memory issues...\n"
		cat $dir/1.Filtering/01.assembly.fa | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename }'
		arr=($(ls -lS | sort -k 5 -n | awk '{print $9}' | awk 'NF' | egrep .fa$ | tac))
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" parallel --verbose --joblog megablast_parallel_log_out_2.txt -j $blast_blocks_size megablast -W 40 -F F -a $((cores / blast_blocks_size *2)) -m 8 -e 1e-80 -d $dir/1.Filtering/01.assembly.fa -i {} -o comp.self1.{}.blast ::: ${arr[@]} &> megablast_parallel_log_out.txt
		awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' megablast_parallel_log_out_2.txt > tmp && mv tmp megablast_parallel_log_out_2.txt
		# Control if failed due to any reason (likely too much RAM...)
		if [ "$(awk '{ print $9 }' megablast_parallel_log_out_2.txt | grep -c "9")" -gt 0 ]; then
			echo -e "\nRunning the megablast of the following samples again giving more resources:"
			echo $(awk '$9 == 9 {print}' megablast_parallel_log_out_2.txt | awk '{ print $24 }')
			arr=($(awk '$9 == 9 {print}' megablast_parallel_log_out_2.txt | awk '{ print $24 }')); length_arr=${#arr[@]}; export length_arr
			cores_split=$((cores / blast_blocks_size_large))
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" parallel --verbose --joblog megablast_parallel_log_out_2_2.txt -j $blast_blocks_size_large megablast -W 40 -F F -a $cores_split -m 8 -e 1e-80 -d $dir/1.Filtering/01.assembly.fa -i {} -o comp.self1.{}.blast ::: ${arr[@]} &> megablast_parallel_log_out2.txt
			awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' megablast_parallel_log_out_2_2.txt > tmp && mv tmp megablast_parallel_log_out_2_2.txt
		fi
		if [ "$(awk '{ print $9 }' megablast_parallel_log_out_2_2.txt | grep -c "9")" -gt 0 ]; then
			echo -e "\nThe MegaBLAST execution failed for some of the contigs, likely due to excessive RAM usage, please double check manually the log megablast_parallel_log_out_2_2.txt or look for support. Exiting for now...\n"
			exit 1
		fi
		cat *.fa.blast | awk -v overlap_identity="$overlap_identity" -v overlap_length="$overlap_length" '$3>overlap_identity && $4>overlap_length && $1 != $2' > comp.self1.blast; rm *.fa *.fa.blast
	else
		megablast -W 40 -F F -a $cores -m 8 -e 1e-80 -d $dir/1.Filtering/01.assembly.fa -i $dir/1.Filtering/01.assembly.fa | awk -v overlap_identity="$overlap_identity" -v overlap_length="$overlap_length" '$3>overlap_identity && $4>overlap_length && $1 != $2' > comp.self1.blast
	fi
	ILRA.addLengthBlast.pl $dir/1.Filtering/01.assembly.fa $dir/1.Filtering/01.assembly.fa comp.self1.blast &> /dev/null

	#### 2a. Delete contained contigs
	# We want the query to be always the smaller one
	echo -e "\n\nSTEP 2a: Delete contained contigs starting..."; echo -e "Current date/time: $(date)\n"
	awk -v overlap_identity="$overlap_identity" -v overlap_length="$overlap_length" '$3>overlap_identity && $4>overlap_length && $13 < $14' comp.self1.blast.length | ILRA.getOverlap.pl > 02.ListContained.txt
	cat 02.ListContained.txt | awk -v overlap_fraction="$overlap_fraction" '$4>overlap_fraction' | cut -f 1 > List.Contained.fofn
	echo -e "\n### Excluded contigs that are contained in others: (MegaBLAST/ILRA.getOverlap.pl/ILRA.deleteContigs.pl)" >> ../Excluded.contigs.fofn
	cat List.Contained.fofn >> ../Excluded.contigs.fofn
	# Now delete the contigs from the MegaBLAST first output, also filter for just hits > Xkb > XX%
	echo "Filtering the MegaBLAST output for hits > $overlap_length pb and >$overlap_identity% identity. Please change the parameters if needed"
	ILRA.deleteContigs.pl List.Contained.fofn $dir/1.Filtering/01.assembly.fa 02.assembly.fa
	cat comp.self1.blast.length | awk -v overlap_identity="$overlap_identity" -v overlap_length="$overlap_length" '$3> overlap_identity && $4 > overlap_length' | ILRA.deleteEntryinBlast.pl List.Contained.fofn > Blast.merge.blast.length
debug="all"
fi
if [[ $debug == "all" || $debug == "step2b" ]]; then
	#### 2b. Find overlaps
	# Get a bam files to do the merge
	echo -e "\n\nSTEP 2b: Find overlaps starting..."; echo -e "Current date/time: $(date)\n"
	cd $dir/2.MegaBLAST
	if [ -f $illuminaReads\_1.fastq.gz ]; then
		echo -e "SMALT parameters are: k-mer size=20, step-size=3, insert size range="$InsertsizeRange". Please check SMALT help and change manually within the pipeline (section 2b) these parameters if needed"
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" ILRA.runSMALT_ver2.sh 02.assembly.fa 20 3 $illuminaReads\_1.fastq.gz $illuminaReads\_2.fastq.gz first $InsertsizeRange $cores 9 &> ILRA.runSMALT_ver2.sh_log_out.txt
		echo -e "Check out the log of ILRA.runSMALT_ver2.sh in the file ILRA.runSMALT_ver2.sh_log_out.txt"
	# Find overlaps
		ILRA.findoverlaps_ver3.pl Blast.merge.blast.length first.bam 02.assembly.fa OUT $merging_coverage_threshold $merging_coverage_deviation 1> ILRA.findoverlaps_ver3.pl_log_out.txt 2> ILRA.findoverlaps_ver3.pl_log_out_perl_warnings_errors.txt
		echo -e "\nCheck out the log of ILRA.findoverlaps_ver3.pl in the files log_OUT.txt, results_OUT.txt, ILRA.findoverlaps_ver3.pl_log_out.txt, and ILRA.findoverlaps_ver3.pl_log_out_perl_warnings_errors.txt"
	# Save the contigs
		echo -e "\n### New contigs merging overlapping contigs(ILRA.findoverlaps_ver3.pl):" >> ../Excluded.contigs.fofn
		if [ $(cat log_OUT.txt | grep ^sequences: | awk '{ gsub("sequences: " , ""); print }') -eq 0 ]; then
			echo "Nothing filtered by ILRA.findoverlaps_ver3.pl"
		else
			cat results_OUT.txt >> ../Excluded.contigs.fofn
		fi
	# Cleaned assembly
		if [ $contigs_illumina_filter == "yes" ]; then
			echo -e "\nRemoving overlapping contigs based on Illumina short reads...\n"
			echo -e "\n### Contigs not covered: (ILRA.findoverlaps_ver3.pl)" >> ../Excluded.contigs.fofn
			cat notcovered_OUT.fasta | awk 'BEGIN {RS = ">"; FS = "\n"; ORS = ""} $2 {print ">"$0}' | grep ">" >> ../Excluded.contigs.fofn
			cp mergedseq_OUT.fasta 03.assembly.fa
		else
			echo -e "\nSkipping removing overlapping contigs based on Illumina short reads because of uneven coverage...\n"
			cat mergedseq_OUT.fasta notcovered_OUT.fasta > 03.assembly.fa
		fi
		rm $(find $dir/2.MegaBLAST/ -name "first.bam*") $(find $dir/2.MegaBLAST/ -name "Genome*") # Cleaning
	else
	# Bypass if Illumina reads not provided
		echo -e "Illumina reads NOT PROVIDED and no filtering by ILRA.findoverlaps_ver3.pl\n"
		cp 02.assembly.fa 03.assembly.fa
	fi
	echo -e "\nAfter this step:"; assembly-stats 02.assembly.fa | head -n 2; assembly-stats 03.assembly.fa | head -n 2
	echo -e "\nSTEP 2: DONE"; echo -e "Current date/time: $(date)"
	time2=`date +%s`; echo -e "STEP 2 time (secs): $((time2-time1))"; echo -e "STEP 2 time (hours): $(echo "scale=2; $((time2-time1))/3600" | bc -l)"
debug="all"
fi

#### 3. ABACAS2
if [[ $debug == "all" || $debug == "step3" ]]; then
	test_installation_dependencies.sh
	if [[ ! -d "$dir/2.MegaBLAST" ]]; then
		mkdir -p $dir/2.MegaBLAST; ln -sf $assembly $dir/2.MegaBLAST/03.assembly.fa
	fi
	time1=`date +%s`
	echo -e "\n\nSTEP 3: ABACAS2 starting..."; echo -e "Current date/time: $(date)\n"
	mkdir -p $dir/3.ABACAS2; cd $dir/3.ABACAS2; rm -rf *
	if [ "$doAbacas2" -eq 1 ]; then
		ABA_CHECK_OVERLAP=0; export ABA_CHECK_OVERLAP; ABA_COMPARISON=nucmer; export ABA_COMPARISON; ABA_SPLIT_PARTS=$abacas2_split; export ABA_SPLIT_PARTS; Min_Alignment_Length=1000; Identity_Cutoff=98
		echo -e "ABACAS2 parameters are:\nABA_CHECK_OVERLAP=0\nABA_COMPARISON=nucmer\nMin_Alignment_Length=1000\nIdentity_Cutoff=98\nBlast=$abacas2_blast\nsplit_parts=$ABA_SPLIT_PARTS\nPlease check ABACAS2 help and change manually within the pipeline (section 3) these parameters if needed"
		abacas2.nonparallel.sh $reference $dir/2.MegaBLAST/03.assembly.fa $cores $Min_Alignment_Length $Identity_Cutoff $abacas2_blast 1> abacas_log_out.txt 2> abacas_log_out_warnings_errors.txt
		rm -rf $(find $dir/3.ABACAS2/ -name "Ref.*") $(find $dir/3.ABACAS2/ -name "Res.*") $(find $dir/3.ABACAS2/ -name "*.coords") $dir/3.ABACAS2/Reference # Cleaning		
		echo -e "\nCheck out the log of abacas2.nonparallel.sh in the files abacas_log_out.txt and abacas_log_out_warnings_errors.txt"
	# Break  and delete N's
		fastaq trim_Ns_at_end Genome.abacas.fasta 03b.assembly.fa
	# Save the contigs
		# echo -e "\n### Overlapping contigs_2 (ABACAS2)" >> ../Excluded.contigs.fofn # Uncomment if ABA_CHECK_OVERLAP=1
		# zcat *.overlap_report.gz >> ../Excluded.contigs.fofn # Uncomment if ABA_CHECK_OVERLAP=1
		echo -e "\n### Final sequences not mapped to the reference: (ABACAS2, these are contained into the bin and are listed here but not removed)" >> ../Excluded.contigs.fofn
		grep ">" Res.abacasBin.fna >> ../Excluded.contigs.fofn
		echo -e "\n### Final sequences mapped to the reference: (ABACAS2, these are reordered and merged into new contigs corresponding to the reference)" >> ../Excluded.contigs.fofn
		for i in $(ls | grep .contigs.gff); do
			echo "# "$(echo $i | sed 's/^[^.]*.\([^.]*\)..*/\1/')":" >> ../Excluded.contigs.fofn
			cat $i | grep contig | awk '{ print $9 }' | sed 's/^[^"]*"\([^"]*\)".*/\1/' | sed 's/.*=//' >> ../Excluded.contigs.fofn
		done
	else
	# Bypass if necessary:
		echo -e "Reference genome NOT PROVIDED and ABACAS2 not executed. STEP 3 for reordering and renaming is skipped\n"
		ln -sf $dir/2.MegaBLAST/03.assembly.fa 03b.assembly.fa
	fi
	echo -e "\nAfter this step:"; assembly-stats 03b.assembly.fa | head -n 2
	echo -e "\nSTEP 3: DONE"; echo -e "Current date/time: $(date)"
	time2=`date +%s`; echo -e "STEP 3 time (secs): $((time2-time1))"; echo -e "STEP 3 time (hours): $(echo "scale=2; $((time2-time1))/3600" | bc -l)"
debug="all"
fi

#### 4. Correction via Illumina short reads
#### iCORN2 or Pilon:
if [[ $debug == "all" || $debug == "step4" || $debug == "step4i" ]]; then
	test_installation_dependencies.sh
	if [[ ! -d "$dir/3.ABACAS2" ]]; then
		mkdir -p $dir/3.ABACAS2; ln -sf $assembly $dir/3.ABACAS2/03b.assembly.fa
	fi
	if [[ $pilon == "no" ]]; then
		time1=`date +%s`
		echo -e "\n\nSTEP 4: iCORN2 starting..."; echo -e "Current date/time: $(date)\n"
		mkdir -p $dir/4.iCORN2; cd $dir/4.iCORN2
		if [[ $debug != "step4i" ]]; then # Deal with whether to remove everything because we are resuming step 4 from the beginning or only some iterations
			rm -rf *
		fi
		if [ -f $illuminaReads\_1.fastq.gz ]; then
			iCORN2_fragmentSize=500
			echo "The iCORN2 fragment size used is iCORN2_fragmentSize=$iCORN2_fragmentSize. Please check iCORN2 help and change manually within the pipeline (section 4) if needed"
			echo -e "Check out the log of icorn2.serial_bowtie2.sh in the files icorn2.serial_bowtie2.sh_log_out.txt and ../7.Stats/07.iCORN2.final_corrections.results.txt"
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" icorn2.serial_bowtie2.sh $illuminaReads $iCORN2_fragmentSize $dir/3.ABACAS2/03b.assembly.fa 1 $number_iterations_icorn $blocks_size $cores $java_memory $parts_icorn2_split $low_mem $low_spa &> icorn2.serial_bowtie2.sh_log_out.txt
			for i in $(find $dir -type d -name "ICORN2_*"); do cd $i && tar cf out_all.tar *.out && rm -rf *.out; done; cd $dir/4.iCORN2 # Cleaning
			if [ "$(find . -name "larger_contigs2.txt" | wc -l)" -gt 0 ]; then
				echo -e "\nPlease note iCORN2 is likely going to fail if any of the sequences being processed is larger than or around 60Mb... The program tried to handle this and divide the input sequences accordingly but failed...\n"
				echo -e "\n\nThe size of at least one individual contig or chromosome within the input sequence is ~60Mb... iCORN2 used Pilon to correct some sequences. If allowed to continue, there would be segmentation errors later on in SNP-o-matic\n\n"
			fi
	# Cleaned assembly:
			ln -sf ICORN2.03b.assembly.fa.[$(expr $number_iterations_icorn + 1)] 04.assembly.fa
			if [ "$(find $dir/4.iCORN2 -name "*pilon.fa*" | wc -l)" -eq 0 ]; then
	# Summary of iCORN2 results:
				mkdir -p $dir/7.Stats; icorn2.collectResults.pl $PWD > ../7.Stats/07.iCORN2.final_corrections.results.txt
				echo -e "### Total SNP: " >> ../7.Stats/07.iCORN2.final_corrections.results.txt; cat ../7.Stats/07.iCORN2.final_corrections.results.txt | awk '{sum+=$6;} END{print sum;}' >> ../7.Stats/07.iCORN2.final_corrections.results.txt
				echo -e "### Total INS: " >> ../7.Stats/07.iCORN2.final_corrections.results.txt; cat ../7.Stats/07.iCORN2.final_corrections.results.txt | awk '{sum+=$7;} END{print sum;}' >> ../7.Stats/07.iCORN2.final_corrections.results.txt
				echo -e "### Total DEL: " >> ../7.Stats/07.iCORN2.final_corrections.results.txt; cat ../7.Stats/07.iCORN2.final_corrections.results.txt | awk '{sum+=$8;} END{print sum;}' >> ../7.Stats/07.iCORN2.final_corrections.results.txt
				echo -e "### Total HETERO: " >> ../7.Stats/07.iCORN2.final_corrections.results.txt; cat ../7.Stats/07.iCORN2.final_corrections.results.txt | awk '{sum+=$9;} END{print sum;}' >> ../7.Stats/07.iCORN2.final_corrections.results.txt
				echo -e "\nA preview of the correction by iCORN2 is:"
				cat ../7.Stats/07.iCORN2.final_corrections.results.txt
			else
				echo -e "\nSome sequences were corrected by iCORN2 with Pilon instead of SNP-o-matic, please check manually if required, but global summary of the corrections is:"
				egrep "^Corrected SNP-o-matic|^Corrected Pilon" icorn2.serial_bowtie2.sh_log_out.txt
			fi
			# Cleaning:
			for ((i=1;$i<=$number_iterations_icorn;i++)); do
				sed -s -e $'$a\\\n' $(ls | egrep .$i.o$) > log.$i.out
				rm -rf $(ls | egrep .$i.o$) tmp_dir
				cd $dir/4.iCORN2/ICORN2_$i; ls -l &> list_files_prev_cleaning.txt; rm *.bam *.bai *.fa *.fai *.dict 
				cd $dir/4.iCORN2/
			done			
		else
	# Bypass if Illumina reads not provided
			echo -e "Illumina reads NOT PROVIDED and iCORN2 is not executed. STEP 4 for correction using short reads is skipped\n"
			ln -sf $dir/3.ABACAS2/03b.assembly.fa 04.assembly.fa
		fi
	elif [[ $pilon == "yes" ]]; then
		type pilon >/dev/null 2>&1 || { echo >&2 "I require to execute 'pilon' but it's not installed or available in the PATH, please install or make a wrapper of the java call so ILRA can execute 'pilon'. Aborting..."; exit 1; }
		time1=`date +%s`
		echo -e "\n\nSTEP 4: Pilon starting..."; echo -e "Current date/time: $(date)\n"
		mkdir -p $dir/4.pilon; cd $dir/4.pilon; rm -rf *
		if [ -f $illuminaReads\_1.fastq.gz ]; then
			for i in $(seq 1 1 $number_iterations_icorn); do
				if [ "$i" -eq 1 ]; then
					mkdir -p $dir/4.pilon/iter_1
					\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" bowtie2-build --threads $cores $dir/3.ABACAS2/03b.assembly.fa genome_iter1 &> bowtie_log_out.txt 2> bowtie_log_out.txt
					\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" bowtie2 -t -x genome_iter1 -p $cores -X 1200 --very-sensitive -N 1 -L 31 --rdg 5,2 -1 $illuminaReads\_1.fastq.gz -2 $illuminaReads\_2.fastq.gz 2>> bowtie_log_out.txt | samtools sort -l 9 -m 2G -@ $cores --write-index -o "ill_reads1.bam##idx##ill_reads1.bam.bai" &>> bowtie_log_out.txt
					\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" pilon --genome $dir/3.ABACAS2/03b.assembly.fa --bam ill_reads1.bam --output genome_pilon1 --outdir $dir/4.pilon/iter_1 --changes --vcf --tracks &>> pilon_log_out.txt
				else
					mkdir -p $dir/4.pilon/iter_$i
					\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" bowtie2-build --threads $cores $dir/4.pilon/iter_[$(expr $i - 1)]/genome_pilon[$(expr $i - 1)].fasta genome_iter$i &> bowtie_log_out.txt
					\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" bowtie2 -t -x genome_iter$i -p $cores -X 1200 --very-sensitive -N 1 -L 31 --rdg 5,2 -1 $illuminaReads\_1.fastq.gz -2 $illuminaReads\_2.fastq.gz 2>> bowtie_log_out.txt | samtools sort -l 9 -m 2G -@ $cores --write-index -o "ill_reads$i.bam##idx##ill_reads$i.bam.bai" &>> bowtie_log_out.txt
					\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" pilon --genome $dir/4.pilon/iter_[$(expr $i - 1)]/genome_pilon[$(expr $i - 1)].fasta --bam ill_reads$i.bam --output genome_pilon$i --outdir $dir/4.pilon/iter_$i --changes --vcf --tracks &>> pilon_log_out.txt
				fi
			done
	# Cleaned assembly:
			echo -e "Check out the log of pilon in the file pilon_log_out.txt and *.changes in the subfolders. The numbers of changes by pilon in the $number_iterations_icorn iterations are:"
			wc -l $(find . -name "*.changes")
			ln -sf $dir/4.pilon/iter_$number_iterations_icorn/genome_pilon$number_iterations_icorn.fasta 04.assembly.fa
		else
	# Bypass if Illumina reads not provided
			echo -e "Illumina reads NOT PROVIDED and pilon is not executed. STEP 4 for correction using short reads is skipped\n"
			ln -sf $dir/3.ABACAS2/03b.assembly.fa 04.assembly.fa
		fi
	fi
echo -e "\nAfter this step:"; assembly-stats 04.assembly.fa | head -n 2
echo -e "\nSTEP 4: DONE"; echo -e "Current date/time: $(date)"
time2=`date +%s`; echo -e "STEP 4 time (secs): $((time2-time1))"; echo -e "STEP 4 time (hours): $(echo "scale=2; $((time2-time1))/3600" | bc -l)"
debug="all"
fi

#### 5. Circlator for organelles
if [[ $debug == "all" || $debug == "step5" ]]; then
	test_installation_dependencies.sh
	if [[ ! -d "$dir/4.iCORN2" ]] && [[ ! -d "$dir/4.pilon" ]]; then
		mkdir -p $dir/4.iCORN2_pilon; ln -sf $assembly $dir/4.iCORN2_pilon/04.assembly.fa
	fi
	time1=`date +%s`
	if [ ! -z "$correctedReads" ]; then
		echo -e "\n\nSTEP 5: Circlator starting..."; echo -e "Current date/time: $(date)\n"
		mkdir -p $dir/5.Circlator; cd $dir/5.Circlator; rm -rf *
		if grep -q -E "$seqs_circl" $(find $dir -name "04.assembly.fa"); then
	# Map the corrected reads (minimap2 can be used and it's faster, but meryl+winnowmap is currently recommended to map long reads against repetitive sequences)
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" meryl count k=15 threads=$cores output merylDB $(find $dir -name "04.assembly.fa") &> meryl_count_log_out.txt
			meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt 2> meryl_print_log_out.txt
			if [ $long_reads_technology == "pacbio" ]; then
				#minimap2 -x map-pb -H -t $cores -d minimap2_index_pacbio.mmi $(find $dir -name "04.assembly.fa") &> mapping_corrected_reads_log_out.txt
				#minimap2 -x map-pb -t $cores -a minimap2_index_pacbio.mmi $correctedReads > Mapped.corrected.04.sam 2>> mapping_corrected_reads_log_out.txt
				winnowmap -W repetitive_k15.txt -ax map-pb -t $cores $(find $dir -name "04.assembly.fa") $correctedReads 2>> mapping_corrected_reads_log_out.txt | samtools view -@ $cores -h -F 260 -F 2048 -o Mapped.corrected.04.sam
			else
				#minimap2 -x map-ont -t $cores -d minimap2_index_nanopore.mmi $(find $dir -name "04.assembly.fa") &> mapping_corrected_reads_log_out.txt
				#minimap2 -x map-ont -t $cores -a minimap2_index_nanopore.mmi $correctedReads > Mapped.corrected.04.sam 2>> mapping_corrected_reads_log_out.txt
				winnowmap -W repetitive_k15.txt -ax map-ont -t $cores $(find $dir -name "04.assembly.fa") $correctedReads 2>> mapping_corrected_reads_log_out.txt | samtools view -@ $cores -h -F 260 -F 2048 -o Mapped.corrected.04.sam
			fi
	# Circlator:
			seq_ids=$(awk '{print $3}' Mapped.corrected.04.sam | grep -E "$seqs_circl" | tail -n +2 | sort | uniq)
			for i in $seq_ids; do
				cat Mapped.corrected.04.sam | grep $i | awk '{ print ">"$1"\n"$10 }' >> ForCirc.reads.fasta
			done
			for i in $seq_ids; do
				samtools faidx $(find $dir -name "04.assembly.fa") $i >> ForCirc.Ref.fasta
			done
			awk 'BEGIN {RS = ">"; FS = "\n"; ORS = ""} {if ($2) print ">"$0}' ForCirc.reads.fasta > ForCirc.reads_2.fasta
			awk 'BEGIN {RS = ">"; FS = "\n"; ORS = ""} {if ($2) print ">"$0}' ForCirc.Ref.fasta > ForCirc.Ref_2.fasta
			echo -e "Check out the log of the mapping of the corrected reads and circlator in the files mapping_corrected_reads_log_out.txt and circlator_log_out.txt"
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" circlator all ForCirc.Ref_2.fasta ForCirc.reads_2.fasta Out.Circ --threads $cores &> circlator_log_out.txt
			rm -rf $dir/5.Circlator/merylDB # Cleaning
	# Delete the plastids/organelles/circular sequences from the current assembly version (04.assembly.fa)
			for i in $seq_ids; do
				echo $i >> List.circular_sequences.fofn
			done
			ILRA.deleteContigs.pl List.circular_sequences.fofn $(find $dir -name "04.assembly.fa") 05.assembly.fa
			if [ ! -s $dir/5.Circlator/Out.Circ/06.fixstart.fasta ]; then
	# Bypass if circlator failed and didn't end
				echo -e "PLEASE be aware that Circlator HAS FAILED. Check its log and output for futher details. This may be due to some problems with the provided reads, or not enough/even coverage for the required contigs. Bypassing..."
				ln -sf $(find $dir -name "04.assembly.fa") 05.assembly.fa
			else
				cat Out.Circ/06.fixstart.fasta >> 05.assembly.fa;
			fi
		else
	# Bypass if the contigs names provided are not found
			ln -sf $(find $dir -name "04.assembly.fa") 05.assembly.fa
			echo -e "The sequence identifiers provided for circularization were not found in the contig names. Circlator NOT EXECUTED. STEP 5 for circularization is skipped"
		fi
		debug="all"
	else
		echo -e "\nCORRECTED READS NOT PROVIDED, Circlator NOT EXECUTED. STEP 5 for circularization is skipped"
		mkdir -p $dir/5.Circlator; cd $dir/5.Circlator
		ln -sf $(find $dir -name "04.assembly.fa") 05.assembly.fa
	fi
echo -e "\nAfter this step:"; assembly-stats 05.assembly.fa | head -n 2
echo -e "\nSTEP 5: DONE"; echo -e "Current date/time: $(date)"
time2=`date +%s`; echo -e "STEP 5 time (secs): $((time2-time1))"; echo -e "STEP 5 time (hours): $(echo "scale=2; $((time2-time1))/3600" | bc -l)"
fi

#### 6. Decontamination/taxonomic classification/final masking and filtering for databases upload, rename sequences
if [[ $debug == "all" || $debug == "step6" ]]; then
	test_installation_dependencies.sh
	if [[ ! -d "$dir/5.Circlator" ]]; then
		mkdir -p $dir/5.Circlator; ln -sf $assembly $dir/5.Circlator/05.assembly.fa
	fi
	if [[ $mode == "taxon" || $mode == "both" ]]; then
		time1=`date +%s`
		echo -e "\n\nSTEP 6: Kraken2 decontamination starting..."; echo -e "Current date/time: $(date)\n"
		mkdir -p $dir/6.Decontamination; cd $dir/6.Decontamination; rm -rf *		
		if [[ $kraken2_fast == "yes" ]]; then
			echo -e "\nPreparing database $databases/$kraken2_databases for fast access in RAM...\n"
			cp -ru $databases/$kraken2_databases /dev/shm/
			echo -e "\nExecuting kraken2...\n"
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" kraken2 --db /dev/shm/$kraken2_databases --threads $cores --memory-mapping --classified-out classification.txt --unclassified-out classification_unknwn.txt --report report.txt --output kraken2_output.txt --use-names $dir/5.Circlator/05.assembly.fa 1> kraken2_log_out.txt 2> kraken2_log_warnings_errors.txt
		elif [[ $kraken2_fast == "no" ]]; then
			echo -e "\nExecuting kraken2...\n"
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" kraken2 --db $databases/$kraken2_databases --threads $cores --classified-out classification.txt --unclassified-out classification_unknwn.txt --report report.txt --output kraken2_output.txt --use-names $dir/5.Circlator/05.assembly.fa 1> kraken2_log_out.txt 2> kraken2_log_warnings_errors.txt
		fi
		echo -e "Please check the files report.txt, kraken2_log_out.txt and kraken2_log_out_warnings_errors.txt"
		if [ -s report.txt ]; then
			echo -e "\nLog of kraken2:"; echo -e "%_contigs_covered\t#_contigs_covered\t#_contigs_directly_assigned\tRank_code\tTaxon_id\tScientific_name"; cat report.txt
			cat report.txt > report_final.txt; echo -e "\n\nNumber of classified contigs at the genus level: $(cat report.txt | awk '$4 == "G" {print $2"\t"$5}' | awk '{s+=$1}END{print s}')" >> report_final.txt
			echo -e "\nTaxonomy IDs at the genus level assigned to the contigs:" >> report_final.txt; echo -e "#Contig\tTaxID\n$(awk '$4 == "G" {print $2"\t"$5}' report.txt)\n" >> report_final.txt
	# Extract contigs classified as different organisms
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" rcf -n $databases/taxdump -k kraken2_output.txt -o recentrifuge_contamination_report.html -e CSV &> rcf_log_out.txt # Add --sequential if problems with multithreading
			re='^[0-9]+$' # Deal with the user providing taxon id or taxon names
			if [[ $taxonid =~ $re ]]; then
				taxon_name=$(taxonkit list --ids $taxonid -n -r --data-dir $databases/taxdump | grep $taxonid)
			else
				taxonid=$(echo $taxonid | taxonkit name2taxid --data-dir $databases/taxdump | head -1 | cut -f2)
				taxon_name=$(taxonkit list --ids $taxonid -n -r --data-dir $databases/taxdump | grep $taxonid)
			fi
			echo -e "\nOrganisms provided: $taxon_name"
			echo -e "\nIf not correct, please provide it explicitely using the -T argument and rerun. Kraken2 output will be filtered to retain that taxa and below"
			echo -e "\nCheck out the logs in the files rcf_log_out.txt and extract_kraken2_log_out.txt"
			echo -e "\nExtracting reads...\n"
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" extract_kraken_reads.py -k kraken2_output.txt -U $dir/5.Circlator/05.assembly.fa -o 06.assembly.fa -t $taxonid -r report.txt --include-children &> extract_kraken2_log_out.txt
	# Include the ones that may be unclassified:
			#if [ -s classification_unknwn.txt ]; then
			#	cat classification_unknwn.txt >> 06.assembly.fa
			#fi
	# Save contigs
			echo -e "\n### Excluded contigs that are not recognized by kraken2 as the species of interest  or that are unclassified:" >> ../Excluded.contigs.fofn
			comm -23 <(cat $dir/5.Circlator/05.assembly.fa | grep ">" | sort) <(cat 06.assembly.fa | grep ">" | sort) >> ../Excluded.contigs.fofn
	# Clean
			pigz --best -p $cores classification.txt
		fi
		if [ -s 06.assembly.fa ]; then
			echo -e "kraken2 and extraction seem to have run fine. Please check the logs, the Excluded.contigs.fofn and the report.txt files to assess the contigs that corresponded to contamination and that have been excluded by ILRA"
		else
	# Bypass if kraken2 failed and didn't end
			echo -e "\nPLEASE BE AWARE that decontamination based on taxonomic classification HAS FAILED, kraken2 may have been killed due to RAM usage, it may have failed due to other reasons... Please check the logs. For now bypassing step 6 so ILRA can continue..."
			ln -sf $dir/5.Circlator/05.assembly.fa 06.assembly.fa
		fi
	# Renaming and making single-line fasta with length in the headers
		if [ -z "$reference" ]; then
			cat 06.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_".$n } else { print }' | ILRA.fasta2singleLine.pl - | awk '/^>/ { if (name) {printf("%s_%d\n%s", name, len, seq)} name=$0; seq=""; len = 0; next} NF > 0 {seq = seq $0 "\n"; len += length()} END { if (name) {printf("%s_%d\n%s", name, len, seq)} }' > ../$name.ILRA.fasta
		else
			cat 06.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_with_ref_".$n } else { print }' | ILRA.fasta2singleLine.pl - | awk '/^>/ { if (name) {printf("%s_%d\n%s", name, len, seq)} name=$0; seq=""; len = 0; next} NF > 0 {seq = seq $0 "\n"; len += length()} END { if (name) {printf("%s_%d\n%s", name, len, seq)} }' > ../$name.ILRA.fasta
		fi
		echo -e "\nAfter this step:"; assembly-stats 06.assembly.fa | head -n 2
		echo -e "\nSTEP 6 KRAKEN2: DONE"; echo -e "Current date/time: $(date)"
		time2=`date +%s`; echo -e "STEP 6 kraken2 time (secs): $((time2-time1))"; echo -e "STEP 6 kraken2 time (hours): $(echo "scale=2; $((time2-time1))/3600" | bc -l)"
	fi
	if [[ $mode == "blast" || $mode == "both" ]]; then
	# Final filtering, masking and reformatting to conform the requirements of DDBJ/ENA/Genbank
		mkdir -p $dir/6.Decontamination/blast_filtering_for_Genbank_upload; cd $dir/6.Decontamination/blast_filtering_for_Genbank_upload; rm -rf *
		time1=`date +%s`
		echo -e "\n\nSTEP 6: Blast decontamination starting..."; echo -e "Current date/time: $(date)\n"
		if [ $mode == "blast" ]; then
	# Renaming and making single-line fasta with length in the headers
			if [ -z "$reference" ]; then
				cat $dir/5.Circlator/05.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_".$n } else { print }' | ILRA.fasta2singleLine.pl - | awk '/^>/ { if (name) {printf("%s_%d\n%s", name, len, seq)} name=$0; seq=""; len = 0; next} NF > 0 {seq = seq $0 "\n"; len += length()} END { if (name) {printf("%s_%d\n%s", name, len, seq)} }' > ../../$name.ILRA.fasta
			else
				cat $dir/5.Circlator/05.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_with_ref_".$n } else { print }' | ILRA.fasta2singleLine.pl - | awk '/^>/ { if (name) {printf("%s_%d\n%s", name, len, seq)} name=$0; seq=""; len = 0; next} NF > 0 {seq = seq $0 "\n"; len += length()} END { if (name) {printf("%s_%d\n%s", name, len, seq)} }' > ../../$name.ILRA.fasta
			fi
		fi
	# BLAST to screen the submitted sequences against: (this and all the blast commands or thresholds are following Genbank requirements, the NCBI's Foreign Contamination Screen)
	cat ../../$name.ILRA.fasta | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename }'
	arr=($(ls -lS | awk '{print $9}' | awk 'NF' | egrep .fa$ | grep -v ref.fa))
	# 1. Common contaminants database that contains vector sequences, bacterial insertion sequences, E. coli and phage genomes...
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" parallel -q --verbose --joblog blast_contam_euks_log_out_2.txt -j $blocks_size blastn -query {} -db $databases/contam_in_euks.fa -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads $((cores / blocks_size)) -out {}.contam_in_euks_genbank.out ::: ${arr[@]} &> blast_contam_euks_log_out.txt
		awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' blast_contam_euks_log_out_2.txt > tmp && mv tmp blast_contam_euks_log_out_2.txt
		cat *.fa.contam_in_euks_genbank.out | egrep -v "^#" | awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' > contam_in_euks_genbank.out; rm *.fa.contam_in_euks_genbank.out
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" parallel -q --verbose --joblog blast_contam_proks_log_out_2.txt -j $blocks_size blastn -query {} -db $databases/contam_in_prok.fa -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads $((cores / blocks_size)) -out {}.contam_in_proks_genbank.out ::: ${arr[@]} &> blast_contam_proks_log_out.txt
		awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' blast_contam_proks_log_out_2.txt > tmp && mv tmp blast_contam_proks_log_out_2.txt
		cat *.fa.contam_in_proks_genbank.out | egrep -v "^#" | awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' > contam_in_proks_genbank.out; rm *.fa.contam_in_proks_genbank.out
	# 2. A database of adaptors linkers and primers
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" parallel --verbose --joblog vecscreen_contam_euks_log_out_2.txt -j $cores vecscreen -d $databases/adaptors_for_screening_euks.fa -f3 -i {} -o {}.vecscreen_in_euks_genbank ::: ${arr[@]} &> vecscreen_contam_euks_log_out.txt
		awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' vecscreen_contam_euks_log_out_2.txt > tmp && mv tmp vecscreen_contam_euks_log_out_2.txt
		cat *.fa.vecscreen_in_euks_genbank > vecscreen_in_euks_genbank; rm *.fa.vecscreen_in_euks_genbank
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" parallel --verbose --joblog vecscreen_contam_proks_log_out_2.txt -j $cores vecscreen -d $databases/adaptors_for_screening_proks.fa -f3 -i {} -o {}.vecscreen_in_proks_genbank ::: ${arr[@]} &> vecscreen_contam_proks_log_out.txt
		awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' vecscreen_contam_proks_log_out_2.txt > tmp && mv tmp vecscreen_contam_proks_log_out_2.txt
		cat *.fa.vecscreen_in_proks_genbank > vecscreen_in_proks_genbank; rm *.fa.vecscreen_in_proks_genbank
		VSlistTo1HitPerLine.sh suspect=0 weak=0 vecscreen_in_euks_genbank > vecscreen_in_euks_genbank.out
		VSlistTo1HitPerLine.sh suspect=0 weak=0 vecscreen_in_proks_genbank > vecscreen_in_proks_genbank.out
	# 3. A database of mitochondrial genomes
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" parallel --verbose --joblog blast_mito_log_out_2.txt -j $blocks_size blastn -query {} -db $databases/mito -out {}.mito_sequences -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 98.6 -soft_masking true -outfmt 7 -num_threads $((cores / blocks_size)) ::: ${arr[@]} &> blast_mito_log_out.txt
		awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' blast_mito_log_out_2.txt > tmp && mv tmp blast_mito_log_out_2.txt
		cat *.fa.mito_sequences | egrep -v "^#" > mito_sequences; rm *.fa.mito_sequences
	# Deal with several hits: (from the BLAST output get the largest alignments and then the largest % of identity)
		if [ "$(cat mito_sequences | wc -l)" -gt 0 ]; then
			echo "#Hits:" > mito_sequences.out
			grep "Fields: " mito_sequences | uniq | sed -e 's/,/\t/g' | sed 's/.*: //' >> mito_sequences.out
			cat mito_sequences >> mito_sequences.out
			echo -e "\n#Info sequences of the hits:" >> mito_sequences.out
			for i in $(cat mito_sequences | cut -f2 | sort | uniq); do
				blastdbcmd -db $databases/mito -entry $i | grep ">" >> mito_sequences.out
			done
			echo -e "\n#Filtering hits - largest alignments:" >> mito_sequences.out
			for i in $(cat mito_sequences | cut -f1 | sort | uniq); do
				grep -e "^$i" mito_sequences | awk -v aln_length=$(grep -e "^$i" mito_sequences | sort -n -k4 -r | head -n 1 | cut -f 4) '$4 == aln_length' >> mito_sequences.out
			done
			seq_mit=$(awk '/largest/,0' mito_sequences.out | awk 'NR!=1' | sort -n -k4 -r | head -n 1 | cut -f 1)
			echo -e "\n#Contig chosen to be labelled as mitochondrion based on % identity (> 98.6%) and the largest alignments:" $seq_mit >> mito_sequences.out
		fi
		if [ "$(grep -c 'Not a valid version' blast_mito_log_out.txt)" -gt 0 ]; then
			echo -e "\nDOUBLE CHECK WHICH VERSIONS OF BLAST AND DATABASES ARE YOU USING, IT SEEMS THERE'S SOME MISMATCH AND SOMETHING FAILED, CHECK blast_mito_log_out.txt\n"
  		fi
	# 4. A database of ribosomal RNA genes
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" parallel --verbose --joblog blast_rrna_log_out_2.txt -j $blocks_size blastn -query {} -db $databases/rrna -task megablast -template_length 18 -template_type coding -window_size 120 -word_size 12 -xdrop_gap 20 -no_greedy -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 1E-9 -gapextend 2 -gapopen 4 -penalty -4 -perc_identity 95 -reward 3 -soft_masking true -outfmt 7 -num_threads $((cores / blocks_size)) -out {}.rrna_genbank.out ::: ${arr[@]} &> blast_rrna_log_out.txt
		awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' blast_rrna_log_out_2.txt > tmp && mv tmp blast_rrna_log_out_2.txt
		cat *.fa.rrna_genbank.out | egrep -v "^#" | awk '$4>=100' > rrna_genbank.out; rm *.fa.rrna_genbank.out
	# 5. The chromosomes of unrelated organisms. Foreign organisms are those that belong to a different taxonomic group compared to the organism whose sequences are being screened.
	# This is a requirement by the NCBI's Foreign Contamination Screen, but we skipped it here, because ILRA can alternatively use kraken2 to get rid of contamination of other organisms.
	# blastn -query _input_fasta_sequences_ -db _distant_organism_dbs_ -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -min_raw_gapped_score 100 -penalty -5 -perc_identity 98.0 -soft_masking true
	# Process the result and get the final files:
		awk ' NF>2 {print $1"\t"$7"\t"$8} ' *genbank.out | grep -v "#" | bedtools sort -i - > sequences_from_blast_to_mask.bed
		bedtools maskfasta -fi ../../$name.ILRA.fasta -bed sequences_from_blast_to_mask.bed -fo $name.ILRA_masked.fasta -fullHeader
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $name.ILRA_masked.fasta > $name.ILRA.fasta; rm $name.ILRA_masked.fasta
		rm ../../$name.ILRA.fasta
	# Remove if there are other contigs matching the mitochondrion:
		if [ ! -z "$seq_mit" ]; then
			if [ "$(cat mito_sequences | cut -f1 | sort | uniq | grep -v $seq_mit | wc -l)" -gt 0 ]; then
				echo -e "\n### Excluded contigs based on blast against mitochondrion sequences:" >> ../../Excluded.contigs.fofn
				cat mito_sequences | cut -f1 | sort | uniq | grep -v $seq_mit >> ../../Excluded.contigs.fofn
				cat mito_sequences | cut -f1 | sort | uniq | grep -v $seq_mit > mito_sequences_to_remove
				sed -i "s/${seq_mit}/${seq_mit} [location=mitochondrion]/" $name.ILRA.fasta
				samtools faidx $name.ILRA.fasta
				contigs_to_retain=$(awk '{print $1}' $name.ILRA.fasta.fai | grep -v -f mito_sequences_to_remove)
				samtools faidx $name.ILRA.fasta -o ../../$name.ILRA.fasta $contigs_to_retain; rm $name.ILRA.fasta*
			else
				mv $name.ILRA.fasta ../../$name.ILRA.fasta
			fi
		else
			mv $name.ILRA.fasta ../../$name.ILRA.fasta
		fi
		echo -e "\nAfter this step:"; assembly-stats ../../$name.ILRA.fasta | head -n 2
		echo -e "\nSTEP 6 BLAST: DONE"; echo -e "Current date/time: $(date)"
		time2=`date +%s`; echo -e "STEP 6 blast time (secs): $((time2-time1))"; echo -e "STEP 6 blast time (hours): $(echo "scale=2; $((time2-time1))/3600" | bc -l)"
		if [[ $kraken2_fast == "yes" ]]; then
			echo -e "\nRemember to double check or manually execute rm -rf /dev/shm/$kraken2_databases if you have selected kraken2_fast and won't be using ILRA again soon. The kraken2 database $kraken2_databases may be left allocated the RAM memory\n"
		fi	
	fi
	if [ $mode == "light" ]; then
		echo -e "\nSTEP 6 for decontamination is skipped. Light mode activated"; echo -e "Current date/time: $(date)\n"
	# Renaming and making single-line fasta with length in the headers
		if [ -z "$reference" ]; then
			cat $dir/5.Circlator/05.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_".$n } else { print }' | ILRA.fasta2singleLine.pl - | awk '/^>/ { if (name) {printf("%s_%d\n%s", name, len, seq)} name=$0; seq=""; len = 0; next} NF > 0 {seq = seq $0 "\n"; len += length()} END { if (name) {printf("%s_%d\n%s", name, len, seq)} }' > $dir/$name.ILRA.fasta
		else
			cat $dir/5.Circlator/05.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_with_ref_".$n } else { print }' | ILRA.fasta2singleLine.pl - | awk '/^>/ { if (name) {printf("%s_%d\n%s", name, len, seq)} name=$0; seq=""; len = 0; next} NF > 0 {seq = seq $0 "\n"; len += length()} END { if (name) {printf("%s_%d\n%s", name, len, seq)} }' > $dir/$name.ILRA.fasta
		fi
	fi
debug=""
fi

#### 7. Evaluate the assemblies, get telomere sequences counts, GC stats, sequencing depth, converting files...
mkdir -p $dir/data; cd $dir/data; ln -sf $assembly .; ln -sf $illuminaReads\_1.fastq.gz .; ln -sf $illuminaReads\_2.fastq.gz .

if [[ $debug == "all" || $debug == "step7" || $quality_step == "yes" ]]; then
	test_installation_dependencies.sh
	if [[ ! -d "$dir/6.Decontamination" ]] && [[ ! -d "$dir/5.Circlator" ]]; then
		ln -sf $assembly $dir/$name.ILRA.fasta
	fi
	time1=`date +%s`
	echo -e "\n\n\nSTEP 7 starting: Renaming, gathering stats and evaluation, cleaning up files..."; echo -e "Current date/time: $(date)"
	echo -e "This is the final ILRA step and it may take long, but the assembly has already been corrected and won't change more. You can use already the final corrected assembly "$dir"/"$name.ILRA.fasta
	mkdir -p $dir/7.Stats; cd $dir/7.Stats; rm -rf *

	# Evaluating the assemblies:
	echo -e "\nA preview of the final corrections by ILRA is: (full details in the file 07.assembly_stats_original_correction_ILRA.txt)"
	echo "Original:"; assembly-stats $assembly | head -n 3; assembly-stats $assembly > 07.assembly_stats_original_correction_ILRA.txt
	echo "Corrected:"; assembly-stats ../$name.ILRA.fasta | head -n 3; assembly-stats ../$name.ILRA.fasta >> 07.assembly_stats_original_correction_ILRA.txt

	# Getting telomere counts:
	echo -e "\nCheck out the output of the telomeres analyses in the file 07.TelomersSummary.txt"
	# Extracting potential telomeres (1Kb)
	cat ../$name.ILRA.fasta | perl -nle 'if (/>/){ print } else { $l="Left:"; $left=substr($_,0,1000); print "$l\t$left"}' | tr "\n" "\t" | sed 's/>/\n&/g' | sed '1d' > 07.Telomeres_seq_1kb.txt
	cat ../$name.ILRA.fasta | perl -nle 'if (/>/){ print } else { $r="Right:"; $right=substr($_,(length($_)-1000),1000); print "$r\t$right"}' | tr "\n" "\t" | sed 's/>/\n&/g' >> 07.Telomeres_seq_1kb.txt
	# Get contigs with the specified sequence repetitions in those regions
	echo -e "### Sequences with the left telomere sequence "$telomere_seq_1": (within a 1Kb end region, extracted in 07.Telomeres_seq_1kb.txt)" > 07.TelomerContigs.txt
	cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i $telomere_seq_1 | awk '{print $1}' | sort >> 07.TelomerContigs.txt
	echo -e "\n### Sequences with the right telomere sequence "$telomere_seq_2": (within a 1Kb end region, extracted in 07.Telomeres_seq_1kb.txt)" >> 07.TelomerContigs.txt
	cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i $telomere_seq_2 | awk '{print $1}' | sort  >> 07.TelomerContigs.txt
	# Summary of telomere repeats statistics
	echo -e "# Amount of telomere repeats left: "$(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i -c $telomere_seq_1) > 07.TelomeresSummary.txt
	echo -e "# Amount of telomere repeats right: "$(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i -c $telomere_seq_2) >> 07.TelomeresSummary.txt
	echo -e "# Amount of telomere repeats: "$(expr $(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i -c $telomere_seq_1) + $(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i -c $telomere_seq_2)) >> 07.TelomeresSummary.txt
	contigs_both=$(echo $(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i $telomere_seq_1 | awk '{print $1}') $(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i $telomere_seq_2 | awk '{print $1}') | tr ' ' '\n' | sort | uniq -d)
	echo -e "# Amount of telomere repeats at contigs with repeats at both ends: "$(expr $(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep "$contigs_both" | grep -i -c $telomere_seq_1) + $(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep "$contigs_both" | grep -i -c $telomere_seq_2)) >> 07.TelomeresSummary.txt
	echo -e "\n# Amount of contigs with telomere repeats left: "$(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i $telomere_seq_1 | awk '{print $1}' | wc -l) >> 07.TelomeresSummary.txt
	echo -e "# Amount of contigs with telomere repeats right: "$(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i $telomere_seq_2 | awk '{print $1}' | wc -l) >> 07.TelomeresSummary.txt
	if [ -z "$contigs_both" ]; then
		echo -e "# Amount of contigs/potential chromosomes with telomere repeats at both ends: 0" >> 07.TelomeresSummary.txt
		echo -e "\n# Contigs/potential chromosomes with telomere repeats at both ends: (Length, bp)" >> 07.TelomeresSummary.txt
	else
		echo -e "# Amount of contigs/potential chromosomes with telomere repeats at both ends: "$(echo $contigs_both | tr " " "\n" | wc -l) >> 07.TelomeresSummary.txt
		echo -e "\n# Contigs/potential chromosomes with telomere repeats at both ends: (Length, bp)" >> 07.TelomeresSummary.txt
		awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../$name.ILRA.fasta | tr "\n" "\t" | sed 's/>/\n&/g' | grep "$contigs_both" | sort >> 07.TelomeresSummary.txt
	fi

	# Getting GC stats: (https://github.com/oXis/gtool)
	echo -e "\nA preview of the overall GC content is: (full details in the file 07.Contigs_GC_size.txt)"
	gtool.py -sg ../$name.ILRA.fasta | tee 07.Contigs_GC_size.txt; echo -e "\nIndividual contigs:" >> 07.Contigs_GC_size.txt
	for i in $(grep "^>" ../$name.ILRA.fasta | awk 'sub(/^>/, "")'); do
		cat ../$name.ILRA.fasta | gtool.py -sg -c $i - >> 07.Contigs_GC_size.txt
	done
	sed -i 's/stdin//g' 07.Contigs_GC_size.txt

	# Getting sequencing depth and stats: (modified from https://github.com/wudustan/fastq-info)
	echo -e "\nCheck out the sequencing depth and other stats in the file 07.fastq_info_depth_IlluminaReads_assembly.txt"
	if [ -f $illuminaReads\_1.fastq.gz ]; then
		fastq-info.sh $illuminaReads\_1.fastq.gz $illuminaReads\_2.fastq.gz ../$name.ILRA.fasta > 07.fastq_info_depth_IlluminaReads_assembly.txt
		fastqc $illuminaReads\_1.fastq.gz $illuminaReads\_2.fastq.gz -o $PWD --noextract -q -t $cores
	fi

	# Getting whether the organism is eukaryote for the following steps:
	if [[ $mode == "light" || $mode == "blast" ]]; then
		mkdir -p $databases/taxdump && cd $databases/taxdump && wget -q https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip && unzip -o -qq taxdmp.zip && rm taxdmp.zip
		re='^[0-9]+$' # Deal with the user providing taxon id or taxon names
		if [[ $taxonid =~ $re ]]; then
			taxon_name=$(taxonkit list --ids $taxonid -n -r --data-dir $databases/taxdump | grep $taxonid)
		else
			taxonid=$(echo $taxonid | taxonkit name2taxid --data-dir $databases/taxdump | head -1 | cut -f2)
			taxon_name=$(taxonkit list --ids $taxonid -n -r --data-dir $databases/taxdump | grep $taxonid)
		fi
	fi

	if [ "$(echo $taxonid | taxonkit lineage --data-dir $databases/taxdump | grep "Eukaryota" | wc -l)" -gt 0 ]; then
		top_level="E"
	elif [ "$(echo $taxonid | taxonkit lineage --data-dir $databases/taxdump | grep "Prokaryota" | wc -l)" -gt 0 ]; then
		top_level="P"
	fi

	# Comparing with reference genes (QUAST):
	echo -e "\nRunning QUAST..."
	echo -e "Please be aware that for providing reference genes, a GFF file with gene or operon as feature type field, or a bed file (sequence name, start position, end position, gene ID) are accepted"
	echo -e "Please be aware that ILRA is automatically checking if the provided NCBI taxon ID is eukaryotic or not, to use the '--eukaryote' argument in quast.py. If your species is prokaryotic QUAST would also work. In other cases, you may need to manually run quast.py with the argument '--fungus' for fungi or the arguments '--large' and '--memory-efficient' for large genomes. If the taxon ID is not known or present in the NCBI taxonomy databases, this step will be skipped..."
	echo -e "Current NCBI taxon ID: $taxonid"
	cd $dir/7.Stats
	if [ "$top_level"=="E" ]; then
		if [ -z "$reference" ]; then
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta --threads $cores --labels $name --eukaryote &> quast.py_log_out.txt
		else
			if [ -f $illuminaReads\_1.fastq.gz ]; then
				if [ -z "$gff_file" ]; then					
					if [ ! -z "$correctedReads" ]; then
						if [[ $correctedReads == *.fasta.gz || $correctedReads == *.fa.gz ]]; then
							pigz -p $cores -dc $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $(basename $correctedReads).fq.gz --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $(basename $correctedReads).fq.gz --eukaryote &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq.gz
						elif [[ $correctedReads == *.fasta || $correctedReads == *.fa ]]; then
							seqtk seq -F '#' $correctedReads | pigz -c -p $cores --fast > $(basename $correctedReads).fq
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $(basename $correctedReads).fq --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $(basename $correctedReads).fq --eukaryote &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq
						elif [[ $correctedReads == *.fastq.gz || $correctedReads == *.fq.gz || $correctedReads == *.fastq || $correctedReads == *.fq ]]; then
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $correctedReads --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $correctedReads --eukaryote &> quast.py_log_out.txt
							fi
						fi
					else
						\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --eukaryote &> quast.py_log_out.txt
					fi
				else
					if [ ! -z "$correctedReads" ]; then
						if [[ $correctedReads == *.fasta.gz || $correctedReads == *.fa.gz ]]; then
							pigz -p $cores -dc $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $(basename $correctedReads).fq.gz --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $(basename $correctedReads).fq.gz --eukaryote &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq.gz
						elif [[ $correctedReads == *.fasta || $correctedReads == *.fa ]]; then
							seqtk seq -F '#' $correctedReads | pigz -c -p $cores --fast > $(basename $correctedReads).fq
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $(basename $correctedReads).fq --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $(basename $correctedReads).fq --eukaryote &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq
						elif [[ $correctedReads == *.fastq.gz || $correctedReads == *.fq.gz || $correctedReads == *.fastq || $correctedReads == *.fq ]]; then
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $correctedReads --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $correctedReads --eukaryote &> quast.py_log_out.txt
							fi
						fi
					else
						\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --eukaryote &> quast.py_log_out.txt
					fi
				fi
			else
				if [ -z "$gff_file" ]; then					
					if [ ! -z "$correctedReads" ]; then
						if [[ $correctedReads == *.fasta.gz || $correctedReads == *.fa.gz ]]; then
							pigz -p $cores -dc $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pacbio $(basename $correctedReads).fq.gz --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --nanopore $(basename $correctedReads).fq.gz --eukaryote &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq.gz
						elif [[ $correctedReads == *.fasta || $correctedReads == *.fa ]]; then
							seqtk seq -F '#' $correctedReads | pigz -c -p $cores --fast > $(basename $correctedReads).fq
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pacbio $(basename $correctedReads).fq --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --nanopore $(basename $correctedReads).fq --eukaryote &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq
						elif [[ $correctedReads == *.fastq.gz || $correctedReads == *.fq.gz || $correctedReads == *.fastq || $correctedReads == *.fq ]]; then
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pacbio $correctedReads --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --nanopore $correctedReads --eukaryote &> quast.py_log_out.txt
							fi
						fi
					else
						\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --eukaryote &> quast.py_log_out.txt
					fi
				else
					if [ ! -z "$correctedReads" ]; then
						if [[ $correctedReads == *.fasta.gz || $correctedReads == *.fa.gz ]]; then
							pigz -p $cores -dc $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pacbio $(basename $correctedReads).fq.gz --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --nanopore $(basename $correctedReads).fq.gz --eukaryote &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq.gz
						elif [[ $correctedReads == *.fasta || $correctedReads == *.fa ]]; then
							seqtk seq -F '#' $correctedReads | pigz -c -p $cores --fast > $(basename $correctedReads).fq
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pacbio $(basename $correctedReads).fq --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --nanopore $(basename $correctedReads).fq --eukaryote &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq
						elif [[ $correctedReads == *.fastq.gz || $correctedReads == *.fq.gz || $correctedReads == *.fastq || $correctedReads == *.fq ]]; then
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pacbio $correctedReads --eukaryote &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --nanopore $correctedReads --eukaryote &> quast.py_log_out.txt
							fi
						fi
					else
						\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --eukaryote &> quast.py_log_out.txt
					fi
				fi
			fi
		fi			
	echo -e "Check out the file quast.py_log_out.txt and the QUAST report within the folder 7.Stats/quast_results"
	elif [ "$top_level"=="P" ]; then
		if [ -z "$reference" ]; then
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta --threads $cores --labels $name &> quast.py_log_out.txt
		else
			if [ -f $illuminaReads\_1.fastq.gz ]; then
				if [ -z "$gff_file" ]; then					
					if [ ! -z "$correctedReads" ]; then
						if [[ $correctedReads == *.fasta.gz || $correctedReads == *.fa.gz ]]; then
							pigz -p $cores -dc $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $(basename $correctedReads).fq.gz &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $(basename $correctedReads).fq.gz &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq.gz
						elif [[ $correctedReads == *.fasta || $correctedReads == *.fa ]]; then
							seqtk seq -F '#' $correctedReads | pigz -c -p $cores --fast > $(basename $correctedReads).fq
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $(basename $correctedReads).fq &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $(basename $correctedReads).fq &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq
						elif [[ $correctedReads == *.fastq.gz || $correctedReads == *.fq.gz || $correctedReads == *.fastq || $correctedReads == *.fq ]]; then
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $correctedReads &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $correctedReads &> quast.py_log_out.txt
							fi
						fi
					else
						\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz &> quast.py_log_out.txt
					fi
				else
					if [ ! -z "$correctedReads" ]; then
						if [[ $correctedReads == *.fasta.gz || $correctedReads == *.fa.gz ]]; then
							pigz -p $cores -dc $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $(basename $correctedReads).fq.gz &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $(basename $correctedReads).fq.gz &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq.gz
						elif [[ $correctedReads == *.fasta || $correctedReads == *.fa ]]; then
							seqtk seq -F '#' $correctedReads | pigz -c -p $cores --fast > $(basename $correctedReads).fq
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $(basename $correctedReads).fq &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $(basename $correctedReads).fq &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq
						elif [[ $correctedReads == *.fastq.gz || $correctedReads == *.fq.gz || $correctedReads == *.fastq || $correctedReads == *.fq ]]; then
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --pacbio $correctedReads &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz --nanopore $correctedReads &> quast.py_log_out.txt
							fi
						fi
					else
						\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pe1 $illuminaReads\_1.fastq.gz --pe2 $illuminaReads\_2.fastq.gz &> quast.py_log_out.txt
					fi
				fi
			else
				if [ -z "$gff_file" ]; then					
					if [ ! -z "$correctedReads" ]; then
						if [[ $correctedReads == *.fasta.gz || $correctedReads == *.fa.gz ]]; then
							pigz -p $cores -dc $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pacbio $(basename $correctedReads).fq.gz &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --nanopore $(basename $correctedReads).fq.gz &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq.gz
						elif [[ $correctedReads == *.fasta || $correctedReads == *.fa ]]; then
							seqtk seq -F '#' $correctedReads | pigz -c -p $cores --fast > $(basename $correctedReads).fq
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pacbio $(basename $correctedReads).fq &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --nanopore $(basename $correctedReads).fq &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq
						elif [[ $correctedReads == *.fastq.gz || $correctedReads == *.fq.gz || $correctedReads == *.fastq || $correctedReads == *.fq ]]; then
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --pacbio $correctedReads &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --nanopore $correctedReads &> quast.py_log_out.txt
							fi
						fi
					else
						\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name &> quast.py_log_out.txt
					fi
				else
					if [ ! -z "$correctedReads" ]; then
						if [[ $correctedReads == *.fasta.gz || $correctedReads == *.fa.gz ]]; then
							pigz -p $cores -dc $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pacbio $(basename $correctedReads).fq.gz &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --nanopore $(basename $correctedReads).fq.gz &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq.gz
						elif [[ $correctedReads == *.fasta || $correctedReads == *.fa ]]; then
							seqtk seq -F '#' $correctedReads | pigz -c -p $cores --fast > $(basename $correctedReads).fq
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pacbio $(basename $correctedReads).fq &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --nanopore $(basename $correctedReads).fq &> quast.py_log_out.txt
							fi
							rm $(basename $correctedReads).fq
						elif [[ $correctedReads == *.fastq.gz || $correctedReads == *.fq.gz || $correctedReads == *.fastq || $correctedReads == *.fq ]]; then
							if [ $long_reads_technology == "pacbio" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --pacbio $correctedReads &> quast.py_log_out.txt
							elif [ $long_reads_technology == "ont2d" ]; then
								\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --nanopore $correctedReads &> quast.py_log_out.txt
							fi
						fi
					else
						\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name &> quast.py_log_out.txt
					fi
				fi
			fi
		fi
	echo -e "Check out the file quast.py_log_out.txt and the QUAST report within the folder 7.Stats/quast_results"
	fi

	# Asessing genome completeness (BUSCO):
	echo -e "\nRunning BUSCO..."
	echo -e "Please be aware that ILRA is automatically checking if the provided NCBI taxon ID is eukaryotic or not to use the automatic lineage selection in BUSCO. This may not work and you would need to run manually BUSCO with the final corrected assembly, providing as the argument '--lineage_dataset' the BUSCO dataset closest to your species in 'https://busco-data.ezlab.org/v4/data/lineages/', which would be automatically downloaded. If the taxon ID is not known or present in the NCBI taxonomy databases, the automatic mode will be executed"
	echo -e "Current NCBI taxon ID: $taxon_name"
	echo -e "Check out the file busco_log_out.txt and the BUSCO reports within the folder 7.Stats/busco_results"
	if [ ! -z "$quality_step_busco_db" ]; then
		echo -e "Current BUSCO db: $quality_step_busco_db"
		echo -e "Running BUSCO with --lineage_dataset $quality_step_busco_db"
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" busco -i ../$name.ILRA.fasta -o $name -m genome -f -c $cores --lineage_dataset $quality_step_busco_db --tar --out_path $dir/7.Stats/busco_results --download_path $dir/7.Stats/busco_results/downloads &> busco_log_out.txt
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" busco -i $assembly -o $(echo $(basename $assembly) | sed 's,.fa.*,,g')"_preILRA" -m genome -f -c $cores --lineage_dataset $quality_step_busco_db --tar --out_path $dir/7.Stats/busco_results --download_path $dir/7.Stats/busco_results/downloads &> busco_log_out.txt
	else
		if [ "$top_level"=="E" ]; then
			echo -e "Running BUSCO for eukaryotes in the mode '--auto-lineage-euk'"
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" busco -i ../$name.ILRA.fasta -o $name -m genome -f -c $cores --auto-lineage-euk --tar --out_path $dir/7.Stats/busco_results --download_path $dir/7.Stats/busco_results/downloads &> busco_log_out.txt
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" busco -i $assembly -o $(echo $(basename $assembly) | sed 's,.fa.*,,g')"_preILRA" -m genome -f -c $cores --auto-lineage-euk --tar --out_path $dir/7.Stats/busco_results --download_path $dir/7.Stats/busco_results/downloads &> busco_log_out.txt
		elif [ "$top_level"=="P" ]; then
			echo -e "Running BUSCO for prokaryotes in the mode '--auto-lineage-prok'"
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" busco -i ../$name.ILRA.fasta -o $name -m genome -f -c $cores --auto-lineage-prok --tar --out_path $dir/7.Stats/busco_results --download_path $dir/7.Stats/busco_results/downloads &> busco_log_out.txt
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" busco -i $assembly -o $(echo $(basename $assembly) | sed 's,.fa.*,,g')"_preILRA" -m genome -f -c $cores --auto-lineage-prok --tar --out_path $dir/7.Stats/busco_results --download_path $dir/7.Stats/busco_results/downloads &> busco_log_out.txt
		else
			echo -e "Running BUSCO in the automatic mode, '--auto-lineage'"
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" busco -i ../$name.ILRA.fasta -o $name -m genome -f -c $cores --auto-lineage --tar --out_path $dir/7.Stats/busco_results --download_path $dir/7.Stats/busco_results/downloads &> busco_log_out.txt
			\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" busco -i $assembly -o $(echo $(basename $assembly) | sed 's,.fa.*,,g')"_preILRA" -m genome -f -c $cores --auto-lineage --tar --out_path $dir/7.Stats/busco_results --download_path $dir/7.Stats/busco_results/downloads &> busco_log_out.txt
		fi
	fi
	rm -rf $dir/7.Stats/busco_results/downloads
	busco_cogeqc_plots.R $dir/7.Stats/busco_results
	cat $(find . -name busco.log) | grep "|" | grep -v "#"

	# Visualisation of synteny and structural rearrangements against reference (plotsr):
	if [ ! -z "$reference" ]; then
		echo -e "\nRunning plotsr..."
		mkdir -p $dir/7.Stats/plotsr_results; cd $dir/7.Stats/plotsr_results/
		# minimap2 -ax asm5 -t $cores --eqx ../$name.ILRA.fasta $reference | samtools sort -l 9 -m 2G -@ $cores --write-index -o "$dir/7.Stats/plotsr_results/plotsr.bam##idx##$dir/7.Stats/plotsr_results/plotsr.bam.bai"
		# (minimap2 can be used and it's faster, but meryl+winnowmap is currently recommended to map against repetitive sequences)
		ln -sf ../../$name.ILRA.fasta assembly_ILRA.fasta; ln -sf $reference ref.fa # only the matching contigs can be shown with this visualization, so prepare copies renaming the contigs
		paste <(for f in $(grep ">" ref.fa | sort); do grep ">" assembly_ILRA.fasta | sort | grep ${f:1} | sed -r 's/^.{1}//'; done) <(for f in $(grep ">" ref.fa | sort); do if [ $(echo $(grep ">" assembly_ILRA.fasta | sort | grep ${f:1} | sed -r 's/^.{1}//') | sed '/^$/d' | wc -l) -gt 0 ]; then echo ${f:1}; fi; done) | awk ' NF==2 {print $0} ' > contigs_pairs.txt
		seqkit replace -p "(.+)" -r '{kv}' -k contigs_pairs.txt -j $cores -m to_remove --quiet assembly_ILRA.fasta | awk '!/to_remove/' RS=">" ORS=">" | sed 's,>$,,g' > assembly_ILRA_reduced_comp_ref.fasta
		ref_to_rm=$(comm -23 <(grep ">" ref.fa | sed 's,^>,,g') <(cat contigs_pairs.txt | awk '{print $2}'))
		for f in $(comm -23 <(cat ref.fa | grep ">") <(echo $ref_to_rm | tr ' ' '\n' | sed 's/^/>/')); do ILRA.fasta2singleLine.pl ref.fa | grep -A1 $f >> ref_reduced_comp_assembly.fasta; done
		samtools faidx ref_reduced_comp_assembly.fasta; \time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" picard CreateSequenceDictionary -R ref_reduced_comp_assembly.fasta -O ref_reduced_comp_assembly.fasta.dict &> indexing_ref_log_out.txt # Prepare reference to fix cigar
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" meryl count k=19 threads=$cores output merylDB assembly_ILRA_reduced_comp_ref.fasta &> meryl_count_log_out.txt
		meryl print greater-than distinct=0.9998 merylDB > repetitive_k19.txt 2> meryl_print_log_out.txt
		winnowmap -W repetitive_k19.txt -ax asm5 -t $cores ref_reduced_comp_assembly.fasta assembly_ILRA_reduced_comp_ref.fasta 2> winnowmap_log_out.txt | samtools view -@ $cores -h -F 260 -F 2048 -o plotsr.sam
		javarkit_samfixcigar=$(for path in $(echo ${PATH//:/ } | cut -d" " -f1,2,3,4,5); do find $path -name samfixcigar.jar; done | head -1) # Look for the jar file in the first 5 folders in $PATH (which if installed with conda will capture ILRA_env/bin)
		$(dirname $javarkit_samfixcigar)/java -jar $javarkit_samfixcigar -o plotsr_fix.sam -R ref_reduced_comp_assembly.fasta plotsr.sam # Updated java version and the precise pathway to the jar are required
		# Get structural rearrangements (syri) and plotsr
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" syri -c plotsr_fix.sam -r ref_reduced_comp_assembly.fasta -q assembly_ILRA_reduced_comp_ref.fasta -F B --nc $cores --seed 1 --dir . --prefix plotsr_ &> syri_log_out.txt # syri fails if there are M in the bam cigar string, hence the reformat above
		echo -e "#file\tname\ttags\n$PWD/ref_reduced_comp_assembly.fasta\treference\tlw:1.5\n$PWD/assembly_ILRA_reduced_comp_ref.fasta\tassembly\tlw:1.5" > genomes.txt
		plotsr --sr plotsr_syri.out --genomes genomes.txt -o plotsr_assembly_reference_plot.pdf
		rm -rf $(ls | egrep -v ".pdf$|.out$|.fasta$|.txt$") $dir/7.Stats/plotsr_results/merylDB # Cleaning
	fi

	# Visualization of syntenic relationships on multiple contigs (NGenomeSyn):
	if [ ! -z "$reference" ]; then
		mkdir -p $dir/7.Stats/NGenomeSyn; cd $dir/7.Stats/NGenomeSyn/
		a=../../$name.ILRA.fasta
		b=$reference # Get individual contigs and iterate
		cat $a | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) "_first.fa")} print $0 >> filename close(filename) }'
		cat $b | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) "_second.fa")} print $0 >> filename close(filename) }'
		for f in $(ls | egrep .fa-.*); do mv $f $(echo $f | sed 's,.fa-.*,.fa,g'); done
		for f in $(ls | egrep _first.fa$); do
			for i in $(ls | egrep _second.fa$); do
				perl $(dirname $(which GetTwoGenomeSyn.pl))/NGenomeSyn-1.41/bin/GetTwoGenomeSyn.pl -InGenomeA $f -InGenomeB $i -OutPrefix ${f%.*}vs${i%.*} -MappingBin minimap2 -MinLenA $(assembly-stats $f | grep sum | sed 's/^sum = //g;s/, .*//g') -MinLenB $(assembly-stats $i | grep sum | sed 's/^sum = //g;s/, .*//g') -NumThreads $cores &>> GetTwoGenomeSyn.pl_log_out.txt
			done
		done
		rm $(ls | egrep .fa$)
		tar -cf NGenomeSyn_full_out.tar . &>> GetTwoGenomeSyn.pl_log_out.txt
		rm $(ls | egrep -v '.png$|.tar$|_log_out.txt$')
	fi

	# Assembly metric visualisations (https://github.com/rjchallis/assembly-stats)
	mkdir -p $dir/7.Stats/asm2stats_results; cd $dir/7.Stats/asm2stats_results
	asm2stats.pl ../../$name.ILRA.fasta > output.json
	asm2stats.minmaxgc.pl ../../$name.ILRA.fasta > output.minmaxgc.json
	# https://github.com/blobtoolkit/blobtoolkit/wiki/Snail-plot#command-line
 	mkdir -p $dir/7.Stats/blobdir_results; cd $dir/7.Stats/blobdir_results
  	blobtools create --fasta $(ls -d ../../* | egrep .fasta$) --busco $(find ../busco_results/ -name "full_table.tsv" | grep -v auto_lineage) $PWD &> blobtools.log
   	blobtools view --plot --view snail $PWD

	# Visualization of telomeres and corrected reads coverage (tapestry)
	if [ ! -z "$correctedReads" ]; then
		mkdir -p $dir/7.Stats/tapestry_results; cd $dir/7.Stats/tapestry_results/
		if [[ $correctedReads == *.fasta.gz || $correctedReads == *.fa.gz ]]; then
			pigz -p $cores -dc $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
			weave -a ../../$name.ILRA.fasta -r $(basename $correctedReads).fq.gz -t TTAGGG CTTATT $telomere_seq_1 $telomere_seq_2 -o out -c $cores &> weave_log_out.txt
		elif [[ $correctedReads == *.fastq.gz || $correctedReads == *.fq.gz ]]; then
			weave -a ../../$name.ILRA.fasta -r $correctedReads -t TTAGGG CTTATT $telomere_seq_1 $telomere_seq_2 -o out -c $cores &> weave_log_out.txt
		elif [[ $correctedReads == *.fasta || $correctedReads == *.fa ]]; then
			cat $correctedReads | seqtk seq -F '#' - | pigz -c -p $cores --fast > $(basename $correctedReads).fq.gz
			weave -a ../../$name.ILRA.fasta -r $(basename $correctedReads).fq.gz -t TTAGGG CTTATT $telomere_seq_1 $telomere_seq_2 -o out -c $cores &> weave_log_out.txt
		fi
		cd out; rm -rf $(ls | egrep -v '.tsv$|.html$')
	fi

	time2=`date +%s`; echo -e "\nSTEP 7 time (secs): $((time2-time1))"; echo -e "STEP 7 time (hours): $(echo "scale=2; $((time2-time1))/3600" | bc -l)"
	echo -e "\n\nSTEP 7: DONE"; echo -e "Current date/time: $(date)"
else
	echo "The STEP 7 for assessing the quality of the corrected assembly, gathering sequences, analyzing telomeres... have been skipped"
	echo -e "Current date/time: $(date)"
fi

# Converting files to minimize space and tidying up
echo -e "\nConverting and compressing final files..."
if [ -f $dir/2.MegaBLAST/first.bam ]; then
	samtools view -@ $cores -T $dir/1.Filtering/01.assembly.fa -C -o $dir/2.MegaBLAST/first.bam.cram $dir/2.MegaBLAST/first.bam
fi
if [[ -d "$dir/4.iCORN2" ]]; then
	for f in $(find $dir/4.iCORN2 -type f -name "*.bam"); do
		a=${f%/*}
		samtools view -@ $cores -T $dir/4.iCORN2/ICORN2.03b.assembly.fa.${a##*_} -C -o "${f}".cram "${f}"
	done
fi
if [[ -d "$dir/4.pilon" ]]; then
	for i in $(seq 1 1 $number_iterations_icorn); do
		\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" samtools view -@ $cores -T $dir/3.ABACAS2/03b.assembly.fa -C -o $dir/4.pilon/ill_reads$i.bam.cram ill_reads$i.bam &> $dir/4.pilon/ill_reads$i.bam_log_out.txt
	done
fi
if [ -f $dir/5.Circlator/Mapped.corrected.04.sam ]; then
	samtools view -@ $cores -T $(find $dir -name "04.assembly.fa") -C -o $dir/5.Circlator/Mapped.corrected.04.sam.cram $dir/5.Circlator/Mapped.corrected.04.sam	
fi
if [ -f $dir/5.Circlator/Out.Circ/01.mapreads.bam ]; then
	samtools view -@ $cores -T $dir/5.Circlator/ForCirc.Ref_2.fasta -C -o $dir/5.Circlator/Out.Circ/01.mapreads.bam.cram $dir/5.Circlator/Out.Circ/01.mapreads.bam
fi
for i in $(find $dir -regex '.*\(.fq$\|.fastq$\|.fa$\|.fasta$\|.fna$\|.plot$\|.vcf$\|ICORN2.03b.assembly.fa.[1-10]\)$' | grep -v $name.ILRA.fasta | grep -v "01.assembly.fa" | grep -v "03.assembly.fa" | grep -v "03b.assembly.fa" | grep -v "04.assembly.fa" | grep -v "05.assembly.fa"); do
	pigz -f -p $cores --best $i
done
echo -e "\nMain alignment files have been converted to cram for long-term storage. If needed, for converting compressed .cram files back to .bam apply the command: samtools view -@ $cores -T filename.fasta -b -o output.bam input.cram (check out samtools view statements within ILRA.sh to get the fasta file used)"
# Cleaning up:
if [ "$(find $dir -regex '.*\(.bam$\|.sam$\)$' | wc -l)" -gt 0 ]; then
	rm $(find $dir -regex '.*\(.bam$\|.sam$\)$')
fi

mkdir -p $dir/all_logs; cd $dir
for f in $(find . -type f -name "*out*"); do
	cp $f $dir/all_logs/$(echo $f | sed -e 's,/,_,g' -e 's,^._,,g')
done
rm $(find $dir -regex '.*\(.nhr$\|.nin$\|.nsq$\|.fai$\|.bai$\)$')

echo -e "\n\n\n\n\n(List of output files)" > $dir/all_output_files_log_out.txt
ls -Rslh $dir &> $dir/all_output_files_log_out.txt

echo -e "\nILRA IS DONE"; echo -e "Current date/time: $(date)\n"
echo -e "Original assembly file: "$assembly
echo -e "Final corrected assembly file: "$dir"/"$name.ILRA.fasta
echo -e "Excluded contigs file: "$dir"/Excluded.contigs.fofn"
echo -e "To assess the running please look for the expected output and log files within the step folders, particularly '*log_out*' files"
end=`date +%s`; runtime=$((end-start))
echo -e "\nCores: $cores"
echo -e "Final runtime (hours): $((runtime / 3600))"
echo -e "($runtime secs, check out time for each step after the statements 'STEP X: DONE' or in the logs the output of GNU's time)"
