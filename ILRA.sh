#!/bin/bash
# Authors: José Luis Ruiz & Thomas Dan Otto
# License: GNU General Public License v3.0
# https://github.com/ThomasDOtto/ILRA

########## ILRA
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
		-h*) echo "usage: ILRA.sh [options]
		-h | -help # Type this to get help
		-a | -assembly # Name of the long reads assembly to correct (FASTA format)
		-f | -filter_contig_size # Size threshold to filter the contigs (bp)
		-I | -Illumina_reads # Root name of the paired-end Illumina reads (FASTQ format, must be gzipped)
		-R | -Range_insert_size # Insert size range of the Illumina reads to use in mapping by SMALT (bp)
		-i | -iterations_iCORN2 # Number of iterations to perform in iCORN2
		-r | -reference # Reference file (full pathway, FASTA format)
		-g | -gff_file # Reference annotation file (full pathway, GFF format)
		-c | -corrected_reads # Corrected long reads (FASTQ format, can be gzipped)
		-s | -seq_circularize_1 # Regex pattern to find in the contig names and circularize
		-S | -Seq_circularize_2 # Regex pattern to find in the contig names and circularize
		-L | -Long_reads_technology # Technology of long reads sequencing (pb/ont)
		-T | -Taxon_ID # NCBI taxon ID to extract in the decontamination step
		-e | -ending_telomere_seq_1 # Telomere-associated sequence to search in the first strand
		-E | -Ending_telomere_seq_2 # Telomere-associated sequence to search in the second strand
		-o | -output # Output folder (full pathway)
		-n | -name # Base name of the output file
		-t | -threads # Number of cores to use in multithreaded steps
		-m | -mode # Add 'taxon' to execute decontamination based on taxonomic classification by Centrifuge, add 'blast' to execute decontamination based on BLAST against databases as requested by the DDBJ/ENA/Genbank submission, add 'both' to execute both approaches, and add 'light' to execute ILRA in light mode and skip these steps (default)" && exit 1;;
		-a*) assembly=${arguments[index]} ;;
		-o*) dir=${arguments[index]} ;;
		-c*) correctedReads=${arguments[index]} ;;
		-n*) name=${arguments[index]} ;;
		-r*) reference=${arguments[index]} ;;
		-I*) illuminaReads=${arguments[index]} ;;
		-s*) seqs_circl_1=${arguments[index]} ;;
		-S*) seqs_circl_2=${arguments[index]} ;;
		-i*) number_iterations_icorn=${arguments[index]} ;;
		-t*) cores=${arguments[index]} ;;
		-f*) contigs_threshold_size=${arguments[index]} ;;
		-R*) InsertsizeRange=${arguments[index]} ;;
		-T*) taxonid=${arguments[index]} ;;
		-g*) gff_file=${arguments[index]} ;;
		-e*) telomere_seq_1=${arguments[index]} ;;
		-E*) telomere_seq_2=${arguments[index]} ;;
		-L*) seq_technology=${arguments[index]} ;;
		-m*) mode=${arguments[index]} ;;
	esac
done
export name; export telomere_seq_1; export telomere_seq_2

echo -e "Check help with parameter '-h' or the read '/path/to/ILRA/test_data/README.txt' for full example usage"
echo -e "In case you want to test ILRA with a smaller subset of example reads, check the subfolder /path/to/ILRA/test_data"


##### Checking Arguments / Variables:
echo -e "\nI'm now quickly checking and showing the arguments that are going to be used in the ILRA run..."
if [[ $assembly == /* ]]; then
	echo -e "Assembly provided correctly"
elif [ -z "$assembly" ]; then
  echo -e "Assembly not provided"
	exit 1
else
	echo -e "Assuming the assembly is in the current pathway..."
	assembly=$PWD/$assembly
fi

if [[ $dir == /* ]]; then
	echo -e "Setting up the directory to work"
else
	dir=$PWD/$dir
fi

if [[ $correctedReads == /* ]]; then
	echo -e "Corrected long reads provided correctly"
elif [ -z "$correctedReads" ]; then
  echo -e "Corrected reads not provided"
	exit 1
else
	echo -e "Assuming the corrected reads are in the current pathway..."
	correctedReads=$PWD/$correctedReads
fi

if [[ $illuminaReads == /* ]]; then
	echo -e "Illumina short reads provided correctly"
elif [ -z "$illuminaReads" ]; then
	echo -e "Illumina reads not detected. Some steps of the pipeline will not be executed"
else
	echo -e "Trying to look for the Illumina reads in the current pathway"
	illuminaReads=$PWD/$illuminaReads
fi

if [ -f $illuminaReads\_1.fastq.gz ]; then
  echo "Good, ILRA is detecting the naming required for the Illumina reads: _1.fastq.gz and _2.fastq.gz. Dealing with them now..."
	mkdir -p $dir/1.Filtering; cd $dir/1.Filtering
	pigz -d -f -k -c -p $cores $illuminaReads\_1.fastq.gz > $dir/1.Filtering/"${illuminaReads##*/}"\_1.fastq
	pigz -d -f -k -c -p $cores $illuminaReads\_2.fastq.gz > $dir/1.Filtering/"${illuminaReads##*/}"\_2.fastq
	illuminaReads=$dir/1.Filtering/"${illuminaReads##*/}"
else
  echo -e "Please be aware of the naming required for the Illumina reads: _1.fastq.gz and _2.fastq.gz\n"
  echo "ILRA is not detecting the Illumina reads files. If you want to use Illumina reads for polishing, please check naming, paths and that the files do exist. Otherwise ignore this if you don't want to use Illumina reads. ILRA will skip some steps accordingly"
	doAbacas2=0
fi

if [ -z "$reference" ]; then
	echo "No reference genome provided. ABACAS2 will be skipped ..."
	doAbacas2=0
else
	doAbacas2=1
fi

if [ -z "mode" ]; then
	mode="light"
	echo "Ligth mode activated, steps for decontamination and preparation for online databases steps will be skipped ..."
else
	echo "ILRA execution mode (-m) is: " $mode
fi

echo -e "Final arguments used:"
if [ -z "$number_iterations_icorn" ]; then
	number_iterations_icorn=3
	echo "Number of iCORN2 iterations: "$number_iterations_icorn
fi

if [ -z "$cores" ] ; then
	cores=40
	echo "Number of threads: "$cores
fi

if [ -z "$seqs_circl_1" ] ; then
	seqs_circl_1="MT|M76611"
	echo "Seqs to circularize: "$seqs_circl_1
fi

if [ -z "$seqs_circl_2" ] ; then
	seqs_circl_2="API"
	echo "Seqs to circularize: "$seqs_circl_2
fi

if [ -z "$contigs_threshold_size" ] ; then
	contigs_threshold_size=5000
	echo "Length threshold to discard contigs (bp): "$contigs_threshold_size
fi

if [ -z "$InsertsizeRange" ] ; then
	InsertsizeRange=800
	echo "Insert size range Illumina short reads (bp): "$InsertsizeRange
fi

if [ -z "$taxonid" ] ; then
	taxonid=5833
	echo "NCBI taxon id to keep in decontamination step: "$taxonid
fi

if [ $seq_technology == "pb" ] ; then
	long_reads_technology="pacbio"
	echo "Long reads sequencing technology: PacBio"
	echo "If Oxford Nanopore were used, please provide the argument 'ont'"
elif [ $seq_technology == "ont" ] ; then
	long_reads_technology="ont2d"
	echo "Long reads sequencing technology: ONT"
	echo "If PacBio were used, please provide the argument 'pb'"
else
	long_reads_technology="pacbio"
	echo "Using by default long reads sequencing technology: PacBio"
	echo "If this needs to be change please provide 'pb' or 'ont' as the last argument"
fi
if [ -z "$telomere_seq_1" ] ; then
	telomere_seq_1="CCCTAAACCCTAAACCCTAAA"
fi
if [ -z "$telomere_seq_2" ] ; then
	telomere_seq_2="TTTAGGGTTTAGGGTTTAGGG"
fi
echo -e "assembly="$assembly
echo -e "dir="$dir
echo -e "correctedReads="$correctedReads
echo -e "name="$name
echo -e "reference="$reference
echo -e "illuminaReads="$illuminaReads
echo -e "cores="$cores
# ICORN2_THREADS=$(expr $cores / 2); export ICORN2_THREADS
ICORN2_THREADS=$cores; export ICORN2_THREADS
SMALT_PARAMETER="-n "$cores; export SMALT_PARAMETER
echo -e "cores_iCORN2="$ICORN2_THREADS
echo -e "seqs_circl_1="$seqs_circl_1
echo -e "seqs_circl_2="$seqs_circl_2
echo -e "number_iterations_icorn="$number_iterations_icorn
echo -e "contigs_threshold_size="$contigs_threshold_size
echo -e "InsertsizeRange="$InsertsizeRange
echo -e "taxonid="$taxonid
echo -e "The telomere sequences used are:\nLeft:\t" $telomere_seq_1"\nRight:\t" $telomere_seq_2


##### Checking the installed software:
# PATH with paths to all tools must be properly set: export PATH=$PATH:... in the bashrc or bash_profile files in HOME directory
echo -e "\nPlease be aware that many scripts (bash, perl...) within the ILRA folder are used, and you may need to manually change the interpreter in the corresponding shebang (first line #!) statements so everything works in your system. Many software and dependencies are also required and are being automatically checked. The pipeline will exit if any required software is not found in the variable PATH of your system, which you likely need to change accordingly... You may also need to make scripts executable (chmod) and to make source or export so that the PATH variable and others are available for all scripts"
echo -e "\nIf ILRA keeps working, congrats, all the required software is available"
type ILRA.runSMALT_ver2.sh >/dev/null 2>&1 || { echo >&2 "I require the ILRA folder and all files within to be in the PATH. Aborting..."; exit 1; }
type formatdb >/dev/null 2>&1 || { echo >&2 "I require formatdb but it's not installed or available in the PATH. Aborting..."; exit 1; }
type megablast >/dev/null 2>&1 || { echo >&2 "I require megablast but it's not installed or available in the PATH. Aborting..."; exit 1; }
type samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed or available in the PATH. Aborting..."; exit 1; }
type smalt >/dev/null 2>&1 || { echo >&2 "I require smalt but it's not installed or available in the PATH. Aborting..."; exit 1; }
type abacas2.nonparallel.sh >/dev/null 2>&1 || { echo >&2 "I require ABACAS2 and dependencies but something is not installed or available in the PATH. Aborting..."; exit 1; }
type fastaq >/dev/null 2>&1 || { echo >&2 "I require Fastaq but it's not installed or available in the PATH. Aborting..."; exit 1; }
type icorn2.serial_bowtie2.sh >/dev/null 2>&1 || { echo >&2 "I require iCORN2 and dependencies but something is not installed or available in the PATH. Aborting..."; exit 1; }
type bowtie2 >/dev/null 2>&1 || { echo >&2 "I require Bowtie2 but it's not installed or available in the PATH. Aborting..."; exit 1; }
type java >/dev/null 2>&1 || { echo >&2 "I require java, particularly v1.7.0 for iCORN2, but it's not installed or available in the PATH. Aborting..."; exit 1; }
type minimap2 >/dev/null 2>&1 || { echo >&2 "I require minimap2 but it's not installed or available in the PATH. Aborting..."; exit 1; }
type circlator >/dev/null 2>&1 || { echo >&2 "I require Circlator but it's not installed or available in the PATH. Aborting..."; exit 1; }
type gtool.py >/dev/null 2>&1 || { echo >&2 "I require gtool.py for getting GC content but it's not installed or available in the PATH. Aborting..."; exit 1; }
type assembly-stats >/dev/null 2>&1 || { echo >&2 "I require assembly-stats but it's not installed or available in the PATH. Aborting..."; exit 1; }
type fastq_info.sh >/dev/null 2>&1 || { echo >&2 "I require fastq_info.sh but it's not installed or available in the PATH. Aborting..."; exit 1; }
type pigz >/dev/null 2>&1 || { echo >&2 "I require pigz but it's not installed or available in the PATH. Aborting..."; exit 1; }
type fastqc >/dev/null 2>&1 || { echo >&2 "I require fastqc but it's not installed or available in the PATH. Aborting..."; exit 1; }
type quast.py >/dev/null 2>&1 || { echo >&2 "I require quast but it's not installed or available in the PATH. Aborting..."; exit 1; }
type makeblastdb >/dev/null 2>&1 || { echo >&2 "I require makeblastdb but it's not installed or available in the PATH. Aborting..."; exit 1; }
type blastn >/dev/null 2>&1 || { echo >&2 "I require blastn but it's not installed or available in the PATH. Aborting..."; exit 1; }
type bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools but it's not installed or available in the PATH. Aborting..."; exit 1; }


if [ $mode == "taxon" ] || [ $mode == "both" ] ; then
	databases=$(dirname $0)/databases; mkdir -p $databases
	##### Checking the required software for decontamination steps and the installed databases:
	type centrifuge >/dev/null 2>&1 || { echo >&2 "I require centrifuge but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type retaxdump >/dev/null 2>&1 || { echo >&2 "I require recentrifuge but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type rcf >/dev/null 2>&1 || { echo >&2 "I require recentrifuge but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type rextract >/dev/null 2>&1 || { echo >&2 "I require recentrifuge but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type vecscreen >/dev/null 2>&1 || { echo >&2 "I require vecscreen but it's not installed or available in the PATH. Aborting..."; exit 1; }
	type VSlistTo1HitPerLine.awk >/dev/null 2>&1 || { echo >&2 "I require VSlistTo1HitPerLine.awk but it's not installed or available in the PATH. Aborting..."; exit 1; }

	echo -e "\nPlease note some local databases for the decontamination step (Centrifuge and blast according to DDBJ/ENA/Genbank requirements) are needed"
	echo -e "These databases are large, so please be aware that the RAM memory usage at the step 6 of ILRA may reach hundreds of GBs (100-150GB for P. falciparum, more depending on the genome assembly size)"
	echo -e "ILRA is now going to give you instructions so the databases are downloaded and placed in the corresponding folders (main folder where you have placed ILRA folder, under the directory databases that has been automatically created). ILRA will exit until these steps are performed:"
	echo -e "The NCBI nucleotide non-redundant sequences database (64GB) has to be downloaded and uncompressed by the user"
	echo -e "It can be downloaded from NCBI or from the Centrifuge's webpage. The commands would be: cd /path/ILRA_folder/databases/ && wget https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz && tar -xvzf nt_2018_3_3.tar.gz"
	echo -e "Alternatively, please execute: cd /path/to/ILRA/databases/ && wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz && tar -xvzf nt.*.tar.gz"
	echo -e "The databases names.dmp and nodes.dmp has to be downloaded by the user executing the command from Recentrifuge: cd /path/ILRA_folder/databases/ && retaxdump"
	echo -e "Alternatively, please execute: mkdir -p /path/to/ILRA/databases/taxdump && cd /path/to/ILRA/databases/taxdump && wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip && unzip taxdmp.zip"
	if [[ -f $databases/nt.1.cf ]] && [[ -f $databases/nt.2.cf ]] && [[ -f $databases/nt.3.cf ]] && [[ -f $databases/nt.4.cf ]] && [[ -f $databases/taxdump/names.dmp ]] && [[ $databases/taxdump/nodes.dmp ]] ; then
	  echo -e "Good, ILRA is detecting all of the required databases "
	else
	  echo -e "ILRA is not detecting the required databases to decontaminate and you are not in the light mode, so the pipeline is exiting. Please double check the instructions printed in the log"
		exit 1
	fi
elif [ $mode == "blast" ] || [ $mode == "both" ] ; then
	echo -e "Several databases for conforming to DDBJ/ENA/Genbank requirements are needed, please execute:"
	echo -e "cd /path/to/ILRA/databases/"
	echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz && pigz -d -k -c -p $cores contam_in_euks.fa.gz | makeblastdb -in - -dbtype nucl"
	echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_prok.fa && makeblastdb -in contam_in_prok.fa -dbtype nucl"
	echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_euks.fa && formatdb -p F -i adaptors_for_screening_euks.fa"
	echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_proks.fa && formatdb -p F -i adaptors_for_screening_proks.fa"
	echo -e "wget https://ftp.ncbi.nlm.nih.gov/blast/db/mito.tar.gz && tar -xvzf mito.tar.gz"
	echo -e "wget https://ftp.ncbi.nlm.nih.gov/pub/kitts/rrna.gz && pigz -d -k -c -p $cores rrna.gz | makeblastdb -in - -dbtype nucl"
	if [[ $databases/contam_in_euks.fa ]] && [[ $databases/contam_in_prok.fa ]] && [[ $databases/adaptors_for_screening_euks.fa ]] && [[ $databases/adaptors_for_screening_proks.fa ]] && [[ $databases/mito.ndb ]] && [[ $databases/taxdb.btd ]] && [[ $databases/rrna ]] ; then
	  echo -e "Good, ILRA is detecting all of the required databases "
	else
	  echo -e "ILRA is not detecting the required databases to decontaminate and you are not in the light mode, so the pipeline is exiting. Please double check the instructions printed in the log"
		exit 1
	fi
fi


##### ILRA execution STEPS:
# 1.  Discard smaller contigs
# 2.  megablast -F F
# 2a. Delete contained contigs
# 2b. Find overlaps
# 3.  ABACAS2 with overlap checking (resolve trivial gaps, reorder contigs, get the chromosome names, rename the sequences...)
# 4.  iCORN2 - error correction
# 5.  Circularization or organelles or any sequence reqired
# 6.  Decontamination. Taxonomic classification (Centrifuge/recentrifuge)/final filtering for databases upload
# 7.  Rename sequences, evaluate the assemblies, get telomere sequences counts, GC stats, sequencing depth, converting files...

#### 1. Discard contigs smaller than a threshold:
echo -e "\n\nSTEP 1: Size filtering starting..."; echo -e "Current date/time: $(date)\n"
echo -e "### Excluded contigs based on length threshold: (ILRA.removesmalls.pl)" > Excluded.contigs.fofn
perl -S ILRA.removesmalls.pl $contigs_threshold_size $assembly | sed 's/|/_/g' > 01.assembly.fa
mv Excluded.contigs.fofn ../Excluded.contigs.fofn
echo -e "\n\nSTEP 1: DONE"; echo -e "Current date/time: $(date)\n"


#### 2. MegaBLAST
echo -e "\n\nSTEP 2: MegaBLAST starting..."; echo -e "Current date/time: $(date)\n"
mkdir -p $dir/2.MegaBLAST; cd $dir/2.MegaBLAST
formatdb -p F -i $dir/1.Filtering/01.assembly.fa
megablast -W 40 -F F -a $cores -m 8 -e 1e-80 -d $dir/1.Filtering/01.assembly.fa -i $dir/1.Filtering/01.assembly.fa | awk '$3>98 && $4>500 && $1 != $2' > comp.self1.blast
ILRA.addLengthBlast.pl $dir/1.Filtering/01.assembly.fa $dir/1.Filtering/01.assembly.fa comp.self1.blast &> /dev/null
# 2a. Delete contained contigs
# We want the query to be always the smaller one
echo -e "\n\nSTEP 2a: Delete contained contigs starting..."; echo -e "Current date/time: $(date)\n"
awk '$3>99 && $4>500 && $13 < $14' comp.self1.blast.length | ILRA.getOverlap.pl > 02.ListContained.txt
cat 02.ListContained.txt | awk '$4>90' | cut -f 1 > List.Contained.fofn
echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Excluded contigs that are contained in others: (MegaBLAST/ILRA.getOverlap.pl/ILRA.deleteContigs.pl)" >> ../Excluded.contigs.fofn
cat List.Contained.fofn >> ../Excluded.contigs.fofn
# Now delete the contigs from the MegaBLAST first output, also filter for just hits > 2kb > 99%
echo "Filtering the MegaBLAST output for hits > 2Kb and >99% identity. Please change manually within the pipeline (section 2a) these parameters if needed"
ILRA.deleteContigs.pl List.Contained.fofn $dir/1.Filtering/01.assembly.fa 02.assembly.fa
cat comp.self1.blast.length | awk '$3> 99 && $4 > 2000' | ILRA.deleteEntryinBlast.pl List.Contained.fofn > Blast.merge.blast.length
# 2b. Find overlaps
# Get a bam files to do the merge
echo -e "\n\nSTEP 2b: Find overlaps starting..."; echo -e "Current date/time: $(date)\n"
if [ -f $illuminaReads\_1.fastq ]; then
	echo -e "SMALT parameters are: k-mer size=20, step-size=3, insert size range="$InsertsizeRange". Please check SMALT help and change manually within the pipeline (section 2b) these parameters if needed"
	ILRA.runSMALT_ver2.sh 02.assembly.fa 20 3 $illuminaReads\_1.fastq $illuminaReads\_2.fastq first $InsertsizeRange $cores &> ILRA.runSMALT_ver2.sh_log_out.txt
	echo -e "Check out the log of ILRA.runSMALT_ver2.sh in the file ILRA.runSMALT_ver2.sh_log_out.txt"
# Find overlaps
	ILRA.findoverlaps_ver3.pl Blast.merge.blast.length first.bam 02.assembly.fa OUT 1> ILRA.findoverlaps_ver3.pl_log.txt 2> ILRA.findoverlaps_ver3.pl_log_perl_warnings.txt
	echo -e "\nCheck out the log of ILRA.findoverlaps_ver3.pl in the files log_OUT.txt and results_OUT.txt"
# Save the contigs
	echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Contigs not covered: (ILRA.findoverlaps_ver3.pl)" >> ../Excluded.contigs.fofn
	cat notcovered_OUT.fasta | awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' | grep ">" >> ../Excluded.contigs.fofn
	echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### New contigs merging overlapping contigs(ILRA.findoverlaps_ver3.pl):" >> ../Excluded.contigs.fofn
	if [  $(cat log_OUT.txt | grep ^sequences: | awk '{ gsub("sequences: " , "") ; print }') -eq 0 ] ; then
		echo "No filtering ILRA.findoverlaps_ver3.pl"
	else
		cat results_OUT.txt >> ../Excluded.contigs.fofn
	fi
# Cleaned assembly
	cp mergedseq_OUT.fasta 03.assembly.fa
else
# Bypass if Illumina reads not provided
	echo -e "Illumina reads not provided and no filtering by ILRA.findoverlaps_ver3.pl\n"
	cp 02.assembly.fa 03.assembly.fa
fi
echo -e "\n\nSTEP 2: DONE"; echo -e "Current date/time: $(date)\n"


#### 3. ABACAS2
echo -e "\n\nSTEP 3: ABACAS2 starting..."; echo -e "Current date/time: $(date)\n"
mkdir -p $dir/3.ABACAS2; cd $dir/3.ABACAS2
if [ "$doAbacas2" -eq 1 ] ; then
	ABA_CHECK_OVERLAP=0; export ABA_CHECK_OVERLAP; Min_Alignment_Length=1000; Identity_Cutoff=98; doBLAST_for_ACT_inspection=1
	echo "ABACAS2 parameters are: ABA_CHECK_OVERLAP=0, Min_Alignment_Length=1000, Identity_Cutoff=98. Please check ABACAS2 help and change manually within the pipeline (section 3) these parameters if needed"
	abacas2.nonparallel.sh $reference ../2.MegaBLAST/03.assembly.fa $Min_Alignment_Length $Identity_Cutoff $doBLAST_for_ACT_inspection 1> abacas_log_out.txt 2> abacas_log_warnings_errors.txt
	echo -e "\nCheck out the log of abacas2.nonparallel.sh in the files abacas_log_out.txt and abacas_log_warnings_errors.txt"
# Break  and delete N's
	fastaq trim_Ns_at_end Genome.abacas.fasta 03b.assembly.fa
# Save the contigs
	# echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Overlapping contigs_2 (ABACAS2)" >> ../Excluded.contigs.fofn # Uncomment if ABA_CHECK_OVERLAP=1
	# zcat *.overlap_report.gz >> ../Excluded.contigs.fofn # Uncomment if ABA_CHECK_OVERLAP=1
	echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Final sequences not mapped to the reference: (ABACAS2)" >> ../Excluded.contigs.fofn
	cat Res.abacasBin.fna | grep ">" >> ../Excluded.contigs.fofn
	echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Final sequences mapped to the reference: (ABACAS2)" >> ../Excluded.contigs.fofn
	for i in $(ls | grep .contigs.gff); do
		echo "# "$(echo $i | sed 's/^[^.]*.\([^.]*\)..*/\1/')":" >> ../Excluded.contigs.fofn
		cat $i | grep contig | awk '{ print $9 }' | sed 's/^[^"]*"\([^"]*\)".*/\1/' | sed 's/.*=//' >> ../Excluded.contigs.fofn
	done
else
# Bypass if Illumina reads not provided
	echo -e "Illumina reads or reference genome not provided and ABACAS2 not executed\n"
	ln -s 03.assembly.fa 03b.assembly.fa
fi
echo -e "\n\nSTEP 3: DONE"; echo -e "Current date/time: $(date)\n"


#### 4. iCORN2
echo -e "\n\nSTEP 4: iCORN2 starting..."; echo -e "Current date/time: $(date)\n"
mkdir -p $dir/4.iCORN2; cd $dir/4.iCORN2
if [ -f $illuminaReads\_1.fastq ]; then
	echo -e "Be aware that a particular version of Java, 1.7.0, is required for iCORN2. Users must look for and install 'jdk1.7.0_80'. Any iCORN2 error may be due to incorrect Java version, including newer\n"
	iCORN2_fragmentSize=500
	echo "The iCORN2 fragment size used is iCORN2_fragmentSize="$iCORN2_fragmentSize. "Please check iCORN2 help and change manually within the pipeline (section 4) if needed"
	echo -e "Check out the log of icorn2.serial_bowtie2.sh in the files icorn2.serial_bowtie2.sh_log_out.txt and ../7.Stats/07.iCORN2.final_corrections.results.txt. But a preview of the corrections is:"
	icorn2.serial_bowtie2.sh $illuminaReads $iCORN2_fragmentSize $dir/3.ABACAS2/03b.assembly.fa 1 $number_iterations_icorn &> icorn2.serial_bowtie2.sh_log_out.txt
# Cleaned assembly:
	ln -s ICORN2.03b.assembly.fa.[$(expr $number_iterations_icorn + 1)] 04.assembly.fa
# Summary of iCORN2 results:
	mkdir -p $dir/7.Stats; icorn2.collectResults.pl $PWD > ../7.Stats/07.iCORN2.final_corrections.results.txt
	echo -e "### Total SNP: " >> ../7.Stats/07.iCORN2.final_corrections.results.txt; cat ../7.Stats/07.iCORN2.final_corrections.results.txt | awk '{sum+=$6;} END{print sum;}' >> ../7.Stats/07.iCORN2.final_corrections.results.txt
	echo -e "### Total INS: " >> ../7.Stats/07.iCORN2.final_corrections.results.txt; cat ../7.Stats/07.iCORN2.final_corrections.results.txt | awk '{sum+=$7;} END{print sum;}' >> ../7.Stats/07.iCORN2.final_corrections.results.txt
	echo -e "### Total DEL: " >> ../7.Stats/07.iCORN2.final_corrections.results.txt; cat ../7.Stats/07.iCORN2.final_corrections.results.txt | awk '{sum+=$8;} END{print sum;}' >> ../7.Stats/07.iCORN2.final_corrections.results.txt
	echo -e "### Total HETERO: " >> ../7.Stats/07.iCORN2.final_corrections.results.txt; cat ../7.Stats/07.iCORN2.final_corrections.results.txt | awk '{sum+=$9;} END{print sum;}' >> ../7.Stats/07.iCORN2.final_corrections.results.txt
	cat ../7.Stats/07.iCORN2.final_corrections.results.txt
else
# Bypass if Illumina reads not provided
	echo -e "Illumina reads not provided and iCORN2 not executed\n"
	ln -s $dir/3.ABACAS2/03b.assembly.fa 04.assembly.fa
fi
echo -e "\n\nSTEP 4: DONE"; echo -e "Current date/time: $(date)\n"


#### 5. Circlator for organelles
echo -e "\n\nSTEP 5: Circlator starting..."; echo -e "Current date/time: $(date)\n"
mkdir -p $dir/5.Circlator; cd $dir/5.Circlator
if grep -q -E "$seqs_circl_1|$seqs_circl_2" $dir/4.iCORN2/04.assembly.fa; then
# Map the corrected reads
	if [ long_reads_technology="pacbio" ]; then
		minimap2 -x map-pb -H -t $cores -d minimap2_index_pacbio.mmi $dir/4.iCORN2/04.assembly.fa &> mapping_corrected_reads_log_out.txt
		minimap2 -x map-pb -t $cores -a minimap2_index_pacbio.mmi $correctedReads > Mapped.corrected.04.sam 2>> mapping_corrected_reads_log_out.txt
	else
		minimap2 -x map-ont -t $cores -d minimap2_index_nanopore.mmi $dir/4.iCORN2/04.assembly.fa &> mapping_corrected_reads_log_out.txt
		minimap2 -x map-ont -t $cores -a minimap2_index_nanopore.mmi $correctedReads > Mapped.corrected.04.sam 2>> mapping_corrected_reads_log_out.txt
	fi
# Circlator:
	seq_ids=$(awk '{print $3}' Mapped.corrected.04.sam | grep -E "$seqs_circl_1|$seqs_circl_2" | tail -n +2 | sort | uniq)
	for i in $seq_ids; do
		cat Mapped.corrected.04.sam | grep $i | awk '{ print ">"$1"\n"$10 }' >> ForCirc.reads.fasta
	done
	for i in $seq_ids; do
		samtools faidx $dir/4.iCORN2/04.assembly.fa $i >> ForCirc.Ref.fasta
	done
	awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' ForCirc.reads.fasta > ForCirc.reads_2.fasta
	awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' ForCirc.Ref.fasta > ForCirc.Ref_2.fasta
	echo -e "Check out the log of the mapping of the corrected reads and circlator in the files mapping_corrected_reads_log_out.txt and circlator_log_out.txt"
	circlator all ForCirc.Ref_2.fasta ForCirc.reads_2.fasta Out.Circ --threads $cores &> circlator_log_out.txt
# Delete the plastids/organelles/circular sequences from the current assembly version (04.assembly.fa)
	for i in $seq_ids; do
		echo $i >> List.circular_sequences.fofn
	done
	ILRA.deleteContigs.pl List.circular_sequences.fofn $dir/4.iCORN2/04.assembly.fa 05.assembly.fa
	if [ ! -s $dir/5.Circlator/Out.Circ/06.fixstart.fasta ]; then
# Bypass if circlator failed and didn't end
		echo -e "Circlator has failed. Check its log and output for futher details. One likely reason is that there were some problems with the provided reads, and not enough/even coverage for the required contigs. Bypassing..."
		ln -s -f $dir/4.iCORN2/04.assembly.fa 05.assembly.fa
	else
	cat $PWD/Out.Circ/06.fixstart.fasta >> 05.assembly.fa;
	fi
else
# Bypass if the contigs names provided are not found
	ln -s $dir/4.iCORN2/04.assembly.fa 05.assembly.fa
	echo -e "The sequence identifiers provided for circularization were not found in the contig names. Circlator not executed"
fi
echo -e "\n\nSTEP 5: DONE"; echo -e "Current date/time: $(date)\n"


#### 6. Decontamination/taxonomic classification/final masking and filtering for databases upload
if [ $mode == "taxon" ] || [ $mode == "both" ] ; then
	echo -e "\n\nSTEP 6: Centrifuge and decontamination starting..."; echo -e "Current date/time: $(date)\n"
	mkdir -p $dir/6.Decontamination; cd $dir/6.Decontamination
	echo -e "\nLog of Centrifuge:"
	centrifuge -f -x $databases/nt -U $dir/5.Circlator/05.assembly.fa -p $cores --report-file report.txt -S classification.txt --min-hitlen 100
	cat report.txt
	# Extract contigs classified as different organisms
	rcf -n $databases/taxdump -f classification.txt -o recentrifuge_contamination_report.html -e CSV &> rcf_log_out.txt # Add --sequential if problems with multithreading
	perl -S fasta_to_fastq.pl $dir/5.Circlator/05.assembly.fa ? > 05.assembly.fa.fq # assuming default "fake" quality 30 (? symbol, see https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm)
	source $(echo $PATH | awk '{gsub(/:/,"\n",$0)}1' | egrep ILRA | uniq)/ILRA_exclude_taxons_recentrifuge.txt
	echo -e "Check out the log of Recentrifuge in the files rcf_log_out.txt and rextract_log_out.txt"
	rextract -f classification.txt -i $taxonid $TAXONS_TO_EXCLUDE -n $databases/taxdump -q 05.assembly.fa.fq &> rextract_log_out.txt
	sed -n '1~4s/^@/>/p;2~4p' *.fastq > 06.assembly.fa
	# Save contigs
	echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Excluded contigs that are not recognized by Centrifuge as the species of interests: (Check the output of Centrifuge in Step 6)" >> ../Excluded.contigs.fofn
	comm -23 <(cat $dir/5.Circlator/05.assembly.fa | grep ">" | sort) <(cat 06.assembly.fa | grep ">" | sort) >> ../Excluded.contigs.fofn
	if [ -s 06.assembly.fa ]
	then
		echo -e "Centrifuge seems to have run fine. Please check the report file to assess the contigs that corresponded to contamination. ILRA has helped with this, but it is recommended to run Centrifuge/Recentrifuge on the raw sequencing reads to decontaminate and then reassemble and rerun ILRA, if possible"
		# Renaming and making single-line fasta
		if [ -z "$reference" ]; then
			cat $dir/6.Decontamination/06.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_".$n } else { print }' | ILRA.fasta2singleLine.pl - > ../$name.ILRA.fasta
		else
			cat $dir/6.Decontamination/06.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_with_ref_".$n } else { print }' | ILRA.fasta2singleLine.pl - > ../$name.ILRA.fasta
		fi
	else
		# Bypass if Centrifuge failed and didn't end
		echo -e "\nBypassing decontamination step because Centrifuge has been killed due to RAM usage or it failed due to other reasons..."
		ln -f -s $dir/5.Circlator/05.assembly.fa 06.assembly.fa
	fi
	echo -e "\n\nSTEP 6 Centrifuge: DONE"; echo -e "Current date/time: $(date)\n"
elif [ $mode == "blast" ] || [ $mode == "both" ] ; then
	# Final filtering, masking and reformatting to conform the requirements of DDBJ/ENA/Genbank
	mkdir -p $dir/6.Decontamination/Genbank_upload; cd $dir/6.Decontamination/Genbank_upload
	echo -e "\nPlease note the assembly .ILRA_GenBank.fasta is the final assembly with some conversions such as masking or editing the header of mitochondrial sequences, together with other requirements to ease the upload of the genomes to DDBJ/ENA/Genbank"
	# BLAST to screen the submitted sequences against:
	# 1. Common contaminants database that contains vector sequences, bacterial insertion sequences, E. coli and phage genomes...
	blastn -query ../../$name.ILRA.fasta -db $databases/contam_in_euks.fa -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads $cores | awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' > contam_in_euks_genbank.out
	blastn -query ../../$name.ILRA.fasta -db $databases/contam_in_prok.fa -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads $cores | awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' > contam_in_proks_genbank.out
	# 2. A database of adaptors linkers and primers
	vecscreen -d $databases/adaptors_for_screening_euks.fa -f3 -i ../../$name.ILRA.fasta -o vecscreen_in_euks_genbank
	vecscreen -d $databases/adaptors_for_screening_proks.fa -f3 -i ../../$name.ILRA.fasta -o vecscreen_in_proks_genbank
	VSlistTo1HitPerLine.awk suspect=0 weak=0 vecscreen_in_euks_genbank > vecscreen_in_euks_genbank.out
	VSlistTo1HitPerLine.awk suspect=0 weak=0 vecscreen_in_proks_genbank > vecscreen_in_proks_genbank.out
	# 3. A database of mitochondrial genomes
	blastn -query ../../$name.ILRA.fasta -db $databases/mito -out mito_sequences -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 98.6 -soft_masking true -outfmt 7 -num_threads $cores
	awk '$4>=120' mito_sequences > mito_sequences.out
	# 4. A database of ribosomal RNA genes
	blastn -query ../../$name.ILRA.fasta -db $databases/rrna -task megablast -template_length 18 -template_type coding -window_size 120 -word_size 12 -xdrop_gap 20 -no_greedy -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 1E-9 -gapextend 2 -gapopen 4 -penalty -4 -perc_identity 95 -reward 3 -soft_masking true -outfmt 7 -num_threads $cores | awk '$4>=100' > rrna_genbank.out
	# 5. The chromosomes of unrelated organisms. Foreign organisms are those that belong to a different taxonomic group compared to the organism whose sequences are being screened.
	# This is a requirement by the databases but we skipped it here, because ILRA can alternatively use Centrifuge/Recentrifuge to get rid of contamination of other organisms.
	# Process the result and get the final files:
	awk ' NF>2 {print $1"\t"$7"\t"$8} ' *genbank.out | grep -v "#" | bedtools sort -i - > sequences_from_blast_to_mask.bed
	bedtools maskfasta -fi ../../$name.ILRA.fasta -bed sequences_from_blast_to_mask.bed -fo $name.ILRA_masked.fasta -fullHeader
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $name.ILRA_masked.fasta > $name.ILRA_GenBank.fasta
	for i in $(awk ' NF>2 {print $1} ' mito_sequences.out | grep -v "#" | uniq); do
		sed -i "s/${i}/${i} [location=mitochondrion]/" ../../$name.ILRA_GenBank.fasta
	done
	mv $name.ILRA_GenBank.fasta ../../
	echo -e "\n\nSTEP 6 BLAST: DONE"; echo -e "Current date/time: $(date)\n"
else
	echo -e "\n\nLight mode activated. STEP 6: Skipped"; echo -e "Current date/time: $(date)\n"
	# Renaming and making single-line fasta
	if [ -z "$reference" ]; then
		cat $dir/5.Circlator/05.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_".$n } else { print }' | ILRA.fasta2singleLine.pl - > ../$name.ILRA.fasta
	else
		cat $dir/5.Circlator/05.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_with_ref_".$n } else { print }' | ILRA.fasta2singleLine.pl - > ../$name.ILRA.fasta
	fi
fi


#### 7. Rename sequences, evaluate the assemblies, get telomere sequences counts, GC stats, sequencing depth, converting files...
echo -e "\n\nSTEP 7: Renaming, gathering stats and evaluation starting..."; echo -e "Current date/time: $(date)\n"
mkdir -p $dir/7.Stats; cd $dir/7.Stats
# Evaluating the assemblies:
assembly-stats $assembly | head -n 3
assembly-stats $assembly > 07.assembly_stats_original_correction_ILRA.txt
assembly-stats ../$name.ILRA.fasta | head -n 3
assembly-stats ../$name.ILRA.fasta >> 07.assembly_stats_original_correction_ILRA.txt
# Getting telomere counts:
# Extracting potential telomeres (1Kb)
cat ../$name.ILRA.fasta | perl -nle 'if (/>/){ print } else { $l="Left:"; $left=substr($_,0,1000); print "$l\t$left"}' | tr "\n" "\t" | sed 's/>/\n&/g' | sed '1d' > 07.Telomeres_seq_1kb.txt
cat ../$name.ILRA.fasta | perl -nle 'if (/>/){ print } else { $r="Right:"; $right=substr($_,(length($_)-1000),1000); print "$r\t$right"}' | tr "\n" "\t" | sed 's/>/\n&/g' >> 07.Telomeres_seq_1kb.txt
# Get contigs with the specified sequence repetitions in those regions
echo -e "### Sequences with the left telomere sequence "$telomere_seq_1": (within a 1Kb end region, extracted in 07.Telomeres_seq_1kb.txt)" > 07.TelomerContigs.txt
cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i $telomere_seq_1 | awk '{print $1}' >> 07.TelomerContigs.txt
echo -e "\n### Sequences with the right telomere sequence "$telomere_seq_2": (within a 1Kb end region, extracted in 07.Telomeres_seq_1kb.txt)" >> 07.TelomerContigs.txt
cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i $telomere_seq_2 | awk '{print $1}' >> 07.TelomerContigs.txt
# Summary of telomere repeats statistics
echo -e "# Amount of telomere repeats left: "$(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i -c $telomere_seq_1) > 07.TelomersSummary.txt
echo -e "# Amount of telomere repeats right: "$(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i -c $telomere_seq_2) >> 07.TelomersSummary.txt
echo -e "# Amount of telomere repeats: "$(expr $(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i -c $telomere_seq_1) + $(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i -c $telomere_seq_2)) >> 07.TelomersSummary.txt
contigs_both=$(echo $(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i $telomere_seq_1 | awk '{print $1}') $(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i $telomere_seq_2 | awk '{print $1}') | tr ' ' '\n' | sort | uniq -d)
echo -e "# Amount of telomere repeats at contigs with repeats at both ends: "$(expr $(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep "$contigs_both" | grep -i -c $telomere_seq_1) + $(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep "$contigs_both" | grep -i -c $telomere_seq_2)) >> 07.TelomersSummary.txt
echo -e "\n# Amount of contigs with telomere repeats left: "$(cat 07.Telomeres_seq_1kb.txt | grep 'Left' | grep -i $telomere_seq_1 | awk '{print $1}' | wc -l) >> 07.TelomersSummary.txt
echo -e "# Amount of contigs with telomere repeats right: "$(cat 07.Telomeres_seq_1kb.txt | grep 'Right' | grep -i $telomere_seq_2 | awk '{print $1}' | wc -l) >> 07.TelomersSummary.txt
if [ -z "$contigs_both" ] ; then
	echo -e "# Amount of contigs/potential chromosomes with telomere repeats at both ends: 0" >> 07.TelomersSummary.txt
	echo -e "\n# Contigs/potential chromosomes with telomere repeats at both ends: (Length, bp)" >> 07.TelomersSummary.txt
else
	echo -e "# Amount of contigs/potential chromosomes with telomere repeats at both ends: "$(echo $contigs_both | tr " " "\n" | wc -l) >> 07.TelomersSummary.txt
	echo -e "\n# Contigs/potential chromosomes with telomere repeats at both ends: (Length, bp)" >> 07.TelomersSummary.txt
	awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../$name.ILRA.fasta | tr "\n" "\t" | sed 's/>/\n&/g' | grep "$contigs_both" >> 07.TelomersSummary.txt
fi
# Getting GC stats:
gtool.py -sg ../$name.ILRA.fasta | tee 07.Contigs_GC_size.txt; echo -e "\nIndividual contigs:" >> 07.Contigs_GC_size.txt
for i in $(grep "^>" ../$name.ILRA.fasta | awk 'sub(/^>/, "")'); do
	cat ../$name.ILRA.fasta | gtool.py -sg -c $i - >> 07.Contigs_GC_size.txt
done
sed -i 's/stdin//g' 07.Contigs_GC_size.txt
# Getting sequencing depth and stats:
# https://github.com/wudustan/fastq-info
if [ -f $illuminaReads\_1.fastq ]; then
	fastq_info_2.sh $illuminaReads\_1.fastq $illuminaReads\_2.fastq > 07.fastq_info_depth_IlluminaReads.txt
	illumina_read_length=$(head -n 2 $illuminaReads\_1.fastq |tail -n 2|wc -c)
	echo "Illumina_read_length =" $illumina_read_length
	fastq-info.sh -r $illumina_read_length $illuminaReads\_1.fastq $illuminaReads\_2.fastq ../$name.ILRA.fasta > 07.fastq_info_depth_IlluminaReads_assembly.txt
	fastqc $illuminaReads\_1.fastq $illuminaReads\_2.fastq -o $PWD --noextract -q -t $cores
fi
# Comparing with reference genes:
echo -e "Running QUAST. Please be aware that for providing reference genes, a GFF file with gene or operon as feature type field, or a bed file (sequence name, start position, end position, gene id) are accepted"
if [ -z "$reference" ] ; then
	echo -e "Running QUAST without reference..."
	quast.py ../$name.ILRA.fasta --threads $cores --labels $name --space-efficient --eukaryote &> quast.py_log_out.txt
else
	if [ -f $illuminaReads\_1.fastq ]; then
		if [ -z "$gff_file" ] ; then
			echo -e "Running QUAST with reference and structural variants calling and processing mode..."
			quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --space-efficient --pe1 $illuminaReads\_1.fastq --pe2 $illuminaReads\_2.fastq --eukaryote &> quast.py_log_out.txt
		else
			echo -e "Running QUAST with reference, with gff file and with structural variants calling and processing mode..."
			quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --space-efficient --pe1 $illuminaReads\_1.fastq --pe2 $illuminaReads\_2.fastq --eukaryote &> quast.py_log_out.txt
		fi
	else
		if [ -z "$gff_file" ] ; then
			echo -e "Running QUAST with reference..."
			quast.py ../$name.ILRA.fasta -r $reference --threads $cores --labels $name --space-efficient --eukaryote &> quast.py_log_out.txt
		else
			echo -e "Running QUAST with reference and with gff file..."
			quast.py ../$name.ILRA.fasta -r $reference -g $gff_file --threads $cores --labels $name --space-efficient --eukaryote &> quast.py_log_out.txt
		fi
	fi
fi
# Converting files to minimize space
if [ -f $dir/2.MegaBLAST/first.bam ]; then
	samtools view -@ $cores -T $dir/1.Filtering/01.assembly.fa -C -o $dir/2.MegaBLAST/first.bam.cram $dir/2.MegaBLAST/first.bam &> $dir/2.MegaBLAST/first.bam.cram_log_out.txt
fi
if [ -f $dir/4.iCORN2/ICORN2_1/out.sorted.markdup.bam ]; then
	cd $dir/4.iCORN2
	for f in $( find . -type f -name "*.bam" ); do
		samtools view -@ $cores -T $dir/3.ABACAS2/03b.assembly.fa -C -o "${f}".cram "${f}" &> $dir/4.iCORN2/"${f}.cram_log_out.txt"
	done
fi
if [ -f $dir/5.Circlator/Mapped.corrected.04.sam ]; then
	samtools view -@ $cores -T $dir/4.iCORN2/04.assembly.fa -C -o $dir/5.Circlator/Mapped.corrected.04.sam.cram $dir/5.Circlator/Mapped.corrected.04.sam &> $dir/5.Circlator/Mapped.corrected.04.sam.cram_log_out.txt
	samtools view -@ $cores -T $dir/5.Circlator/ForCirc.Ref_2.fasta -C -o $dir/5.Circlator/Out.Circ/01.mapreads.bam.cram $dir/5.Circlator/Out.Circ/01.mapreads.bam &> $dir/5.Circlator/01.mapreads.bam.cram_log_out.txt
fi
echo -e "\nIf needed, for converting compressed .cram files back to .bam apply the command: samtools view -@ $cores -T filename.fasta -b -o output.bam input.cram"
# Cleaning up:
cd $dir; rm $(find . -regex ".*\.\(bam\|sam\)"); rm $illuminaReads\_1.fastq; rm $illuminaReads\_2.fastq

echo -e "\n\nSTEP 7: DONE"; echo -e "Current date/time: $(date)\n"


echo -e "\nILRA IS DONE"; echo -e "Current date/time: $(date)\n"
end=`date +%s`; runtime=$((end-start))
echo -e "\nCores: $cores"
echo -e "Final runtime (secs): $runtime"
