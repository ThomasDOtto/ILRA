# ILRA
Improvement of Long Read Assemblies (ILRA) is a pipeline to help in the post-assembly process (finishing of genomes) with steps such as cleaning and merging contigs, correcting homopolymer tracks, circularizing plastids and quality control.

## Installation
ILRA is based on several standard tools and novel scripts. We suggest two alternatives for installation. Please choose one of:


1) The fastest option is to use the folder 'external_software', which contain some of the required software, a script to install dependencies (mainly through miniconda and downloading binaries), and a suggestion of the PATH to be set. To install and setup everything required to run ILRA, please execute:
```
git clone https://github.com/ThomasDOtto/ILRA
cd ILRA/external_software/ILRA/
bash ILRA_installation.sh 2>&1 | tee ILRA_installation.sh.log # Check log to ensure successful installation of dependencies
```
This should work systems-wide if you already have miniconda3 installed, and also install miniconda3 in the same directory if not available. Please make sure that your PATH is as clean as possible and that you have not activated other conda environment before installation (i.e. you are in the 'base' environment or without that conda feature activated). The versions of the tools installed by conda are frozen ('.yml' files), so please do not install other versions to avoid conflicts. You may need to manually change installed versions or please open an [issue](https://github.com/ThomasDOtto/ILRA/issues) for any installation-related aspect. The full instalation above, including 'git clone' of the repository, took ~40 minutes.

If not executing the 'light' mode that skips the decontamination step in ILRA (by default, or argument '-m light'), many databases will be required and must be also installed. ILRA is going to try and do it automatically, and instructions are going to be provided after a first execution.

After installation and prior any execution, please make sure that you have not activated any other environment, and manually export the PATH required by executing:
```
source ILRA/external_software/ILRA/path_to_source
```


2) The less recommended option is to manually install the required software.
If you want to manually install the software, check out the files 'external_software/ILRA/ILRA_installation.sh' and in the '.yml' files (i.e. frozen conda environments) the list of required tools, which must be in the PATH when running. Please be aware that many scripts (bash, perl...) within the ILRA folder are used, and you may need to manually change the interpreter in the corresponding shebang statements (i.e. first line #!) so everything works in your system. The software and dependencies are automatically checked and the pipeline will exit if any required software is not found in the variable PATH of your system, which you likely need to change accordingly. You may also need to make scripts executable ('chmod' command) and to make source or export so that the PATH variable and others are available for all scripts.


## ILRA arguments
Please find a guided minimal example in the [test_data](https://github.com/ThomasDOtto/ILRA/tree/main/test_data) subfolder and a brief ILRA walktrough below.

Parameters are not positional and all are going to be automatically handled by the pipeline (i.e. if you did not provide a required parameter, the pipeline may exit asking for it or use default values if possible). For example, if the argument '-f' is not provided, contigs shorter than 5,000 bp by default will be filtered out. 

In general, from an assembly as input (argument '-a'), ILRA is going to provide a polished assembly as output with the argument '-n' as name (the file 'name.ILRA.fasta', in the folder provided by the argument '-o').

Please do provide or not the arguments '-C', '-c', '-R' and '-I' to indicate whether to use short reads to perform error correction (i.e. iCORN2 / Pilon iteratively, with argument '-i' as number of iterations), and to find and filter out overlapping contigs (if coverage by Illumina short reads is even, argument 'F'). Please do provide the long reads sequencing technology used with the argument '-L'.

Depending on whether you provided a reference genome (argument '-r'), reordering and renaming of the contigs (ABACAS2 and the argument '-B' to perform blasting) is going to be skipped, and assessment by QUAST would be run without the reference. Similarly, the availability of a reference annotation (argument '-g') would determine the mode to run QUAST, or the presence of certain names in the contigs marking the sequences to circularize (arguments '-s' and '-S') would mean that Circlator is executed or not. ILRA would automatically resume any previous attempt or incomplete run, and the debug mode (argument '-d') makes possible to resume the execution of ILRA from a particular step. The argument '-p' determine whether Pilon should be used for short reads correction instead of iCORN2 (default 'no'). The argument '-q' determines whether a rather computationally expensive extra step for assessing the quality and completeness of the corrected assembly (i.e., QUAST, BUSCO, gathering sequences, looking in the telomeres for the sequenes provided by the arguments '-e' and '-E'...) is included (default 'yes'). The argument 'M' is required to control the maximum RAM used, and the argument '-l' activates a 'low memory' mode at the expense of more processing time. The arguments '-b', '-P' and '-A' provide the number of parts to split the sequences and to simultaneously process in ILRA, iCORN2/Pilon and ABACAS2, respectively. The argument '-t' provides the number of cores to be used when multithreading is possible.

Finally, ILRA can be run in alternative modes (argument '-m'): 
* '-m taxon':	To perform decontamination based on taxonomic classification, which would be more computationally expensive.
* '-m blast':	To perform decontamination and formatting for online submission based on blasting against databases, which would be less computationally expensive.
* '-m both':	To perform both. 
* '-m light':	To skip decontamination and expedite the process (default if argument not provided).
The location of the downloaded databases to be used in the decontamination step should be provided with the argument '-D' and the location of the kraken2 database should be provided with the argument '-k'. The taxon id of the organism of interest that should be kept when filtering must be provided with the argument '-T'. A first execution of ILRA will provide comprehensive isntructions to perform the installation of the databases.

Please refer to the help page for futher details:
```
ILRA.sh -h

usage: ILRA.sh [options]
		-h | -help # Type this to get help
		-a | -assembly # Name of the long reads assembly to correct (FASTA format, can be gzipped)
		-f | -filter_contig_size # Size threshold to filter the contigs (bp)
		-F | -Do you have even Illumina reads coverage? If yes contigs in the long reads assembly not covered by Illumina short reads by a certain threshold will be filtered out ('no' /'yes' by default)
		-I | -Illumina_reads # Root name of the paired-end Illumina reads (FASTQ format, must be gzipped and name root_1.fastq.gz and root_2.fastq.gz prior execution)
		-C | -Correction_Illumina_reads # Whether illumina reads are provided and all steps of correction should be performed ('no' /'yes' by default)
		-R | -Range_insert_size # Insert size range of the Illumina reads to use in mapping by SMALT (bp)
		-i | -iterations_iCORN2 # Number of iterations to perform in iCORN2
		-r | -reference # Reference file (full pathway, FASTA format)
		-g | -gff_file # Reference annotation file (full pathway, GFF format)
		-c | -corrected_reads # Corrected long reads for circulatization of the contigs containing the strings by -s and -S (FASTQ format, can be gzipped)
		-s | -seq_circularize_1 # Regex pattern to find in the contig names and circularize
		-S | -Seq_circularize_2 # Regex pattern to find in the contig names and circularize
		-L | -Long_reads_technology # Technology of long reads sequencing (pb/ont)
		-T | -Taxon_ID # NCBI taxon ID or a scientific name
		-e | -ending_telomere_seq_1 # Telomere-associated sequence to search in the first strand
		-E | -Ending_telomere_seq_2 # Telomere-associated sequence to search in the second strand
		-o | -output # Output folder (absolute pathway)
		-n | -name # Base name of the output file
		-t | -threads # Number of cores to use in multithreaded steps
		-d | -debug_step # For debug, step to remove the content of the corresponding folder and resume a failed run ('step1', 'step2a', 'step2b', 'step3', 'step4', 'step4i', 'step5', 'step6', or 'step7')
		-D | -databases # Folder for storing the databases for the decontamination step (by default, 'databases' under ILRA main folder)
		-K | -Kraken2_fast_mode # Kraken2 fast mode, consisting on copying the Kraken2 database to /dev/shm (RAM) so execution is faster ('yes' /'no' by default)
		-k | -Kraken2_databases # Folder within the folder databases (-D) containing the database used by Kraken2 (by default, 'standard_eupathdb_48_kraken2_db')
		-b | -block_size # Block size for parallel processing (by default, 10)
		-p | -pilon # Whether to use pilon instead of iCORN2 ('yes'/'no' by default)
		-P | -parts_icorn2_split # Number of parts to split the input sequences of iCORN2 before processing them (0 by default, which means no splitting)
		-A | -abacas2_split # Number of parts to split and process in parallel in ABACAS2 (by default the argument -b, block_size, but may be necessary to decrease due to memory issues)
		-B | -abacas2_blast # Whether to do blast within ABACAS2 to compare with the reference and display in ACT (1 by default, which means blasting, or 0)
		-q | -quality_assesment # Whether to execute a final step for assessing the quality of the corrected assembly, gathering sequences, analyzing telomeres... etc ('no'/'yes' by default)
		-M | -java_memory # Max Java memory (heap space) to be used ('XXg', by default 200g=200GB used)
		-l | -low_memory # Activate low memory mode for iCORN2 ('yes'/'no' by default)
		-m | -mode # Add 'taxon' to execute decontamination based on taxonomic classification by kraken2, add 'blast' to execute decontamination based on BLAST against databases as requested by the DDBJ/ENA/Genbank submission, add 'both' to execute both approaches, and add 'light' to execute ILRA in light mode and skip these steps (default)
```

## Comments
We used ILRA to improve many genomes. The novel Plasmodium de novo assemblies included in the article are in the folder 'PlasmodiumGenomes' within this repository, in Zenodo ([10.5281/zenodo.7516750](https://doi.org/10.5281/zenodo.7516750)) and NCBI ([PRJNA714074](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA714074) / [PRJNA757237](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA757237)).

ILRA is a continuation of the IPA project (https://github.com/ThomasDOtto/IPA). We felt to rename the tool as our updated version works with every long-read technology.
  
Please cite this [reference](https://doi.org/10.1093/bib/bbad248) and our [Zenodo](https://doi.org/10.5281/zenodo.7516750) tag when using ILRA for your publications:

> Ruiz, J. L., Reimering, S., Escobar-Prieto, J. D., Brancucci, N. M., Echeverry, D. F., Abdi, A. I., ... & Otto, T. D. (2023). From contigs towards chromosomes: automatic improvement of long read assemblies (ILRA). Briefings in Bioinformatics, bbad248.
```
@article{ruiz2023contigs,
  title={From contigs towards chromosomes: automatic improvement of long read assemblies (ILRA)},
  author={Ruiz, Jos{\'e} Luis and Reimering, Susanne and Escobar-Prieto, Juan David and Brancucci, Nicolas MB and Echeverry, Diego F and Abdi, Abdirahman I and Marti, Matthias and G{\'o}mez-D{\'\i}az, Elena and Otto, Thomas D},
  journal={Briefings in Bioinformatics},
  pages={bbad248},
  year={2023},
  publisher={Oxford University Press}
}
```

## Support
ILRA is under active maintenance. Please report any [issue](https://github.com/ThomasDOtto/ILRA/issues) or [contact us](mailto:joseluis.ruiz@csic.es?subject=[GitHub]%20Source%20ILRA%20Support) to request support.
