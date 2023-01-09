# ILRA
Improvement of Long Read Assemblies (ILRA) is a pipeline to help in the post-assembly process (finishing of genomes) by cleaning and merging contigs, correcting homopolymer tracks and circularizing plastids. 

## Installation
ILRA is based on several standard tools and novel scripts. We suggest three different alternatives for installation. Please choose one of:


1) Users without a bioinformatics setup can use our Linux virtual machine (https://q-r.to/ILRA_VM). You will need to install VirtualBox (https://www.oracle.com/virtualization/technologies/vm/downloads/virtualbox-downloads.html), set up Ubuntu x64 and mount the downloaded disc (.vdi). See the file 'VM.install.pdf' for further help. The username is 'bioinfo' and the password 'Glasgow2020'.
Please be aware that the ILRA version within the VM may be outdated. Please double check and update if necessary following option 2.


2) The fastest option is to use the folder 'external_software', which contain some of the required software, a script to install dependencies (mainly through miniconda), and a suggestion of the PATH to be set. To install and setup everything required to run ILRA, please execute:
```
git clone https://github.com/ThomasDOtto/ILRA
cd ILRA/external_software/ILRA/
bash ILRA_installation.sh 2>&1 | tee ILRA_installation.sh.log # Check log to ensure successful installation of dependencies
```

This should work if you already have miniconda3 installed, and also install miniconda3 if not available. Plese keep in mind that in the ILRA_installation.sh script most of the versions of the tools installed by conda are not frozen, so please open an issue or you may need to manually change installed versions if conda fails with new dependencies problems.

If not executing the light mode that skips the decontamination step in ILRA (by default, or argument '-m light'), many databases will be required and must be also installed. ILRA is going to try and do it automatically, or to provide instructions after a first execution.


3) The last and less recommended option is to manually install the required software.
If you want to manually install the software, check out in the file external_software/ILRA/ILRA_installation.sh the list of required tools, which must be in the PATH when running. Please be aware that many scripts (bash, perl...) within the ILRA folder are used, and you may need to manually change the interpreter in the corresponding shebang statements (first line #!) so everything works in your system. The software and dependencies are automatically checked and the pipeline will exit if any required software is not found in the variable PATH of your system, which you likely need to change accordingly. You may also need to make scripts executable ('chmod' command) and to make source or export so that the PATH variable and others are available for all scripts.


## Quick start
```
source ILRA/external_software/ILRA/path_to_source # To set up the PATH if you have followed option 2 for installation above
# 'Light' mode skipping decontamination:
ILRA.sh -a $PWD/test_data/assembly_Pf_test.fasta -o $PWD/test_data/out_ILRA_test -c $PWD/test_data/corrected_reads_Pf_test_subset.fastq.gz -n test -r $PWD/test_data/PlasmoDB-47_Pfalciparum3D7_Genome_core_PMID_29862326.fasta -I $PWD/test_data/Illumina_short_reads_Pf_test_subset -t 4 -g $PWD/test_data/PlasmoDB-50_Pfalciparum3D7.gff -L pb -q no 2>&1 | tee -a $PWD/test_data/out_ILRA_test_log.txt
# 'Both' mode with decontamination based on kraken2 and blast:
ILRA.sh -a $PWD/test_data/assembly_Pf_test.fasta -o $PWD/test_data/out_ILRA_test_m_both -c $PWD/test_data/corrected_reads_Pf_test_subset.fastq.gz -n test -r $PWD/test_data/PlasmoDB-47_Pfalciparum3D7_Genome_core_PMID_29862326.fasta -I $PWD/test_data/Illumina_short_reads_Pf_test_subset -t 4 -g $PWD/test_data/PlasmoDB-50_Pfalciparum3D7.gff -L pb -m both -q no 2>&1 | tee -a $PWD/test_data/out_ILRA_test_m_both_log.txt
```
The test run will take around 5 minutes in 'light' mode and around 10 minutes in 'both' mode using 4 cores.

Please go through the output file 'out_ILRA_test_log.txt' to get the details on the pipeline processing steps and final output. For the test run, the last step of quality control (using QUAST, BUSCO, gtools, telomere analyses...) is disconnected by using '-q no'.



## ILRA arguments
Please refer to the help page for futher details:
```
ILRA.sh -h
```
Parameters are not positional. If you did not provide a required parameter, the pipeline may exit or use default values if possible (check the help, the log after execution, or the 'Arguments / Variables' section in the ILRA main script 'ILRA.sh').

In general, from an assembly as input (argument '-a'), ILRA is going to provide a polished assembly as output (the file 'Name.ILRA.fasta').

Please do provide or not the arguments '-C' and '-I' to indicate whether to use short reads to perform error correction (iCORN2) and to find overlapped contigs (ILRA.findoverlaps_ver3.pl).

Depending on whether you provided a reference genome (argument '-r'), reordering and renaming of the contigs (ABACAS2) are going to be skipped and assessment by QUAST would be run without the reference. Similarly, the availability of a reference annotation (argument '-g') would determine the mode to run QUAST. The debug mode (argument '-d' makes possible to resumen the execution of ILRA from a particular step). The argument '-p' determine whether Pilon is used for short reads correction (default 'no') and the argument '-q' determines whether a final extra step for assessing the quality and completeness of the corrected assembly (QUAST, BUSCO, gathering sequences, analyzing telomeres...) is included (default 'yes').

More arguments... XXX

Finally, ILRA can be run in alternative modes (argument '-m'): 
* '-m taxon': To perform decontamination based on taxonomic classification, which would be more computationally expensive.
* '-m blast': To perform decontamination and formatting for online submission based on blasting against databases, which would be less computationally expensive.
* '-m both': To perform both. 
* '-m light': To skip decontamination and expedite the process (default if argument not provided).



## Comments
We used ILRA to improve many genomes. The novel Plasmodium de novo assemblies included in the article are in the folder 'PlasmodiumGenomes', in Zenodo (doi:10.5281/zenodo.7516750) and NCBI (accession:XXX).

ILRA is a continuation of the IPA project (https://github.com/ThomasDOtto/IPA). We felt to rename the script as our updated version works with every long-read technology.
  
Please cite this reference and our Zenodo tag when using ILRA for your publications:

> From contigs to chromosomes: automatic Improvement of Long Read Assemblies (ILRA)
> 
> José L Ruiz, Susanne Reimering, Mandy Sanders, Juan David Escobar-Prieto, Nicolas M. B. Brancucci, Diego F. Echeverry, Abdirahman I. Abdi, Matthias Marti, Elena Gomez-Diaz, Thomas D. Otto
> 
> bioRxiv 2021.07.30.454413; doi: https://doi.org/10.1101/2021.07.30.454413
```
@article {Ruiz2021.07.30.454413,
	author = {Ruiz, Jos{\'e} L and Reimering, Susanne and Sanders, Mandy and Escobar-Prieto, Juan David and Brancucci, Nicolas M. B. and Echeverry, Diego F. and Abdi, Abdirahman I. and Marti, Matthias and Gomez-Diaz, Elena and Otto, Thomas D.},
	title = {From contigs to chromosomes: automatic Improvement of Long Read Assemblies (ILRA)},
	elocation-id = {2021.07.30.454413},
	year = {2021},
	doi = {10.1101/2021.07.30.454413},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Recent advances in long read technologies not only enable large consortia to aim to sequence all eukaryotes on Earth, but they also allow many laboratories to sequence their species of interest. Although there is a promise to obtain {\textquoteright}perfect genomes{\textquoteright} with long read technologies, the number of contigs often exceeds the number of chromosomes significantly, containing many insertion and deletion errors around homopolymer tracks. To overcome these issues, we implemented ILRA to correct long reads-based assemblies, a pipeline that orders, names, merges, and circularizes contigs, filters erroneous small contigs and contamination, and corrects homopolymer errors with Illumina reads. We successfully tested our approach to assemble the genomes of four novel Plasmodium falciparum samples, and on existing assemblies of Trypanosoma brucei and Leptosphaeria spp. We found that correcting homopolymer tracks reduced the number of genes incorrectly annotated as pseudogenes, but an iterative correction seems to be needed to reduce high numbers of homopolymer errors. In summary, we described and compared the performance of a new tool, which improves the quality of long read assemblies. It can be used to correct genomes of a size of up to 300 Mb.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2021/08/01/2021.07.30.454413},
	eprint = {https://www.biorxiv.org/content/early/2021/08/01/2021.07.30.454413.full.pdf},
	journal = {bioRxiv}
}
```


