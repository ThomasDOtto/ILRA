# ILRA
Improvement of Long Read Assemblies (ILRA) is a pipeline to help in the post-assembly process (finishing) by cleaning and merging contigs, correcting homopolymer tracks and circularizing plastids. 

## Installation
ILRA is based on several standard tools and novel scripts. We suggest three different alternatives for installation. Please choose one of:


1) Users without a bioinformatics setup can use our Linux virtual machine (https://q-r.to/ILRA_VM). You will need to install VirtualBox (https://www.oracle.com/virtualization/technologies/vm/downloads/virtualbox-downloads.html), set up Ubuntu x64 and mount the downloaded disc (.vdi). See the file 'VM.install.pdf' for further help. The username is 'bioinfo' and the password 'Glasgow2020'.



2) The fastest option is to use the folder 'external_software', which contains a wrapper script to install the required software and PATH to set. To install and setup everything required to run ILRA, please execute:
```
git clone https://github.com/ThomasDOtto/ILRA
cd ILRA
bash external_software/ILRA/finish_installation.sh | tee -a external_software/ILRA/finish_installation_log.txt
```



3) The last and less recommended option is to manually install the required software (please find the list and further details in the 'INSTALL' file).



## Quick start
```
cd ILRA
source external_software/ILRA/path_to_source # To set up the PATH if you have followed option 2 for installation above
# Light mode:
ILRA.sh -a $PWD/test_data/assembly_Pf_test.fasta -o $PWD/test_data/out_ILRA_test -c $PWD/test_data/corrected_reads_Pf_test_subset.fastq.gz -n test -r $PWD/test_data/PlasmoDB-47_Pfalciparum3D7_Genome_core_PMID_29862326.fasta -I $PWD/test_data/Illumina_short_reads_Pf_test_subset -t 4 -g $PWD/test_data/PlasmoDB-50_Pfalciparum3D7.gff -L pb | tee -a $PWD/test_data/out_ILRA_test_log.txt
# Decontamination based on centrifuge and blast. Both mode:
ILRA.sh -a $PWD/test_data/assembly_Pf_test.fasta -o $PWD/test_data/out_ILRA_test_m_both -c $PWD/test_data/corrected_reads_Pf_test_subset.fastq.gz -n test -r $PWD/test_data/PlasmoDB-47_Pfalciparum3D7_Genome_core_PMID_29862326.fasta -I $PWD/test_data/Illumina_short_reads_Pf_test_subset -t 4 -g $PWD/test_data/PlasmoDB-50_Pfalciparum3D7.gff -L pb -m both | tee -a $PWD/test_data/out_ILRA_test_m_both_log.txt
```
The test run will take around 5 minutes in 'light' mode and around 15 minutes in 'both' mode using 4 cores.

Please go through the output file 'out_ILRA_test_log.txt' to get the details on the pipeline processing and final output. For the test run, QUAST and BUSCO are expected to fail given the input assembly is an artificial subset.



## ILRA arguments
```
ILRA.sh -a <Assembly> -o <Results directory> -c <Long reads corrected reads> -n <Name> -r <Reference genome> -I <Root name of Illumina short reads> -t <Number of cores to use> -s <First sequence name to circularize> -S <Second sequence name to circularize> -i <Number of iterations for iCORN2> -f <Size threshold for discarding contigs> -R <Insert size range for Illumina short reads> -T <NCBI taxonomy id to extract> -g <GFF reference genes annotation file> -L <Long reads sequencing technology> -e <Telomeric sequence left> -E <Telomeric sequence right> -m <Execution mode> -C <Perform error correction by short reads> -d <Step to commence the run> -h <Show help>
```
Parameters are not positional. If you did not provide a required parameter, the pipeline may exit or use default values if possible (check the help, the log after execution, and the 'Arguments / Variables' section in the ILRA main script 'ILRA.sh').

Please refer to the help for futher details:
```
ILRA.sh -h
```
In general, from an assembly as input (argument '-a'), ILRA is going to provide a polished assembly as output (the file 'Name.ILRA.fasta').

Please provide or not the arguments '-C' and '-I' to indicate whether to use short reads to perform error correction (iCORN2) and to find overlapped contigs (ILRA.findoverlaps_ver3.pl).

Depending on whether you provided a reference genome (argument '-r'), reordering and renaming of the contigs (ABACAS2) are going to be skipped and assessment by QUAST would be run without the reference. Similarly, the availability of a reference annotation (argument '-g') would determine the mode to run QUAST. The debug mode (argument '-d' makes possible to resumen the execution of ILRA from a particular step).

Finally, ILRA can be run in alternative modes (argument '-m'): 
* '-m taxon': To perform decontamination based on taxonomic classification, which would be more computationally expensive.
* '-m blast': To perform decontamination and formatting for online submission based on blasting against databases, which would be less computationally expensive.
* '-m both': To perform both. 
* '-m light': To skip decontamination and expedite the process (default if argument not provided).



## Comments
We used ILRA to improve many genomes. The novel Plasmodium de novo assemblies included in the article are in the folder 'PlasmodiumGenomes' and in Zenodo (doi:XXX)

ILRA is a continuation of the IPA (https://github.com/ThomasDOtto/IPA) project. We felt to rename the script as our updated version works with every long-read technology.
  
Please cite this preprint when using ILRA for your publications:

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


