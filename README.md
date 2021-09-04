# ILRA
Improvement of Long Read Assemblies (ILRA) is a pipeline to help in the post-assembly process (finishing) by cleaning and merging contigs, correcting homopolymer tracks and circularizing plastids. We also have the code preinstalled on a virtual Linux machine.  

## Installation
ILRA is based on several standard tools. 

1) We would recommend users without a bioinformatics setup to use our virtual machine (https://q-r.to/ILRA_VM). You will need to install VirtualBox (https://www.oracle.com/virtualization/technologies/vm/downloads/virtualbox-downloads.html), set up Ubuntu x64 and mount the disc. See the pdf VM.install.pdf for further help. The user is "bioinfo" and the password "Glasgow2020". 

Once logged in, you can run the ILRA installation with test_data by typing:
```
~/ILRA/example.sh
```

Please refer to the /test_data/README file or the help ('ILRA.sh -h') for further details. This line will execute the pipeline:
```
# ...  Check test_data/README or ILRA.sh -h
ILRA.sh -a $ASSEMBLY -o $OUTPUT_FOLDER_ILRA -c $CORRECTED_READS -n subset_test -r $REFERENCE -I $ILLU_READS -t $CORES -g $GFF_REF_FILE -L pb
```

Please note the version of ILRA installed in the virtual machine is likely outdated. You can overwrite the version within the virtual machine with the most recent one by first removing the folder of ILRA and then obtaining the most updated version by running:
```
git clone https://github.com/ThomasDOtto/ILRA
```

2) The fastest option is to use the file "external_software.tar.gz", which contains the precompiled binaries and wrapper scripts to install the software required. To install everything required to run ILRA, please execute:
```
git clone https://github.com/ThomasDOtto/ILRA
cd ILRA
wget --no-check-certificate "http://q-r.to/ILRA_soft" -O external_software.tar.gz; tar -xvzf external_software.tar.gz; rm external_software.tar.gz
source external_software/path_to_source; ./external_software/finish_installation.sh
```

3) The last option is to manually install the required software (please find the list and further details in the INSTALL file).


## Quick start
```
cd ILRA
source external_software/path_to_source
ILRA.sh -a test_data/assembly_Pf_test.fasta -o test_data/out_ILRA_test -c test_data/corrected_reads_Pf_test_subset.fastq.gz -n test -r test_data/PlasmoDB-47_Pfalciparum3D7_Genome_core_PMID_29862326.fasta -I test_data/Illumina_short_reads_Pf_test_subset -t 4 -g test_data/PlasmoDB-50_Pfalciparum3D7.gff -L pb | tee -a test_data/test_log.txt
```
The test run will take around 5-10 minutes ('light' mode) and around 10-20 minutes ('both' mode) using 4 cores.


### ILRA arguments:
```
ILRA.sh -a <Assembly> -o <Results directory> -c <Long reads corrected reads> -n <Name of results> -r <Reference genome for ABACAS2> -I <Root name of Illumina short reads> -t <Number of cores to use> -s <First sequence name to circularize> -S <Second sequence name to circularize> -i <Number of iterations for iCORN2> -f <Size threshold for discarding contigs> -R <Insert size range for Illumina short reads> -T <NCBI taxonomy id to extract> -g <GFF reference genes annotation file> -L <Long reads sequencing technology> -e <Telomeric sequence left> -E <Telomeric sequence right> -m <Execution mode> -C <Perform error correction by short reads> -h <Show help>
```
Parameters are not positional. If you did not provide any required parameter, the pipeline will exit or use default values if possible (check the help, the log after execution or the 'Arguments / Variables' section in the pipeline ILRA.sh).

Please check the help page for futher details:
```
ILRA.sh -h
```


## Comments
We used ILRA to improve many genomes. It is a continuation of the IPA (https://github.com/ThomasDOtto/IPA) project. We felt to rename the script as our updated version works with every long-read technology.
  
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


