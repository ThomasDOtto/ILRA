# ILRA
Improvement of Long Read Assemblies (ILRA) is a pipeline to help in the post-assembly process (finishing) by cleaning and merging contigs, correct homopolymer tracks and circulate plastids.   We also have the code preinstalled on a virtual Linux machine.  

## Installation
ILRA is based on several standard tools. Please find the list of software you need to install in the INSTALL.txt file.

We would recommend users without a bioinformatics setup to use our virtual machine - https://drive.google.com/file/d/1jCxlKYF9A4IQfgmGldV3rt87FugrdK9l/view?usp=sharing You will need to install Virtual Box (https://www.oracle.com/virtualization/technologies/vm/downloads/virtualbox-downloads.html) and mount the disc. See the pdf VM.install.pdf for further help. The user is "bioinfo" and the password "Glasgow2020". 
Once logged in, you can test the ILRA installation by typing:
~/ILRA/example.sh
this will execute the command:
"ILRA.sh -a $ASSEMBLY -o $OUTPUT_FOLDER_ILRA -c $CORRECTED_READS -n subset_test -r $REFERENCE -I $ILLU_READS -t $CORES -g $GFF_REF_FILE -L p
