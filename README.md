# ILRA
Improvement of Long Read Assemblies (ILRA) is a pipeline to help in the post-assembly process (finishing) by cleaning and merging contigs, correcting homopolymer tracks and circularizing plastids. We also have the code preinstalled on a virtual Linux machine.  

## Installation
ILRA is based on several standard tools. Please find the list of software you need to install in the INSTALL.txt file.

We would recommend users without a bioinformatics setup to use our virtual machine - https://q-r.to/ILRA_VM. You will need to install Virtual Box (https://www.oracle.com/virtualization/technologies/vm/downloads/virtualbox-downloads.html), set up Ubuntu x64 and mount the disc. See the pdf VM.install.pdf for further help. The user is "bioinfo" and the password "Glasgow2020". 

Once logged in, you can run the ILRA installation with test_data by typing:<BR>
```
~/ILRA/example.sh
```
This will execute the command:
```
ILRA.sh -a $ASSEMBLY -o $OUTPUT_FOLDER_ILRA -c $CORRECTED_READS -n subset_test -r $REFERENCE -I $ILLU_READS -t $CORES -g $GFF_REF_FILE -L pb
```

Please refer to the /test_data/README file or the help ('ILRA.sh -h') for further details. The test run will take around 10 minutes ('light' mode) and around 20 minutes ('both' mode) using 4 cores.

## Comments
We used ILRA to improve many genomes. It is a continuation of the IPA (https://github.com/ThomasDOtto/IPA) project. We felt to rename the script as our updated version works with every long-read technology.
  
Please cite this preprint when using ILRA for your publications:


