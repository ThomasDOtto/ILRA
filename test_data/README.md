## Quick start / Minimal example
```
source ILRA/external_software/ILRA/path_to_source # To get the PATH that must be exported if you have followed the conda-based solution for installation
cd ILRA/test_data
# 'Light' mode skipping decontamination:
ILRA.sh -a assembly_Pf_test.fasta -o out_ILRA_test -c corrected_reads_Pf_test_subset.fastq.gz -n test -r PlasmoDB-47_Pfalciparum3D7_Genome_core_PMID_29862326.fasta -I Illumina_short_reads_Pf_test_subset -t 4 -g PlasmoDB-50_Pfalciparum3D7.gff -L pb -M 32g 2>&1 | tee -a out_ILRA_test_log.txt
```

The output of this test run is contained in the folder 'out_ILRA_test' and the log in the file 'out_ILRA_test_log.txt'. 


Please note that the majority of intermediate files have been removed for the sake of simplicitiy in this walkthrough, keeping only relevant logs and output files as results. The input sequences were a mock subset with 4 contigs. This run in 'light' mode with 4 threads took ~1 hour and required ~30 GB RAM. The majority of runtime (~50 minutes) and RAM usage corresponded to the last step of quality assesment of the already-corrected sequences (which could be ommited with the parameter '-q no').
This was a run in the default 'light mode'. Once several external databases are downloaded, a full run including '-m both' to activate the decontamination step based on taxonomy classification and blasting would require the installation of the databases and ~5 more minutes and ~70 GB of RAM, and one contig would be identified as contamination and discarded

Please go through the output file 'out_ILRA_test_log.txt' and the folder 'out_ILRA_test' to get the details on the pipeline processing steps and the final output. Here, we provide the following overview:

### Pending to fully upload example and folder indexing and create the release

