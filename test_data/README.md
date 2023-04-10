## Quick start / Minimal example
```
source ILRA/external_software/ILRA/path_to_source # To get the PATH that must be exported if you have followed the conda-based solution for installation
cd ILRA/test_data
# 'Light' mode skipping decontamination:
ILRA.sh -a assembly_Pf_test.fasta -o out_ILRA_test -c corrected_reads_Pf_test_subset.fastq.gz -n test -r PlasmoDB-47_Pfalciparum3D7_Genome_core_PMID_29862326.fasta -I Illumina_short_reads_Pf_test_subset -t 4 -g PlasmoDB-50_Pfalciparum3D7.gff -L pb -M 32g 2>&1 | tee -a out_ILRA_test_log.txt
```

The output of this test run is contained in the folder 'out_ILRA_test' and the log in the file 'out_ILRA_test_log.txt'. 


Please note that the majority of intermediate files have been removed for the sake of simplicitiy in this walkthrough, keeping only relevant logs and output files as results. The input sequences were a mock subset with 4 contigs. This run in 'light' mode with 4 threads took ~1 hour and required ~30 GB RAM. The majority of runtime (~50 minutes) and RAM usage corresponded to the last step of quality assesment of the already-corrected sequences (which could be ommited with the parameter '-q no').
This was a run in the default 'light mode'. Once several external databases are downloaded, a full run including '-m both' to activate the decontamination step based on taxonomy classification and blasting would require the installation of the databases and ~5 more minutes and ~70 GB of RAM. One contig would be identified as contamination and discarded. Please note that if you want to rerun the code above, the folder would already exist so please change the directory name in the argument '-o' or remove the folder with the output to fully rerun.

Please go through the output file 'out_ILRA_test_log.txt' and the folder 'out_ILRA_test' to get the details on the pipeline processing steps and the final output. Here, we provide the following overview:

0. Root directory
   - Excluded.contigs.fofn - Information on the contigs that have been merged or discarded during execution
   - all_output_files_log_out.txt - List of all the files contained in the directories and subfolders
   - test.ILRA.fasta - ILRA-corrected sequences
1. Filtering
   - 01.assembly.fa - sequences corrected after step 1
   - formatdb.log - Log of the indexing for BLAST
2. MegaBLAST
   - 02.assembly.fa.gz
   - 03.assembly.fa
   - mergedseq_OUT.fasta.gz
   - notcovered_OUT.fasta.gz
   - contained_OUT.fasta.gz - Sequences pre- and post-correction in this step, and the contigs identified to be contained in other, not covered by short reads, or merged
   - List.Contained.fofn
   - 02.ListContained.txt 
   - comp.self1.blast.length - List of contigs identified to be contained in others and information on the BLAST performed
   - ILRA.findoverlaps_ver3.pl_log_out.txt
   - ILRA.findoverlaps_ver3.pl_log_out_perl_warnings_errors.txt
   - ILRA.runSMALT_ver2.sh_log_out.txt
   - megablast_parallel_log_out.txt 
   - megablast_parallel_log_out_2.txt
   - results_OUT.txt 
   - log_OUT.txt - Various logs containing the output of the execution of multiple software and scripts in this step

     - Second nested list item
