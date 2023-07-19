## Quick start / Minimal example
```
source ILRA/external_software/ILRA/path_to_source # To get the PATH that must be exported if you have followed the conda-based solution for installation, please make sure that no other conda environment is activated beforehand
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
   - 02.assembly.fa.gz - Sequences pre- and post-correction in this step, and the contigs identified to be contained in other, not covered by short reads, or merged
   - 03.assembly.fa
   - mergedseq_OUT.fasta.gz
   - notcovered_OUT.fasta.gz
   - contained_OUT.fasta.gz
   
   - List.Contained.fofn - List of contigs identified to be contained in others and information on the BLAST performed
   - 02.ListContained.txt
   - comp.self1.blast.length
   
   - ILRA.findoverlaps_ver3.pl_log_out.txt - Various logs containing the output of the execution of multiple software and scripts in this step
   - ILRA.findoverlaps_ver3.pl_log_out_perl_warnings_errors.txt
   - ILRA.runSMALT_ver2.sh_log_out.txt
   - megablast_parallel_log_out.txt
   - megablast_parallel_log_out_2.txt
   - results_OUT.txt
   - log_OUT.txt
3. ABACAS2
   - comp - Folder containing the output files from BLAST against reference. Meant for visualization of each contig with ACT (Artemis) and others
   - 03b.assembly.fa - Sequences post-correction in this step
   - Res.abacasBin.fna.gz - Sequences placed in the bin by ABACAS2 because they could not be ordered when comparing with the reference
   - Res.Pf3D7_01_v3.contigs.gff - Coverage plots and gff files for all the contigs in the corrected sequences
   - Res.Pf3D7_01_v3.refCoverage.plot.gz
   - Res.Pf3D7_API_v3.contigs.gff
   - Res.Pf3D7_API_v3.refCoverage.plot.gz
   - Res.Pf_M76611.contigs.gff
   - Res.Pf_M76611.refCoverage.plot.gz
   
   - abacas2.doTilingGraph.pl_parallel_log_out.txt - Various logs containing the output of the execution of multiple software and scripts in this step
   - abacas2.doTilingGraph.pl_parallel_log_out_2.txt
   - abacas_log_out.txt
   - abacas_log_out_warnings_errors.txt
   - formatdb_parallel_log_out.txt
   - formatdb_parallel_log_out_2.txt
   - megablast_parallel_log_out.txt
4. iCORN2
   - ICORN2.03b.assembly.fa.1 - Original sequences to be corrected in this step
   - ICORN2.03b.assembly.fa.2 - Sequences corrected after first iCORN2 iteration
   - ICORN2.03b.assembly.fa.3 - Sequences corrected after second iCORN2 iteration
   - ICORN2.03b.assembly.fa.4 - Sequences corrected after third iCORN2 iteration

   - log.1.out - Log for the correction of each input sequence first iteration
   - log.2.out - Log for the correction of each input sequence second iteration
   - log.3.out - Log for the correction of each input sequence third iteration

   - icorn2.mapper_log_out.txt - Various logs containing the output of the execution of multiple software and scripts in this step
   - icorn2.serial_bowtie2.sh_log_out.txt
   - processing_decompressing_log_out.txt
   - processing_decompressing_log_out_2.txt
     
     - ICORN2_1 - Output and intermediate logs for each iCORN2 iteration
     - ICORN2_2
     - ICORN2_3
5. Circlator
   - List.circular_sequences.fofn - List of contigs identified to be circular
   - Mapped.corrected.04.sam.cram - CRAM file containing the aligment of corrected long reads to the sequences to be circularized
   - 05.assembly.fa - Sequences pre- and post-circularization in this step
   - ForCirc.Ref.fasta.gz
   - ForCirc.Ref_2.fasta.gz
   - ForCirc.reads.fasta.gz
   - ForCirc.reads_2.fasta.gz

   - circlator_log_out.txt - Various logs containing the output of the execution of multiple software and scripts in this step
   - mapping_corrected_reads_log_out.txt
   - meryl_count_log_out.txt
   - meryl_print_log_out.txt
   - repetitive_k15.txt
6. Decontamination - Skipped
7. Stats
   - 07.Contigs_GC_size.txt - The log files and results of multiple scripts to evaluate the ILRA-corrected sequences, particularly the telomere regions, and the sequencing reads used
   - 07.TelomerContigs.txt
   - 07.TelomeresSummary.txt
   - 07.Telomeres_seq_1kb.txt
   - 07.assembly_stats_original_correction_ILRA.txt
   - 07.fastq_info_depth_IlluminaReads_assembly.txt

   - Illumina_short_reads_Pf_test_subset_1_fastqc.html - FastQC output of the Illumina short reads
   - Illumina_short_reads_Pf_test_subset_2_fastqc.html

   - busco_log_out.txt - Log of the execution of BUSCO
   - quast.py_log_out.txt - Log of the execution of QUAST
   - NGenomeSyn - Folder containing the logs and output of NgenomeSyn
   - asm2stats_results - Folder containing the logs and output of asm2stats
   - busco_results/test - Folder containing the logs and output of BUSCO
   - plotsr_results - Folder containing the logs and output of plotsr
   - quast_results/results_2023_04_05_10_50_46 - Folder containing the logs and output of QUAST
   - tapestry - Folder containing the logs and output of tapestry
