### ILRA arguments:
```
ILRA.sh -a <Assembly> -o <Results directory> -c <Long reads corrected reads> -n <Name of results> -r <Reference genome for ABACAS2> -I <Root name of Illumina short reads> -t <Number of cores to use> -s <First sequence name to circularize> -S <Second sequence name to circularize> -i <Number of iterations for iCORN2> -f <Size threshold for discarding contigs> -R <Insert size range for Illumina short reads> -T <NCBI taxonomy id to extract> -g <GFF reference genes annotation file> -L <Long reads sequencing technology> -e <Telomeric sequence left> -E <Telomeric sequence right> -m <Execution mode> -h <Show help>
```
Parameters are not positional. If you did not provide any required parameter, the pipeline will exit or use default values if possible (check the help, the log after execution or the 'Arguments / Variables' section in the pipeline ILRA.sh).
### Please check the help page for futher details:
```
ILRA.sh -h
```
### Required software:
Many scripts (bash, perl...) within the ILRA folder are used, and you may need to manually change the interpreter in the corresponding shebang statements (first line #!) so everything works in your system. You may also need to make scripts executable (chmod) and to execute "source" or "export" so that the PATH variable and others are available for all scripts.
Many software and dependencies are also required and are going to be automatically checked. The pipeline will exit if any required software is not found in the variable PATH of your system, which you likely need to change accordingly. Check the log after execution or the 'Checking the installed software and databases' section in the pipeline ILRA.sh.
If not executing the light mode that skips decontamination (default if not specified or '-m light'), the subfolder 'databases' is going to be automatically created and must be populated manually by the user before executing ILRA. ILRA is going to provide instructions in the log after a first execution.

### Suggested execution for ILRA: (replace XXX by a number of cores to use and use the appropriate paths to the files of interest)
```
CORES=XXX
input_folder=/path/to/folder
REFERENCE=$input_folder/PlasmoDB-47_Pfalciparum3D7_Genome_core_PMID_29862326.fasta
GFF_REF_FILE=$input_folder/PlasmoDB-50_Pfalciparum3D7.gff
ILLU_READS=$input_folder/Illumina_short_reads_Pf_test_subset
ASSEMBLY=$input_folder/assembly_Pf_test.fasta
CORRECTED_READS=$input_folder/corrected_reads_Pf_test_subset.fastq.gz
OUTPUT_FOLDER_ILRA=/path/to/folder/out_ILRA_subset_test_data
mkdir -p $OUTPUT_FOLDER_ILRA

ILRA.sh -a $ASSEMBLY -o $OUTPUT_FOLDER_ILRA -c $CORRECTED_READS -n subset_test -r $REFERENCE -I $ILLU_READS -t $CORES -g $GFF_REF_FILE -L pb &> $OUTPUT_FOLDER_ILRA/output_terminal_ILRA_subset_test.txt
```