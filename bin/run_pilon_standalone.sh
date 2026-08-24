#!/bin/bash
set -euo pipefail

ILRA_ROOT=<ILRA_GITHUB_REPO_CLONED>
EXTERNAL_SOFTWARE_DIR=$ILRA_ROOT/external_software
conda_envs_path=$EXTERNAL_SOFTWARE_DIR/ILRA

export PATH=$conda_envs_path/ILRA_env_busco/bin:$conda_envs_path/ILRA_env/bin:$ILRA_ROOT:$ILRA_ROOT/bin:$EXTERNAL_SOFTWARE_DIR:$EXTERNAL_SOFTWARE_DIR/ABACAS2:$EXTERNAL_SOFTWARE_DIR/iCORN2:$HOME/bin:$PATH
export PERL5LIB=$EXTERNAL_SOFTWARE_DIR/ABACAS2:$EXTERNAL_SOFTWARE_DIR/iCORN2
export GATK_LOCAL_JAR=$conda_envs_path/ILRA_env/share/gatk4-4.3.0.0-0/gatk-package-4.3.0.0-local.jar

TEST_DATA=$ILRA_ROOT/test_data
READ_ROOT=$TEST_DATA/Illumina_short_reads_Pf_test_subset      # expects *_1.fastq.gz / *_2.fastq.gz
REFERENCE=$TEST_DATA/assembly_Pf_test.fasta

NUM_ITER=3                    # ILRA default (matches number_iterations_icorn)
CORES=32
JAVA_MEM=64g
export _JAVA_OPTIONS="-Xms10g -Xmx$JAVA_MEM"

WORKDIR=<WHERE_YOU_WANT_YOUR_OUTPUT>
mkdir -p "$WORKDIR"
cd "$WORKDIR"

log() { echo -e "\n[$(date +%H:%M:%S)] $*"; }

current_ref=$REFERENCE

for ((i=1; i<=NUM_ITER; i++)); do
    log "=== ITERATION $i ==="
    iter_dir=$WORKDIR/iter_$i
    mkdir -p "$iter_dir"
    cd "$WORKDIR"

    log "bowtie2-build (iter $i) on $current_ref"
    \time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" \
        bowtie2-build --threads $CORES "$current_ref" genome_iter$i \
        &> bowtie_log_out.txt

    log "bowtie2 map + samtools sort (iter $i)"
    \time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" \
        bowtie2 -t -x genome_iter$i -p $CORES -X 1200 --very-sensitive -N 1 -L 31 --rdg 5,2 \
            -1 "${READ_ROOT}_1.fastq.gz" -2 "${READ_ROOT}_2.fastq.gz" 2>> bowtie_log_out.txt \
        | samtools sort -l 9 -m 2G -@ $CORES --write-index \
            -o "ill_reads$i.bam##idx##ill_reads$i.bam.bai" \
        &>> bowtie_log_out.txt

    log "pilon (iter $i) → $iter_dir/genome_pilon$i.fasta"
    \time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" \
        pilon --genome "$current_ref" --bam ill_reads$i.bam \
            --output genome_pilon$i --outdir "$iter_dir" \
            --threads $CORES --changes --vcf --tracks \
        &>> pilon_log_out.txt

    current_ref=$iter_dir/genome_pilon$i.fasta
done

log "DONE. Changes per iteration:"
wc -l $(find "$WORKDIR" -name "*.changes")

ln -sfn "iter_$NUM_ITER/genome_pilon$NUM_ITER.fasta" "$WORKDIR/final.pilon.fasta"
log "Final corrected assembly: $WORKDIR/final.pilon.fasta"
