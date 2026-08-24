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

ICORN2_FRAGMENT_SIZE=500        # ILRA.sh:666
START_ITER=1
END_ITER=3                      # ILRA.sh default number_iterations_icorn
BLOCKS=4                        # -b / parallel_block_size
CORES=32
JAVA_MEM=64g
SEQ_PARTS=0                     # -P / parts_icorn2_split, ILRA default (no splitting)
LOW_MEM=no
LOW_SPA=no

WORKDIR=<WHERE_YOU_WANT_YOUR_OUTPUT>
mkdir -p "$WORKDIR"
cd "$WORKDIR"

\time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" \
    icorn2.serial_bowtie2.sh \
        "$READ_ROOT" \
        "$ICORN2_FRAGMENT_SIZE" \
        "$REFERENCE" \
        "$START_ITER" "$END_ITER" \
        "$BLOCKS" \
        "$CORES" \
        "$JAVA_MEM" \
        "$SEQ_PARTS" \
        "$LOW_MEM" "$LOW_SPA" \
    &> icorn2.serial_bowtie2.sh_log_out.txt
