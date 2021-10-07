#!/bin/bash
# Get the directory of the script:
EXTERNAL_SOFTWARE_DIR=$(dirname $( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ))
echo "The folder with external software is $EXTERNAL_SOFTWARE_DIR"

#### Setting permissions:
# (https://github.com/wudustan/fastq-info)
chmod 775 $EXTERNAL_SOFTWARE_DIR/fastq-info.sh
chmod 775 $EXTERNAL_SOFTWARE_DIR/fastq_info_2.sh
# (https://github.com/oXis/gtool)
chmod 775 $EXTERNAL_SOFTWARE_DIR/gtool.py
# ftp://ftp.ncbi.nlm.nih.gov/blast/demo/vecscreen
chmod 775 $EXTERNAL_SOFTWARE_DIR/vecscreen
# ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/VSlistTo1HitPerLine.awk
chmod 775 $EXTERNAL_SOFTWARE_DIR/VSlistTo1HitPerLine.sh # This wrapper allows tu use the gawk in the environment
# ILRA scripts
chmod 775 $(dirname $EXTERNAL_SOFTWARE_DIR)/ILRA.sh
cd $(dirname $EXTERNAL_SOFTWARE_DIR)/bin; chmod 775 *
# (https://github.com/satta/ABACAS2, but the copy here is updated and better optimized)
cd $EXTERNAL_SOFTWARE_DIR/ABACAS2; chmod 775 *
# (https://sourceforge.net/projects/icorn/files/, but the copy here is updated and better optimized)
cd $EXTERNAL_SOFTWARE_DIR/iCORN2; chmod 775 *


#### Download and install miniconda:
echo -e "\n\n\nI'm downloading and installing an updated copy of miniconda in the folder external_software/miniconda3 to create ILRA environment\n\n\n"
echo -e "\n\n\nIf for some reason an outdated python or conda is required in your system, please ignore the miniconda3 folder created by this script (which may have failed), go to https://repo.anaconda.com/miniconda/ and download and install manually the corresponding Linux installer\n\n\n"
cd $EXTERNAL_SOFTWARE_DIR; wget "https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh"
mkdir -p Miniconda3; sh Miniconda3-py39_4.10.3-Linux-x86_64.sh -b -f -s -p $EXTERNAL_SOFTWARE_DIR/Miniconda3; rm Miniconda3-py39_4.10.3-Linux-x86_64.sh

#### Install packages in an environment via conda:
echo -e "\n\n\nI'm downloading and installing several packages through conda...\n\n\n"
$EXTERNAL_SOFTWARE_DIR/Miniconda3/bin/conda create -n ILRA_env -y -c conda-forge mamba
source $EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/bin/activate
mamba install -y -c conda-forge pigz gawk curl openmp
mamba install -y -c bioconda blast-legacy blast samtools smalt pyfastaq bowtie2 minimap2 circlator assembly-stats fastqc bedtools pilon bwakit spades mummer4 prodigal recentrifuge
pip install --prefix $EXTERNAL_SOFTWARE_DIR/python_modules quast==5.0.2
# Small manual fix to quast
echo -e "\n\n\nSome errors when installing QUAST may be expected, I'm addressing them...\n\n\n"
cd $(dirname $(find . -name jsontemplate.py)); rm jsontemplate.py; wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/site_packages/jsontemplate/jsontemplate.py"
cd $(dirname $(find . -name misc.py | grep ra_utils)); rm misc.py; wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/ra_utils/misc.py"
cd $(dirname $(find . -name reads_analyzer.py)); rm reads_analyzer.py; "wget https://raw.githubusercontent.com/ablab/quast/master/quast_libs/reads_analyzer.py"
# Centrifuge require outdated dependencies, so new environment
mamba create -n centrifuge -y -c bioconda centrifuge


#### Setting up iCORN2:
# Getting picard:
echo -e "\n\n\nGetting picard (v2.26.2)...\n\n\n"
cd $EXTERNAL_SOFTWARE_DIR/iCORN2
wget https://github.com/broadinstitute/picard/releases/download/2.26.2/picard.jar
# Getting the outdated java version required for iCORN2:
echo -e "\n\n\nGetting JAVA v1.7 for iCORN2...\n\n\n"
cd $EXTERNAL_SOFTWARE_DIR/iCORN2
wget https://files-cdn.liferay.com/mirrors/download.oracle.com/otn-pub/java/jdk/7u80-b15/jdk-7u80-linux-x64.tar.gz; tar -xzf jdk-7u80-linux-x64.tar.gz; rm jdk-7u80-linux-x64.tar.gz
cd $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk1.7.0_80/bin; chmod 775 *

echo "The conda/mamba installation is already installing a java version that's working with picard, so not necessary anymore to get an updated java manually, but keep in mind that any related error may be to not enough updated java version..."
# Getting an updated java version required for picard:
# echo -e "\n\n\nGetting an updated JAVA version (openjdk v16.0.2) for picard...\n\n\n"
# cd $EXTERNAL_SOFTWARE_DIR/iCORN2
# et https://download.java.net/java/GA/jdk16.0.2/d4a915d82b4c4fbb9bde534da945d746/7/GPL/openjdk-16.0.2_linux-x64_bin.tar.gz; tar -xzf openjdk-16.0.2_linux-x64_bin.tar.gz; rm openjdk-16.0.2_linux-x64_bin.tar.gz
#  $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk-16.0.2/bin; chmod 775 *

echo -e "\n\n\nALL DONE\n\n\n"

