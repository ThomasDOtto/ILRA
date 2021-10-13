#!/bin/bash

#### This is a wrapper script to install all ILRA dependencies available through conda and finish installation of others, such as iCORN2...




#### Get the directory of the script:
EXTERNAL_SOFTWARE_DIR=$(dirname $( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ))
echo "The folder with external software is $EXTERNAL_SOFTWARE_DIR"


#### Setting permissions of executables:
# (https://github.com/oXis/gtool)
chmod 775 $EXTERNAL_SOFTWARE_DIR/gtool.py
# ftp://ftp.ncbi.nlm.nih.gov/blast/demo/vecscreen
chmod 775 $EXTERNAL_SOFTWARE_DIR/vecscreen
# ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/VSlistTo1HitPerLine.awk
chmod 775 $EXTERNAL_SOFTWARE_DIR/VSlistTo1HitPerLine.sh # This wrapper allows to use the gawk in the environment when executing VSlistTo1HitPerLine.awk
# ILRA scripts
chmod 775 $(dirname $EXTERNAL_SOFTWARE_DIR)/ILRA.sh
cd $(dirname $EXTERNAL_SOFTWARE_DIR)/bin; chmod 775 *
# (https://github.com/satta/ABACAS2, but the copy here is updated and better optimized)
cd $EXTERNAL_SOFTWARE_DIR/ABACAS2; chmod 775 *
# (https://sourceforge.net/projects/icorn/files/, but the copy here is updated and better optimized)
cd $EXTERNAL_SOFTWARE_DIR/iCORN2; chmod 775 *


#### Download and install Miniconda:
echo -e "\n\n\nI'm downloading and installing an updated copy of Miniconda3 in the folder external_software/miniconda3 to create ILRA environment\n\n\n"
echo -e "\n\n\nIf for some reason an outdated python or conda is required in your system, please ignore the miniconda3 folder created by this script (which may have failed), go to https://repo.anaconda.com/miniconda/ and download and install manually the corresponding Linux installer\n\n\n"
cd $EXTERNAL_SOFTWARE_DIR; wget "https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh"
mkdir -p Miniconda3; sh Miniconda3-py39_4.10.3-Linux-x86_64.sh -b -f -s -p $EXTERNAL_SOFTWARE_DIR/Miniconda3; rm Miniconda3-py39_4.10.3-Linux-x86_64.sh


#### Install packages in an environment via conda:
echo -e "\n\n\nI'm downloading and installing several packages through conda...\n\n\n"
$EXTERNAL_SOFTWARE_DIR/Miniconda3/bin/conda create -n ILRA_env -y -c conda-forge mamba
source $EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/bin/activate
mamba install -y -c conda-forge pigz gawk curl openmp
mamba install -y -c bioconda blast-legacy blast samtools smalt pyfastaq bowtie2 minimap2 circlator assembly-stats fastqc bedtools pilon bwakit spades mummer4 prodigal recentrifuge hmmer
# Small comment to pilon:
echo -e "\n\n\nThe JAVA options used by Pilon by default may be not appropriate for processing large genomes. For example, heap memory may not be enough. Please keep in mind that any pilon-related error may be due to this and require manual tuning. You can manually set up the environmental variable _JAVA_OPTIONS to use other -Xms or -Xmx values if needed"
# export _JAVA_OPTIONS='-Xms1g -Xmx5g'

# Centrifuge requires outdated dependencies, so new environment
echo -e "\n\n\nCentrifuge requires outdated dependencies, so new environment...\n\n\n"
mamba create -n centrifuge -y -c bioconda centrifuge

# Busco requires outdated dependencies, so new environment
echo -e "\n\n\nBUSCO requires outdated dependencies, so new environment...\n\n\n"
mamba create -n busco -y python==3.7
mamba install -n busco -y -c conda-forge -c bioconda busco==5.2.2
$EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/envs/busco/bin/pip install numpy==1.17.3
# Small manual fix to the shebang in BUSCO so it uses the appropriate python and numpy versions 
sed -i "1s%.*%\#\!$EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/envs/busco/bin/python%" $EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/envs/busco/bin/busco

# QUAST requires outdated dependencies, so new environment
echo -e "\n\n\nQUAST requires outdated dependencies, so new environment...\n\n\n"
mamba create -n quast -y -c bioconda quast==5.0.2
# Small manual fix to QUAST
echo -e "\n\n\nFew scripts from QUAST needs to be manually updated, addressing...\n\n\n"
cd $(dirname $(find $EXTERNAL_SOFTWARE_DIR -name jsontemplate.py | grep /envs/quast/)); rm jsontemplate.py; wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/site_packages/jsontemplate/jsontemplate.py"
cd $(dirname $(find $EXTERNAL_SOFTWARE_DIR -name misc.py | grep ra_utils | grep /envs/quast/)); rm misc.py; wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/ra_utils/misc.py"
cd $(dirname $(find $EXTERNAL_SOFTWARE_DIR -name reads_analyzer.py | grep /envs/quast/)); rm reads_analyzer.py; wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/reads_analyzer.py"


#### Setting up iCORN2:
# Getting FASTA Splitter:
echo -e "\n\n\nGetting fasta-splitter.pl (v0.2.6)...\n\n\n"
wget http://kirill-kryukov.com/study/tools/fasta-splitter/files/fasta-splitter-0.2.6.zip; unzip -qq fasta-splitter-0.2.6.zip; rm fasta-splitter-0.2.6.zip; chmod 775 fasta-splitter.pl
$EXTERNAL_SOFTWARE_DIR/iCORN2/fasta-splitter.pl --version
# Getting GATK:
echo -e "\n\n\nGetting GATK (v4.2.2.0)...\n\n\n"
wget https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip; unzip -qq gatk-4.2.2.0.zip; rm gatk-4.2.2.0.zip; chmod 775 gatk-4.2.2.0/gatk; ln -fs gatk-4.2.2.0/gatk gatk
$EXTERNAL_SOFTWARE_DIR/iCORN2/gatk --version
# Getting the outdated java version required for iCORN2:
echo -e "\n\n\nGetting JAVA v1.7 required by iCORN2's SNP caller...\n\n\n"
cd $EXTERNAL_SOFTWARE_DIR/iCORN2
wget https://files-cdn.liferay.com/mirrors/download.oracle.com/otn-pub/java/jdk/7u80-b15/jdk-7u80-linux-x64.tar.gz; tar -xzf jdk-7u80-linux-x64.tar.gz; rm jdk-7u80-linux-x64.tar.gz
cd $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk1.7.0_80/bin; chmod 775 *
$EXTERNAL_SOFTWARE_DIR/iCORN2/jdk1.7.0_80/bin/java -version
# Getting the outdated java version required for GATK
echo -e "\n\n\nGetting JAVA v1.8 required by GATK v4... This one must be in the PATH when executing iCORN2 and it is also enough for all the JAVA-related computation in ILRA\n\n\n"
cd $EXTERNAL_SOFTWARE_DIR/iCORN2
wget https://github.com/adoptium/temurin8-binaries/releases/download/jdk8u302-b08/OpenJDK8U-jdk_x64_linux_hotspot_8u302b08.tar.gz; tar -xzf OpenJDK8U-jdk_x64_linux_hotspot_8u302b08.tar.gz; rm OpenJDK8U-jdk_x64_linux_hotspot_8u302b08.tar.gz
cd $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk8u302-b08/bin; chmod 775 *
$EXTERNAL_SOFTWARE_DIR/iCORN2/jdk8u302-b08/bin/java -version
# Getting picard:
echo -e "\n\n\nSetting up iCORN2 dependencies..."
echo -e "Getting picard (v2.26.2)...\n\n\n"
cd $EXTERNAL_SOFTWARE_DIR/iCORN2
wget https://github.com/broadinstitute/picard/releases/download/2.26.2/picard.jar
$EXTERNAL_SOFTWARE_DIR/iCORN2/jdk8u302-b08/bin/java -jar $EXTERNAL_SOFTWARE_DIR/iCORN2/picard.jar ReorderSam --version


echo -e "\n\n\nALL DONE\n\n\n"

