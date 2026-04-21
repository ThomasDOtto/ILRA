#!/bin/bash

echo -e "\nThis is a wrapper script to install and link all ILRA dependencies available through conda (multiple environments due to dependency conflicts) and others external, such as iCORN2, jvarkit, outdated java...\n" && date


#### 0. Get the directory of the script:
EXTERNAL_SOFTWARE_DIR=$(dirname $( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ))
echo -e "\nThe folder with external software is $EXTERNAL_SOFTWARE_DIR\n"

#### 1. Download and install Miniconda if required:
if type mamba >/dev/null 2>&1; then
    echo "Mamba found. Using mamba for installation..."
    CONDA_INSTALL="mamba"
elif type conda >/dev/null 2>&1; then
    echo "Mamba not found, but conda is available. Using conda for installation..."
    CONDA_INSTALL="conda"
else
    echo >&2 "Conda/Mamba is required... and neither executable was found in the PATH."
    echo >&2 "Conda installation within ILRA folder ($EXTERNAL_SOFTWARE_DIR/miniconda3) is going to happen in 30 seconds."
    echo >&2 "Please kill the process (Ctrl+C) and correct the PATH if not necessary to install miniconda, otherwise the installation will continue in the folder $EXTERNAL_SOFTWARE_DIR/miniconda3 and create there ILRA environments\n\n\n"
    echo -e "\nIf for some reason an outdated python or conda is required in your system, please kill this process and go to https://repo.anaconda.com/miniconda/ to download and install manually the corresponding Linux installer\n"

    # 30 second countdown timer
    secs=30
    while [ $secs -gt 0 ]; do
       echo -ne "$secs\033[0K\r"
       sleep 1
       : $((secs--))
    done

    echo -e "\nProceeding with conda installation...\n"

    # Create directory and download installer without changing working directories
    mkdir -p "$EXTERNAL_SOFTWARE_DIR"
    wget -qO --tries=2 Miniconda3-latest-Linux-x86_64.sh "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    bash Miniconda3-latest-Linux-x86_64.sh -b -f -s -p "$EXTERNAL_SOFTWARE_DIR/miniconda3"
    rm Miniconda3-latest-Linux-x86_64.sh

    # Check if installation was successful
    if [[ ! -f "$EXTERNAL_SOFTWARE_DIR/miniconda3/bin/conda" ]]; then
        echo "Installation of Miniconda3 failed, please check manually... Exiting..."
        exit 1
    else
        # Prepend the newly installed conda to the PATH
        export PATH="$PWD/$EXTERNAL_SOFTWARE_DIR/miniconda3/bin:$PATH"
    fi

    # Install mamba into the base environment to speed up future installations
    conda install -y -q -c conda-forge mamba

    # Set the variable to use mamba for the final step
    CONDA_INSTALL="mamba"
fi

#### 2. Install packages via conda/mamba:
echo -e "\n\nInstalling dependencies through conda/mamba...\n\n"
echo -e "Populating the environments based on the .yml files (keep in mind, this is frozen versions of the software)..."
echo -e "The pathway to conda environments is $EXTERNAL_SOFTWARE_DIR/ILRA/...\n"
$CONDA_INSTALL clean -afyq
echo -e "\n\nInstalling ILRA_env...\n\n" && date
$CONDA_INSTALL env create -c conda-forge -c bioconda --override-channels --channel-priority flexible -yq --file $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_1.yml --prefix $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env && echo -e "\n\nDONE... ILRA_env\n\n" && date || { echo >&2 "FATAL ERROR: Failed to install ILRA_env. Exiting script."; exit 1; }
echo -e "\n\nInstalling ILRA_env_busco...\n\n" && date
$CONDA_INSTALL env create -c conda-forge -c bioconda --override-channels --channel-priority flexible -yq --file $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_2.yml --prefix $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco && echo -e "\n\nDONE... ILRA_env_busco\n\n" && date || { echo >&2 "FATAL ERROR: Failed to install ILRA_env. Exiting script."; exit 1; }
echo -e "\n\nInstalling ILRA_env_quast...\n\n" && date
$CONDA_INSTALL env create -c conda-forge -c bioconda --override-channels --channel-priority flexible -yq --file $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_3.yml --prefix $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_quast && echo -e "\n\nDONE... ILRA_env_quast\n\n" && date || { echo >&2 "FATAL ERROR: Failed to install ILRA_env. Exiting script."; exit 1; }
echo -e "\n\nInstalling ILRA_env_syri...\n\n" && date
$CONDA_INSTALL env create -c conda-forge -c bioconda --override-channels --channel-priority flexible -yq --file $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_4.yml --prefix $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_syri && echo -e "\n\nDONE... ILRA_env_syri\n\n" && date || { echo >&2 "FATAL ERROR: Failed to install ILRA_env. Exiting script."; exit 1; }
$CONDA_INSTALL clean -afyq
echo -e "\n\nInstalling bioconductor cogeqc and blobtoolkit...\n\n" && date
$EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/Rscript -e 'BiocManager::install("cogeqc")' # A package not in conda as of yet

export PATH=$EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env/bin/:$PATH
pip -q install "blobtoolkit[full]"
echo -e "\n\nDONE... cogeqc/blobtoolkit\n\n" && date


#### 3. Setting permissions of software included in the external_software_dir:
echo -e "\n\nSetting permissions in $EXTERNAL_SOFTWARE_DIR...\n\n"
# (https://github.com/oXis/gtool)
chmod 775 $EXTERNAL_SOFTWARE_DIR/gtool.py
# ftp://ftp.ncbi.nlm.nih.gov/blast/demo/vecscreen
chmod 775 $EXTERNAL_SOFTWARE_DIR/vecscreen
# ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/VSlistTo1HitPerLine.awk
chmod 775 $EXTERNAL_SOFTWARE_DIR/VSlistTo1HitPerLine.sh # This wrapper allows to use gawk in the environment when executing VSlistTo1HitPerLine.awk
# ILRA scripts
chmod 775 $(dirname $EXTERNAL_SOFTWARE_DIR)/ILRA.sh
cd $(dirname $EXTERNAL_SOFTWARE_DIR)/bin && chmod 775 *
# (https://github.com/satta/ABACAS2, but the version here is updated and better optimized)
cd $EXTERNAL_SOFTWARE_DIR/ABACAS2 && chmod 775 *
# (https://sourceforge.net/projects/icorn/files/, but the version here is updated and better optimized)
cd $EXTERNAL_SOFTWARE_DIR/iCORN2 && chmod 775 *


#### 4. Manual fixes
echo -e "\n\n\nManually fixing some dependencies and soft linking so only one conda environment has to be activated / only one directory has to be included in the PATH...\n\n\n"
## Fix canu and nucmer detection by circlator and spades 3.15.3 not working with modern version of python3
sed -i 's,Canu,canu,g' $(find $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env -name external_progs.py)
sed -i "s,'nucmer': '3.1','nucmer': '0.0',g" $(find $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env -name external_progs.py)
sed -i 's,collections.Hashable,collections.abc.Hashable,g' $(find $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env -name constructor.py | grep spades | grep pyyaml3)
## Fix busco so it can be used from the main environment (but it uses the python version in the subenvironment)
sed -i "s,/usr/bin/env python3,$EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/python3,g" $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/busco
## Fix the ILRA script to perform busco plots, so it works with the approppriate R version:
sed -i "s,/usr/bin/env Rscript,$EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/Rscript,g" $(dirname $EXTERNAL_SOFTWARE_DIR)/bin/busco_cogeqc_plots.R
## Fix and improve kraken2 2.1.2 within conda (i.e. fix download for database building and improve masking with parallel):
# https://github.com/DerrickWood/kraken2/issues/518, you just have to manually replace 'ftp' by 'https' in the line 46 of the file 'rsync_from_ncbi.pl'
echo -e "\n\n\nI'm fixing some issues with kraken2 v2.1.2... (i.e. fixing database download and improving database masking... please see the content and comments within the installation wrapper script)\n\n\n"
K2_CONDA_PATH=$(find $EXTERNAL_SOFTWARE_DIR/ILRA -name "rsync_from_ncbi.pl" | xargs dirname)
cd $K2_CONDA_PATH
sed -i 's,#^ftp:,#^https:,g' rsync_from_ncbi.pl
# database building is very slow, I incorporate this suggestion to paralelize, improving with multithreading and pipepart: https://github.com/DerrickWood/kraken2/pull/39
set +H
echo -e "#!/bin/bash\n" > temp
echo "function mask_data_chunk () { MASKER=\$1; awk -v RS=\">\" -v FS=\"\n\" -v ORS=\"\" ' { if (\$2) print \">\"\$0 } ' | \$MASKER -in - -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' }; export -f mask_data_chunk" | cat - mask_low_complexity.sh |  sed 's,#!/bin/bash,,g' >> temp
sed -i '/$MASKER -in $file -outfmt fasta/c\      parallel --pipepart -j $KRAKEN2_THREAD_CT -a $file --recstart ">" --block 250M mask_data_chunk $MASKER > $file.tmp' temp
sed -i '/$MASKER -in $target -outfmt fasta/c\    parallel --pipepart -j $KRAKEN2_THREAD_CT -a $target --recstart ">" --block 250M mask_data_chunk $MASKER > $target.tmp' temp
sed -i 's,) { ,) { \n  ,g' temp
sed -i 's,; awk,\n  awk,g' temp
sed -i "s,x/g' },x/g'\n},g" temp
sed -i 's,; export,\nexport,g' temp
rm mask_low_complexity.sh; mv temp mask_low_complexity.sh; chmod 775 mask_low_complexity.sh


#### 5. Soft linking to the main directory:
cd $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env/bin/
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/busco busco
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/metaeuk metaeuk
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/pplacer pplacer
ln -sf $(dirname $(find $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco -name sepp.py))/*.py .
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/run_sepp.py run_sepp.py
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/bowtie2* .
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/seqtk seqtk
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/time time
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_busco/bin/weave weave
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_quast/bin/quast.py quast.py
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_quast/bin/convert convert
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_syri/bin/syri syri
ln -sf $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env_syri/bin/time time
## Required to make previous versions of busco work:
# $EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/envs/busco/bin/pip install numpy==1.17.3
# Small manual fix to the shebang in BUSCO so it uses the appropriate python and numpy versions
# sed -i "1s%.*%\#\!$EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/envs/busco/bin/python%" $EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/envs/busco/bin/busco
# Make it usable from the main environ
## Required to make previous versions of quast work:
# Small manual fix to QUAST
# echo -e "\n\n\nFew scripts from QUAST needs to be manually updated, addressing...\n\n\n"
# cd $(dirname $(find $EXTERNAL_SOFTWARE_DIR -name jsontemplate.py | grep /envs/quast/)) && rm jsontemplate.py && wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/site_packages/jsontemplate/jsontemplate.py"
# cd $(dirname $(find $EXTERNAL_SOFTWARE_DIR -name misc.py | grep ra_utils | grep /envs/quast/)) && rm misc.py && wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/ra_utils/misc.py"
# cd $(dirname $(find $EXTERNAL_SOFTWARE_DIR -name reads_analyzer.py | grep /envs/quast/)) && rm reads_analyzer.py && wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/reads_analyzer.py"


#### 6. Install more software:
## Install bbtools:
echo -e "\n\n\nI'm installing bbtools...\n\n\n"
wget -q --tries=2 https://downloads.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz; tar -xzf BBMap_39.01.tar.gz && rm BBMap_39.01.tar.gz
ln -sf $(find bbmap/ -name "*.sh") .

## Install assembly-stats graphical:
echo -e "\n\n\nI'm installing assembly-stats graphical...\n\n\n"
wget -q --tries=2 https://github.com/rjchallis/assembly-stats/archive/refs/tags/17.02.tar.gz && tar -xzf 17.02.tar.gz && rm 17.02.tar.gz
ln -sf assembly-stats-17.02/pl/asm2stats.* .
chmod 775 asm2stats.*
sed -i 's,/usr/bin/perl -w,/usr/bin/env perl,g' asm2stats.*

## Install NGenomeSyn:
echo -e "\n\n\nI'm installing NGenomeSyn...\n\n\n"
wget -q --tries=2 https://github.com/hewm2008/NGenomeSyn/archive/v1.41.tar.gz && tar -xzf v1.41.tar.gz && rm v1.41.tar.gz
chmod -R 775 NGenomeSyn-1.41/bin/*
ln -sf NGenomeSyn-1.41/bin/GetTwoGenomeSyn.pl .; ln -sf NGenomeSyn-1.41/bin/NGenomeSyn .
sed -i 's,/usr/bin/perl -w,/usr/bin/env perl,g' GetTwoGenomeSyn.pl
sed -i 's,/usr/bin/perl -w,/usr/bin/env perl,g' NGenomeSyn

## Install iCORN2:
echo -e "\n\n\nI'm installing iCORN2...\n\n\n"
cd $EXTERNAL_SOFTWARE_DIR/iCORN2
## Getting SNP-o-matic:
echo -e "I'm installing iCORN2...Getting SNP-o-matic"
ln -sf $(which findknownsnps) findknownsnps
## Getting FASTA Splitter:
echo -e "I'm installing iCORN2...Getting FASTA Splitter"
ln -sf $(which fasta-splitter) fasta-splitter
## Getting GATK:
echo -e "I'm installing iCORN2...Getting GATK"
ln -sf $(which gatk) gatk
## Getting picard:
echo -e "I'm installing iCORN2...Getting Picard"
ln -sf $(find $EXTERNAL_SOFTWARE_DIR/ILRA/ILRA_env -name picard.jar | head -1) picard.jar
ln -sf $(which picard) picard
## Getting the outdated java version required for iCORN2:
echo -e "I'm installing iCORN2...Getting JAVA v1.7 required by iCORN2's SNP caller"
wget -q --tries=2 https://repo.huaweicloud.com/java/jdk/7u80-b15/jdk-7u80-linux-x64.tar.gz && tar -xzf jdk-7u80-linux-x64.tar.gz && rm jdk-7u80-linux-x64.tar.gz
if [[ ! -f $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk1.7.0_80/bin/java ]]; then
  echo "Download/Installation of java failed, please check manually the script and logs..."
  exit 1
fi
cd $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk1.7.0_80/bin && chmod 775 *
## Getting the outdated java version required for GATK
echo -e "I'm installing iCORN2...Getting JAVA v1.8 required by GATK v4. This one must be in the PATH when executing iCORN2"
cd $EXTERNAL_SOFTWARE_DIR/iCORN2
wget -q --tries=2 https://github.com/adoptium/temurin8-binaries/releases/download/jdk8u302-b08/OpenJDK8U-jdk_x64_linux_hotspot_8u302b08.tar.gz && tar -xzf OpenJDK8U-jdk_x64_linux_hotspot_8u302b08.tar.gz && rm OpenJDK8U-jdk_x64_linux_hotspot_8u302b08.tar.gz
if [[ ! -f $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk8u302-b08/bin/java ]]; then
  echo "Download/Installation of java failed, please check manually the script and logs..."
  exit 1
fi
cd $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk8u302-b08/bin && chmod 775 *

## Install jvarkit:
echo -e "\n\n\nI'm installing jvarkit v2021.10.13 to match java in ILRA_env...\n\n\n"
git clone https://github.com/lindenb/jvarkit.git && cd jvarkit && git checkout 2021.10.13
export JAVA_HOME=$EXTERNAL_SOFTWARE_DIR/iCORN2/jdk8u302-b08
export PATH=$JAVA_HOME/bin:$PATH
./gradlew -q samfixcigar
cp dist/samfixcigar.jar $EXTERNAL_SOFTWARE_DIR/samfixcigar.jar
rm -rf jvarkit

echo -e "\n\n\nALL DONE. Please make sure to handle properly PATH environments, see $EXTERNAL_SOFTWARE_DIR/ILRA/path_to_source before ILRA execution\n\n\n" && date


