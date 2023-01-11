#!/bin/bash

#### This is a wrapper script to install all ILRA dependencies available through conda, and others that are external, such as iCORN2.

#### Get the directory of the script:
EXTERNAL_SOFTWARE_DIR=$(dirname $( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ))
echo "The folder with external software is $EXTERNAL_SOFTWARE_DIR"

#### Setting permissions of executables:
echo -e "\n\nSetting permissions of executables...\n\n"
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


#### Download and install Miniconda:
echo -e "\n\nInstalling dependencies through conda...\n\n"
export conda_install=""
type conda >/dev/null 2>&1 || { echo >&2 "Conda is required to be installed... and is not being found in the PATH. Conda installation within ILRA folder is going to happen in 30 seconds. Please kill the process and correct the PATH if not necessary, otherwise let the installation continue..."; export conda_install="yes"; }
secs=$((1 * 30))
while [ $secs -gt 0 ]; do
   echo -ne "$secs\033[0K\r"
   sleep 1
   : $((secs--))
done

if [ "$conda_install" == "yes" ]; then
	echo -e "\n\n\nI'm downloading and installing an updated copy of Miniconda3 in the folder external_software/miniconda3 to create ILRA environment\n\n\n"
	echo -e "\nIf for some reason an outdated python or conda is required in your system, please kill this process and and go to https://repo.anaconda.com/miniconda/ to download and install manually the corresponding Linux installer\n"
	echo -e "\nProceeding with conda installation\n"
	cd $EXTERNAL_SOFTWARE_DIR && wget -q "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
	mkdir -p Miniconda3 && bash Miniconda3-latest-Linux-x86_64.sh -b -f -s -p $EXTERNAL_SOFTWARE_DIR/Miniconda3 && rm Miniconda3-latest-Linux-x86_64.sh
fi

if [[ ! -f $EXTERNAL_SOFTWARE_DIR/Miniconda3/bin/conda ]]; then
   echo "Installation of Miniconda3 failed or still not in the PATH, please check manually... Exiting..."
   exit 1
fi


#### Install packages via conda:
if [ "$conda_install" == "yes" ]; then
   conda_exec="$EXTERNAL_SOFTWARE_DIR/Miniconda3/bin/conda"
else
   conda_exec=$(which conda)
fi
conda_dir=$(dirname $conda_exec | sed 's,/bin,,g'); export PATH=$conda_dir/bin:$PATH
echo -e "\n\n\nI'm downloading and installing several packages through conda in $conda_dir, so ILRA and other dependencies work...\n\n\n"
echo -e "First mamba replacing conda to reduce time, and then:\npigz\ngawk\ncurl\nopenmp\nparallel\nkraken2\nkrakentools\ntaxonkit\nblast-legacy\nblast\nsamtools\nsmalt\npyfastaq\nminimap2\nwinnowmap\nassembly-stats\nfastqc\nbedtools\npilon\nbwakit\nspades\nmummer4\nprodigal\nrecentrifuge\nhmmer\ngatk4\npicard\nplotsr\nseqkit\nfasta-splitter\nsnpomatic\ngit\n\n\n"
conda create -n ILRA_env -y -q 
conda_envs_path=$(conda env list | egrep ILRA_env$ | sed 's,ILRA_env,,g' | sed 's/^ *//g' | tail -1)
echo -e "\n\n\nThe pathway to conda environments is $conda_envs_path...\n\n\n"
echo -e "\n\n\nI'm populating the environment ILRA_env...\n\n\n"
conda install -n ILRA_env -y -q -c conda-forge mamba
source $conda_envs_path/ILRA_env/bin/activate
mamba install -y -q -c conda-forge pigz gawk curl openmp parallel
mamba install -y -q -c bioconda -c conda-forge -c anaconda kraken2==2.1.2 krakentools taxonkit blast-legacy blast samtools==1.16.1 smalt pyfastaq minimap2 winnowmap assembly-stats fastqc bedtools pilon bwakit spades mummer4 prodigal recentrifuge hmmer gatk4 picard plotsr seqkit fasta-splitter snpomatic git
rm -rf $(dirname $conda_envs_path)/pkgs/*;rm -rf $conda_envs_path/ILRA_env/pkgs/*  # Conda creates a huge amount of intermediate files when installing packages, and these can be removed afterwards

#### Small comment to pilon:
echo -e "\n\n\nThe JAVA options used by Pilon by default may be not appropriate for processing large genomes. For example, heap memory may not be enough. Please keep in mind that any pilon-related error may be due to this and require manual tuning. You can manually set up the environmental variable _JAVA_OPTIONS to use other -Xms or -Xmx values if needed. For example, to set minimum 5GB and maximum 250GB, use export _JAVA_OPTIONS='-Xms5g -Xmx250g'"
# export _JAVA_OPTIONS='-Xms1g -Xmx5g'

#### Install circlator
echo -e "\n\n\nI'm installing circlator...\n\n\n"
$conda_envs_path/ILRA_env/bin/pip3 -q install circlator
# Fix canu and nucmer detection by circlator:
sed -i 's,Canu,canu,g' $(find $conda_envs_path/ILRA_env -name external_progs.py)
sed -i "s,'nucmer': '3.1','nucmer': '0.0',g" $(find $conda_envs_path/ILRA_env -name external_progs.py)
# Fix spades 3.15.3 not working with python 3.10
sed -i 's,collections.Hashable,collections.abc.Hashable,g' $(find $conda_envs_path/ILRA_env -name constructor.py | grep spades | grep pyyaml3)

#### Install Bowtie2:
cd $conda_envs_path/ILRA_env/bin
echo -e "\n\n\nI'm installing Bowtie2...\n\n\n"
wget -q https://github.com/BenLangmead/bowtie2/releases/download/v2.5.0/bowtie2-2.5.0-linux-x86_64.zip
unzip -qq bowtie2-2.5.0-linux-x86_64.zip && rm bowtie2-2.5.0-linux-x86_64.zip && ln -sf bowtie2-2.5.0-linux-x86_64/bowtie2 .

#### Install jvarkit:
echo -e "\n\n\nI'm installing jvarkit...\n\n\n"
git clone -q "https://github.com/lindenb/jvarkit.git"; JAVA_HOME=$conda_envs_path/ILRA_env
cd jvarkit && ./gradlew -q samfixcigar && mv dist/samfixcigar.jar ../samfixcigar.jar

#### Install BUSCO (it requires outdated dependencies, so new environment):
echo -e "\n\n\nI'm installing BUSCO...\n\n\n"
mamba create -n busco -y -q -c conda-forge -c bioconda -c anaconda busco=5.4.4
rm -rf $(dirname $conda_envs_path)/pkgs/*;rm -rf $conda_envs_path/ILRA_env/pkgs/*  # Conda creates a huge amount of intermediate files when installing packages, and these can be removed afterwards
# Required to make previous versions of busco work:
# $EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/envs/busco/bin/pip install numpy==1.17.3
# Small manual fix to the shebang in BUSCO so it uses the appropriate python and numpy versions
# sed -i "1s%.*%\#\!$EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/envs/busco/bin/python%" $EXTERNAL_SOFTWARE_DIR/Miniconda3/envs/ILRA_env/envs/busco/bin/busco
# Make it usable from the main environ
cd $conda_envs_path/ILRA_env/bin
ln -sf $conda_envs_path/ILRA_env/envs/busco/bin/busco busco; ln -sf $conda_envs_path/ILRA_env/envs/busco/bin/metaeuk metaeuk; ln -sf $conda_envs_path/ILRA_env/envs/busco/bin/pplacer pplacer; ln -sf $(dirname $(find $conda_envs_path/ILRA_env/envs/busco -name sepp.py))/*.py .; ln -sf $conda_envs_path/ILRA_env/envs/busco/bin/run_sepp.py run_sepp.py
sed -i "s,/usr/bin/env python3,$conda_envs_path/ILRA_env/envs/busco/bin/python3,g" busco

#### Install QUAST (it requires outdated dependencies, so new environment):
echo -e "\n\n\nI'm installing QUAST...\n\n\n"
mamba create -n quast -y -q -c bioconda -c conda-forge -c anaconda quast==5.2.0
rm -rf $(dirname $conda_envs_path)/pkgs/*;rm -rf $conda_envs_path/ILRA_env/pkgs/*  # Conda creates a huge amount of intermediate files when installing packages, and these can be removed afterwards
# Required to make previous versions of quast work:
# Small manual fix to QUAST
# echo -e "\n\n\nFew scripts from QUAST needs to be manually updated, addressing...\n\n\n"
# cd $(dirname $(find $EXTERNAL_SOFTWARE_DIR -name jsontemplate.py | grep /envs/quast/)) && rm jsontemplate.py && wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/site_packages/jsontemplate/jsontemplate.py"
# cd $(dirname $(find $EXTERNAL_SOFTWARE_DIR -name misc.py | grep ra_utils | grep /envs/quast/)) && rm misc.py && wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/ra_utils/misc.py"
# cd $(dirname $(find $EXTERNAL_SOFTWARE_DIR -name reads_analyzer.py | grep /envs/quast/)) && rm reads_analyzer.py && wget "https://raw.githubusercontent.com/ablab/quast/master/quast_libs/reads_analyzer.py"
# Make it usable from the main environ
cd $conda_envs_path/ILRA_env/bin; ln -sf $conda_envs_path/ILRA_env/envs/quast/bin/quast.py quast.py

#### Install syri
echo -e "\n\n\nI'm installing syri...\n\n\n"
mamba create -n syri_env -y -q -c bioconda -c conda-forge -c anaconda syri
rm -rf $(dirname $conda_envs_path)/pkgs/*;rm -rf $conda_envs_path/ILRA_env/pkgs/*  # Conda creates a huge amount of intermediate files when installing packages, and these can be removed afterwards
# Make it usable from the main environ
cd $conda_envs_path/ILRA_env/bin; ln -sf $conda_envs_path/ILRA_env/envs/syri_env/bin/syri syri

#### Install assembly-stats graphical:
echo -e "\n\n\nI'm installing assembly-stats graphical...\n\n\n"
cd $conda_envs_path/ILRA_env/bin
wget -q https://github.com/rjchallis/assembly-stats/archive/refs/tags/17.02.tar.gz && tar -xzf 17.02.tar.gz && rm 17.02.tar.gz
ln -sf assembly-stats-17.02/pl/asm2stats.* .
chmod 775 asm2stats.*
sed -i 's,/usr/bin/perl -w,/usr/bin/env perl,g' asm2stats.*

#### Fix kraken2 2.1.2 within conda (e.g. for database building):
# https://github.com/DerrickWood/kraken2/issues/518, you just have to manually replace 'ftp' by 'https' in the line 46 of the file 'rsync_from_ncbi.pl'
echo -e "\n\n\nI'm fixing some issues with kraken2...(see the content and comments within the script)\n\n\n"
sed -i 's,#^ftp:,#^https:,g' $conda_envs_path/ILRA_env/libexec/rsync_from_ncbi.pl
# database building is very slow, I incorporate this suggestion to paralelize, improving with multithreading and pipepart: https://github.com/DerrickWood/kraken2/pull/39
cd $conda_envs_path/ILRA_env/libexec/
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


#### Install iCORN2:
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
# $EXTERNAL_SOFTWARE_DIR/iCORN2/gatk --version
## Getting picard:
echo -e "I'm installing iCORN2...Getting Picard"
ln -sf $(find $conda_envs_path/ILRA_env -name picard.jar | head -1) picard.jar
ln -sf $(which picard) picard
# picard MarkDuplicates --version
## Getting the outdated java version required for iCORN2:
echo -e "I'm installing iCORN2...Getting JAVA v1.7 required by iCORN2's SNP caller"
wget -q https://files-cdn.liferay.com/mirrors/download.oracle.com/otn-pub/java/jdk/7u80-b15/jdk-7u80-linux-x64.tar.gz && tar -xzf jdk-7u80-linux-x64.tar.gz && rm jdk-7u80-linux-x64.tar.gz
cd $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk1.7.0_80/bin && chmod 775 *
$EXTERNAL_SOFTWARE_DIR/iCORN2/jdk1.7.0_80/bin/java -version
if [[ ! -f $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk1.7.0_80/bin/java ]]; then
  echo "Download/Installation of java failed, please check manually the script and logs..."
  exit 1
fi
## Getting the outdated java version required for GATK
echo -e "I'm installing iCORN2...Getting JAVA v1.8 required by GATK v4. This one must be in the PATH when executing iCORN2"
cd $EXTERNAL_SOFTWARE_DIR/iCORN2
wget -q https://github.com/adoptium/temurin8-binaries/releases/download/jdk8u302-b08/OpenJDK8U-jdk_x64_linux_hotspot_8u302b08.tar.gz && tar -xzf OpenJDK8U-jdk_x64_linux_hotspot_8u302b08.tar.gz && rm OpenJDK8U-jdk_x64_linux_hotspot_8u302b08.tar.gz
cd $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk8u302-b08/bin && chmod 775 *
$EXTERNAL_SOFTWARE_DIR/iCORN2/jdk8u302-b08/bin/java -version
if [[ ! -f $EXTERNAL_SOFTWARE_DIR/iCORN2/jdk8u302-b08/bin/java ]]; then
  echo "Download/Installation of java failed, please check manually the script and logs..."
  exit 1
fi


echo -e "\n\n\nALL DONE\n\n\n"


