env_name="16S_otdn"
conda_packages=" \
snakemake \
vsearch \
bbmap \
seqkit \
samtools \
fastqc \
salmon \
bowtie2 \
julia \
"
#d=`pwd`
#pc=$d'/miniconda3'
#if ! test -d $pc; then
#    wget -nc https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#    chmod a+x Miniconda3-latest-Linux-x86_64.sh
#    ./Miniconda3-latest-Linux-x86_64.sh  -b -p  $pc
#    rm Miniconda3-latest-Linux-x86_64.sh
#fi
#envd=$pc/envs/$env_name 
#export PATH=$pc"/bin:$PATH"
#conda init bash
#source ~/.bashrc
#if ! test -d $envd; then
#    conda create -y -n  $env_name 
#fi
echo "These packages will be installed"
echo $conda_packages
conda create -y -n  $env_name 
conda install -y  -n $env_name -c bioconda -c conda-forge $conda_packages
echo To activate workflow for analysis run command:
echo  conda activate $env_name 

