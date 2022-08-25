### Hydra tool installation

```bash 
cd ~
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-MacOSX-x86_64.sh
sh Anaconda3-2022.05-MacOSX-x86_64.sh
```
#### Afterwards, restart your terminal and run the following commands to setup Conda:
```bash
conda activate
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all -y
```

#### Additional tools are installed with Conda using the following commands:
```bash
conda install minimap2 -y
conda install nanoplot -y
conda install nanofilt -y
conda install diamond -y
conda create --name porechop_env
conda install --name porechop_env porechop python=3 gcc
```
#### Reference host genome data and Diamond NCBI nr database have already been downloaded/built in the following hydra directories:
Reference host genome files found in
```bash
/scratch/wrbu/databases/kneaddataDB/
```
Diamond nr database
```bash
/scratch/wrbu/databases/diamond/
```
Megan mapping file
```bash
/scratch/wrbu/databases/megan/
```
### Prepare files for Hydra job submission
#### create a file containing and the list of sample barcodes uses
Use ```bash nano sample_barcode.list```

#### Contents of hydra analysis job file:
```bash 
# /bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 16
#$ -q mThM.q
#$ -l mres=256G,h_data=16G,h_vmem=16G,himem
#$ -cwd
#$ -j y
#$ -N proc
#$ -o proc.log
#
# ----------------Modules------------------------- #
module load bioinformatics/samtools
module load /home/bourkeb/modulefiles/conda
source activate /home/bourkeb/anaconda3/envs/porechop_env
#
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS
#
# define your working directory
BASE_DIR=/scratch/genomics/bourkeb/redinet/mpala/test
# 
# make directories
mkdir data_porechop data_filt data_HostRem diamond
#
for filename in $(cat ${BASE_DIR}/sample_barcode.list)
do
# COPY AND PASTE THE PATH OF THE DIRECTORY CONTAING YOUR BARCODE SEQUENCE FASTQ FILES
cat /scratch/wrbu/redinet_nanopore/Mpala/2022_06_08_gDNA_leech01/2022_06_08_gDNA_leech01/20220608_1229_MC-112986_FAT07274_7b624019/fastq_pass/${filename}/*.fastq.gz > ${BASE_DIR}/data_concat/${filename}.fastq.gz
#
# remove adaptors with Porechop
porechop -i ${BASE_DIR}/data_concat/${filename}.fastq.gz -o ${BASE_DIR}/data_porechop/${filename}_PC.fastq.gz --format fastq.gz --threads $NSLOTS
#
# filter based on length < 100bp and QS of >9
gunzip -c ${BASE_DIR}/data_porechop/${filename}_PC.fastq.gz | NanoFilt -q 9 -l 100 | gzip > ${BASE_DIR}/data_filt/${filename}_filt.fastq.gz
#
# remove host genome
minimap2 -ax map-ont /scratch/wrbu/databases/kneaddataDB/tick_host/tick_host.fna.gz ${BASE_DIR}/data_filt/${filename}_filt.fastq.gz  -I 16G > ${BASE_DIR}/data_HostRem/${filename}.tick_host.sam
samtools view -f 4 ${BASE_DIR}/data_HostRem/${filename}.tick_host.sam | samtools fastq - | gzip -c - > ${BASE_DIR}/data_HostRem/${filename}_clean.fastq.gz
#
# Plot quality summaries
NanoPlot -t $NSLOTS --fastq ${BASE_DIR}/data_concat/${filename}.fastq.gz -o ${BASE_DIR}/data_concat/${filename}_plot 
NanoPlot -t $NSLOTS --fastq ${BASE_DIR}/data_filt/${filename}_filt.fastq.gz -o ${BASE_DIR}/data_filt/${filename}_plot 
NanoPlot -t $NSLOTS --fastq ${BASE_DIR}/data_HostRem/${filename}_clean.fastq.gz -o ${BASE_DIR}/data_HostRem/${filename}_plot 
#
# Diamond read classification
diamond blastx --db /scratch/wrbu/databases/diamond/nr --out ${BASE_DIR}/diamond/${filename}_DMD --outfmt 100 -q ${BASE_DIR}/data_HostRem/${filename}_clean.fastq.gz \
--threads $NSLOTS -b40 --evalue 1e-9 -F 15 --range-culling --top 10
#
# Re-formatting for Megan software
daa-meganizer --in ${BASE_DIR}/diamond/${filename}_DMD.daa --classify --mapDB /scratch/wrbu/databases/megan/megan-map-Feb2022.db --threads $NSLOTS
done
#
echo = `date` job $JOB_NAME done
```

