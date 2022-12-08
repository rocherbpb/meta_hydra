## Whole Genome Sequencing Metagenomics
![Alt text](https://github.com/rocherbpb/meta_hydra/blob/main/metagenomics1.JPG)
## Taxonomic classification workflow

![Alt text](https://github.com/rocherbpb/meta_hydra/blob/main/classification_workflow.png)
## Prepare files for Hydra job submission
#### Log into Hydra
```bash
ssh username@hydra-login01.si.edu
```
#### Make and move into your working directory
```bash
mkdir /path_to_working_directory/
cd /path_to_working_directory/
```
#### Create a file containing the list of sample barcodes uses
Use ```nano sample_barcode.list``` to create the file and copy/paste list of barcodes used. CNTRL X to save and close.

#### Create a job file:
Use ```nano diamond.job``` to create the file and copy/paste the text below. Make the appropriate changes described. CNTRL X to save and close.
```bash 
# /bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=512G,h_data=64G,h_vmem=64G,himem
#$ -l cpu_arch='!zen'
#$ -cwd
#$ -j y
#$ -N diamond
#$ -o diamond.log
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
BASE_DIR=/home/bourkeb/meta_hydra_training/20220623_gDNA_bigleech_newtick
# 
# define directory of raw data
raw_DIR=/scratch/wrbu/redinet_nanopore/Mpala/2022_06_23_gDNA_bigleechandnewtick1
# make directories
mkdir data_concat data_porechop data_filt data_HostRem diamond
#
for filename in $(cat ${BASE_DIR}/sample_barcode.list)
do
cat ${raw_DIR}/${filename}/*.fastq.gz > ${BASE_DIR}/data_concat/${filename}.fastq.gz
#
# remove adaptors with Porechop
porechop -i ${BASE_DIR}/data_concat/${filename}.fastq.gz -o ${BASE_DIR}/data_porechop/${filename}_PC.fastq.gz --format fastq.gz --threads $NSLOTS
#
# filter based on length < 100bp and QS of >9
gunzip -c ${BASE_DIR}/data_porechop/${filename}_PC.fastq.gz | NanoFilt -q 9 -l 100 | gzip > ${BASE_DIR}/data_filt/${filename}_filt.fastq.gz
#
# remove host genome
minimap2 -ax map-ont /scratch/wrbu/databases/kneaddataDB/mosq_tick_leech_host/mosq_tick_leech_host.fna.gz ${BASE_DIR}/data_filt/${filename}_filt.fastq.gz  -I 16G > ${BASE_DIR}/data_HostRem/${filename}.host.sam
samtools view -f 4 ${BASE_DIR}/data_HostRem/${filename}.host.sam | samtools fastq - | gzip -c - > ${BASE_DIR}/data_HostRem/${filename}_clean.fastq.gz
#
# Plot quality summaries
NanoPlot -t $NSLOTS --fastq ${BASE_DIR}/data_concat/${filename}.fastq.gz -o ${BASE_DIR}/data_concat/${filename}_plot 
NanoPlot -t $NSLOTS --fastq ${BASE_DIR}/data_filt/${filename}_filt.fastq.gz -o ${BASE_DIR}/data_filt/${filename}_plot 
NanoPlot -t $NSLOTS --fastq ${BASE_DIR}/data_HostRem/${filename}_clean.fastq.gz -o ${BASE_DIR}/data_HostRem/${filename}_plot 
#
# Diamond read classification
diamond blastx --db /scratch/wrbu/databases/diamond/nr --out ${BASE_DIR}/diamond/${filename}_DMD --outfmt 100 -q ${BASE_DIR}/data_HostRem/${filename}_clean.fastq.gz \
--threads $NSLOTS -b40 --evalue 1e-9 --long-reads
#
# Re-formatting for Megan software
daa-meganizer --in ${BASE_DIR}/diamond/${filename}_DMD.daa --classify --mapDB /scratch/wrbu/databases/megan/megan-map-Feb2022.db --threads $NSLOTS --minSupport 1 --minPercentIdentity 70 --maxExpected 1.0E-9 --lcaAlgorithm longReads --lcaCoveragePercent 51 --longReads 
done
#
echo = `date` job $JOB_NAME done
```
### Submit Nanopore data to Hydra/Diamond
```bash 
qsub diamond.job
```
#### Job status can be checked using
```bash 
qstat
```
#### On completion, the analysis produces a Diamond .daa file for each barcoded sample 
### Megan analysis
Install Megan software on your desktop from https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html

Download the Diamond .daa files to your desktop using Cyberduck or Filezilla.

Open the downloaded .daa file in Megan (File-> Open...).

Default LCA parameters (Options->Change LCA Parameters..) are as follows: Max Expected = 1E-9, Min Percent Identity=90 (70 for viral exploration), Min Support Percent=0, Min Support=1, LCA Algorithm=longReads, Parse as Long Reads=ON (see figure below). Press “Apply”.

![Alt text](https://github.com/rocherbpb/meta_hydra/blob/main/LCA_params.png)

#### LCA classification of Diamond alignments
Select Rank at the level of Species. 

![Alt text](https://github.com/rocherbpb/meta_hydra/blob/main/megan_LCA.png)

To extract pathogen reads from each pathogen, select each pathogen node and then select File->Extract Reads, and “Include Summarized Reads” to save a sample pathogen file to a sample blast directory. 


#### Blast confirmation of Diamond classification
Create blast job files:

Use ```nano blast-megan.job``` to create the file and copy/paste the text below. Make the appropriate changes described. CNTRL X to save and close.

```bash 
# /bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 1
#$ -q mThM.q
#$ -l mres=300G,h_data=300G,h_vmem=300G,himem
#$ -cwd
#$ -j y
#$ -N blast-megan
#$ -o blast-megan.log
#
# ----------------Modules------------------------- #
module load bioinformatics/blast
module load /home/bourkeb/modulefiles/conda
#
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS

#
if [ -z $1 ]; then
  echo "Please give the input filename as an argument to qsub"
  exit 1
fi

FILENAME=${1}

INPUT_DIR=/scratch/genomics/bourkeb/GEIS/ticks/new_analysis/2022_10_03_cDNA_BulgariaBarcode73_96/blast

blastn \
  -task blastn \
  -db nt \
  -query ${INPUT_DIR}/${FILENAME}.fasta \
  -out ${INPUT_DIR}/${FILENAME}_blast \
  -num_threads $NSLOTS \
  -evalue 1e-9

blast2rma --in ${INPUT_DIR}/${FILENAME}_blast --reads ${INPUT_DIR}/${FILENAME}.fasta --format BlastText --out ${FILENAME}_blast.rma --minSupport 1 --minPercentIdentity 90 --maxExpected 1.0E-9 --lcaAlgorithm longReads --lcaCoveragePercent 51 --longReads
#
echo = `date` job $JOB_NAME done
```
Accept changes and press CNTRL X to save and close.

Create a blast submission file - ```nano submit_blast_jobs.sou``` and copy/paste the text below. Then CNTRL X to save and close.
```bash 
mkdir -p blast_logs
ls *.fasta -1|sed -e 's/\.fasta//' > blast_file.list
for x in $(cat blast_file.list); do 
  qsub -o blast_logs/blast-${x}.log blast-megan.job ${x}
done
```
### Submit the job file to Hydra
```bash 
source submit_blast_jobs.sou
```
#### Job status can be checked using
```bash 
qstat
```
#### View blast alignments in Megan
Download the .rma files to your computer and open them in the Megan software (File-> Open...).  

LCA parameters are again set as in Diamond analysis (Figure 1). Press “Apply”.

#### LCA classification of Blast alignments

![Alt text](https://github.com/rocherbpb/meta_hydra/blob/main/megan_blast_LCA.png)

Read classifications that confirm the initial Diamond pathogen classification are referred to as Blast confirmed pathogen reads. 

