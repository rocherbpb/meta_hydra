### Conda tool installation on Hydra

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
