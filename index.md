## Hydra tool installation

```bash 
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-MacOSX-x86_64.sh
sh Anaconda3-2022.05-MacOSX-x86_64.sh
```
#### Afterwards, restart your terminal and run the following commands to setup Anaconda:
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

## UNDER CONSTRUCTION!

### Specimen 
Specimen .  

### COI 
Raw
