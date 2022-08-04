#!/bin/bash
#SBATCH -A scope
#SBATCH -p scope-shared
#SBATCH -t 1-0:00 
#SBATCH -N 1               # 1 node
#SBATCH -c 9
#SBATCH --mem 5G
#SBATCH -J GTDB_DL
#SBATCH -o slurm/GTDBdl.%A.out 
#SBATCH -e slurm/GTDBdl.%A.out 

PYDB_DIR=/home/jmeppley/repos/py-metagenomics/databases

RELEASE=${1:-207.0}

# activate conda env
source /home/jmeppley/.bash_conda
conda activate ${PYDB_DIR}/env

J=$SLURM_CPUS_PER_TASK
snakemake -s ${PYDB_DIR}/download_gtdb.snake -j $J --rerun-incomplete -p \
	--config seqdb_root=/home/jmeppley/work/seqdbs/seqdbs \
             release=$RELEASE
