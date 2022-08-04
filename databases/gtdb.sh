#!/bin/bash
#SBATCH -A scope
#SBATCH -p scope-shared
#SBATCH -t 2-0:00 
#SBATCH -N 1               # 1 node
#SBATCH -c 2
#SBATCH --mem 10G
#SBATCH -J GTDB_FMT
#SBATCH -o slurm/GTDBfmt.%A.out 
#SBATCH -e slurm/GTDBfmt.%A.out 

PYDB_DIR=/home/jmeppley/repos/py-metagenomics/databases

# activate conda env
source /home/jmeppley/.bash_conda
conda activate ${PYDB_DIR}/env

J=$SLURM_CPUS_PER_TASK
WAIT=30
snakemake -s ${PYDB_DIR}/format_gtdb.snake --rerun-incomplete -w $WAIT \
	--config seqdb_root=/home/jmeppley/work/seqdbs/seqdbs \
	release=95.0 fmt_threads=19 \
        -j 999 -p -k --nt --local-cores $J \
	--cluster-config cluster.yaml \
	--cluster "sbatch --mem {cluster.mem} -c {threads} -t {cluster.t} -A scope -p scope-shared"

