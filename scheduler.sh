#! /bin/bash 
#SBATCH -J scheduler
#SBATCH -t 5-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -p all
#SBATCH --mem=2gb

cd /mnt/work1/users/pughlab/bin/swgs

snakemake --cluster-config slurm/cluster.json \
--profile slurm \
--wrapper-prefix 'file:///mnt/work1/users/pughlab/references/snakemake-wrappers/' \
--use-conda \
--use-singularity \
--jobs 5 \
--rerun-incomplete


