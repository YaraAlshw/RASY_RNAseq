#!/bin/bash
#SBATCH --job-name=make_dmnd
#SBATCH --ntasks=1
#SBATCH --time 4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load DIAMOND/2.0.15-GCCcore-10.2.0

diamond makedb --in /vast/palmer/scratch/skelly/yaa23/reffseq/reffseq_database.fasta -d /vast/palmer/scratch/skelly/yaa23/reffseq/reffseq

diamond makedb --in /vast/palmer/scratch/skelly/yaa23/ntnr/nr -d /vast/palmer/scratch/skelly/yaa23/ntnr/nr_database2