#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --ntasks=1
#SBATCH --time 6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yara.alshwairikh@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

#Script adapted from Grace Vaziri 

echo `hostname`
date

#################################################################
# Index the Genome
#################################################################

# load software
module load GCCcore/12.2.0
module load ncurses/6.3-GCCcore-12.2.0
module load libreadline/8.2-GCCcore-12.2.0
module load Tcl/8.6.12-GCCcore-12.2.0
module load SQLite/3.39.4-GCCcore-12.2.0
module load GMP/6.2.1-GCCcore-12.2.0
module load libffi/3.4.4-GCCcore-12.2.0
module load Python/3.10.8-GCCcore-12.2.0
module load HISAT2

# input/output directories
OUTDIR=hisat2_index2/

GENOME=/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/GCA_028564925.1_aRanSyl1.merge_genomic.fna

hisat2-build -p 8 $GENOME $OUTDIR/GCA_028564925.1_aRanSyl1.merge_genomic
