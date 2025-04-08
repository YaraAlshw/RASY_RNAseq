#!/bin/bash
#SBATCH --job-name=merge_stringtie
#SBATCH --ntasks=1
#SBATCH --time 2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#Script adapted from Grace Vaziri 

#################################################################
# Merge stringtie assemblies
#################################################################

# load software
module load StringTie/2.1.4-GCCcore-10.2.0
module load parallel

# input/output variables
INDIR=transcripts_stringtie2/
OUTDIR=merged_transcripts2/

# stringtie GTF list
ls $INDIR/*.gtf >gtflist.txt

# run stringtie merge
stringtie --merge \
    -p 10 \
    -o $OUTDIR/merged.gtf \
    gtflist.txt
