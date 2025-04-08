#!/bin/bash
#SBATCH --job-name=agat
#SBATCH --ntasks=1
#SBATCH --time 2:00:00
#SBATCH --cpus-per-task=2
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
# extract and translate the coding regions:
#################################################################

# activate environment
module load miniconda
conda activate r_env

agat_sp_extract_sequences.pl -g /vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/merged_transcripts2/merged.gtf \
-f /vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/GCA_028564925.1_aRanSyl1.merge_genomic.fna \
-t exon \
--merge \
-o /vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/agat2/agat2.fasta

