#!/bin/bash
#SBATCH --job-name=htseq_count
#SBATCH --ntasks=1
#SBATCH --time 2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --array=[0-18]%5

#Script adapted from Grace Vaziri 

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID
date
# load software----------------------------------------------
module load HTSeq/0.13.5-foss-2020b-Python-3.8.6
#cd into input directory
cd /vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/alignments2
# output directory
OUTDIR=/vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/htseq_counts2

# bam file bash array
BAM=($(find -name "*.bam"))

BAM1=${BAM[$SLURM_ARRAY_TASK_ID]}
echo $BAM1

# get corresponding sample ID
SAM=$(basename $BAM1 | sed 's/.bam//')
echo $SAM
#execute htseq command
htseq-count -s no -r pos -f bam \
$BAM1 /vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/merged_transcripts2/merged.gtf > $OUTDIR/$SAM.counts


