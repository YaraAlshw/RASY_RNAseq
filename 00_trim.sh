#!/bin/bash
#SBATCH --job-name=2trim_rawreads
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yara.alshwairikh@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
echo hostname
######################################################################
# Trimming of Reads using Trimmomatic
#################################################################

#Script adapted from Grace Vaziri 

mkdir singles2
mkdir trimmed_reads

#cd to directory with raw reads
cd /gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/raw_reads

module load Trimmomatic/0.39

for infile in *_R1_001.fastq.gz
do
base=$(basename ${infile} _R1_001.fastq.gz)
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar $Trimmomatic PE -threads 8 -phred33 ${infile} ${base}_R2_001.fastq.gz \
trim_${base}_R1_001.fastq.gz trim_singles_${base}_R1_001.fastq.gz \
trim_${base}_R2_001.fastq.gz trim_singles_${base}_R2_001.fastq.gz \
ILLUMINACLIP:/home/yaa23/project/RNA_seq_RASY/TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:45
done


# output directory
OUTDIR=/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/grace_code/trimmed_reads
mkdir -p $OUTDIR
#mv outputs to new directory

mv trim_* $OUTDIR

cd /gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/grace_code
mv trim_singles_* /gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/grace_code/singles2
