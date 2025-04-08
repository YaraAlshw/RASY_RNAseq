#!/bin/bash
#SBATCH --job-name=align_hisat2
#SBATCH --ntasks=1
#SBATCH --time 2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --array=[0-19]%2

#Script adapted from Grace Vaziri 

echo hostname
echo `hostname`

#################################################################
# Align reads to genome
#################################################################
module load Python/3.10.8-GCCcore-12.2.0
module load HISAT2
module load SAMtools

INDIR=trimmed_reads
OUTDIR=alignments2 #edit to new directory on 11/26/2024

# this is an array job. 
	# one task will be spawned for each sample
	# for each task, we specify the sample as below
	# use the task ID to pull a single line, containing a single accession number from the accession list
	# then construct the file names in the call to hisat2 as below

INDEX=hisat2_index2/GCA_028564925.1_aRanSyl1.merge_genomic #edit to orrect genome on 11/26/2024

ACCLIST=list.txt

NUM=$(expr ${SLURM_ARRAY_TASK_ID} + 1)

SAMPLE=$(sed -n ${NUM}p $ACCLIST)


# run hisat2
hisat2 \
	-p 16 \
	-x $INDEX \
	-1 $INDIR/${SAMPLE}_L002_R1_001.fastq.gz \
	-2 $INDIR/${SAMPLE}_L002_R2_001.fastq.gz | \
samtools view -@ 8 -S -h -u - | \
samtools sort -@ 8 -T $SAMPLE - >$OUTDIR/$SAMPLE.bam


# index bam files
samtools index -c $OUTDIR/$SAMPLE.bam
