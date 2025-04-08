#!/bin/bash
#SBATCH --job-name=stringtie_assembly
#SBATCH --ntasks=1
#SBATCH --time 6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -e %x_%j.err

hostname
date

#Script adapted from Grace Vaziri 

#################################################################
# Assemble transcripts with stringtie
#################################################################

# load software
module load StringTie/2.1.4-GCCcore-10.2.0
module load parallel

INDIR=alignments2/
OUTDIR=transcripts_stringtie2/ # CHECK THIS PATH 

# accession list
ACCLIST=list.txt
ls alignments/trim_*bam | sed 's/.*\///; s/.bam//' >${ACCLIST}

# run stringtie on all samples, up to 5 in parallel
cat $ACCLIST | \
parallel -j 5 \
    "stringtie \
        -o $OUTDIR/{}.gtf \
        -p 2 \
        $INDIR/{}.bam"
