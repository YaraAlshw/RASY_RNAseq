#!/bin/bash
#SBATCH --job-name=entap_run
#SBATCH --ntasks=1
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=week
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#Script adapted from Grace Vaziri 

##########################################
## EnTap                                ## 
##########################################
module load DIAMOND/2.0.15-GCCcore-10.2.0

# Run EnTAP, flagging bacterial and fungal hits as contaminants, favoring hits to the Rana genus
singularity exec entap.sif EnTAP --runP \
-i /vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/transdecoder_output2/agat2.fasta.transdecoder.pep \
-d /vast/palmer/scratch/skelly/yaa23/reffseq/reffseq.dmnd \
-d /gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/grace_code/uniprot_sprot.dmnd \
-d /vast/palmer/scratch/skelly/yaa23/ntnr/nr_database2.dmnd \
--ontology_source 0  \
--threads 8 \
-c bacteria \
-c fungi \
--taxon Lithobates_sylvaticus

date
