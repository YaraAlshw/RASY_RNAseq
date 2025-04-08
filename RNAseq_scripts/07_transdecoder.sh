#!/bin/bash
#SBATCH --job-name=transdecoder2
#SBATCH --ntasks=1
#SBATCH --time 3:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

#Script adapted from Grace Vaziri 

module load miniconda

conda activate TransDecoder_env

# Set variables for input and output
transcriptome=/vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/agat2/agat2.fasta  # StringTie-produced transcriptome file in FASTA format
output_dir=/vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/transdecoder_output2   # Directory to store TransDecoder results
orf_length=100  # Minimum ORF length in base pairs

# Predict long ORFs using TransDecoder
echo "Running TransDecoder.LongOrfs to identify ORFs..."
TransDecoder.LongOrfs -t $transcriptome --output_dir $output_dir

# Predict the best open reading frames using TransDecoder.Predict with single_best_only
echo "Running TransDecoder.Predict to translate ORFs to peptide sequences..."
TransDecoder.Predict -t $transcriptome --single_best_only --output_dir $output_dir

# Check for transcripts without peptide sequences and flag them
echo "Flagging transcripts without associated peptide sequences..."
missing_peptides_file="${output_dir}/missing_peptides.txt"

# Initialize or overwrite the missing peptides file
echo "Transcript ID" > $missing_peptides_file
for transcript in $(grep ">" $transcriptome | sed 's/>//'); do
    peptide_file="${output_dir}/${transcript}.pep"
    
    if [[ ! -f "$peptide_file" ]]; then
        # If peptide file is not found, flag the transcript
        echo "$transcript" >> $missing_peptides_file
    fi
done

echo "Annotation process complete. Check ${missing_peptides_file} for missing peptide sequences."
