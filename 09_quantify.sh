#!/bin/bash
#SBATCH --job-name=expression_quant2
#SBATCH --ntasks=1
#SBATCH --time 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

#Script adapted from Grace Vaziri 

# load software----------------------------------------------
module load HTSeq/0.13.5-foss-2020b-Python-3.8.6

# Directory paths
ANNOTATION=/vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/merged_transcripts2/merged.gtf
BAM_DIR=/vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/alignments2/sorted3 # Directory containing SORTED HISAT2 mapping BAM files with /csi files
OUTPUT_DIR=/vast/palmer/scratch/skelly/yaa23/RNA_seq_RASY_scratch/grace_code/expression_quant2

# Step: Quantify gene expression using HTSeq-count
echo "Quantifying gene expression with HTSeq..."
for bam_file in "$BAM_DIR"/*.bam; do
    sample_name=$(basename "$bam_file" .bam)
    htseq-count -f bam -r name -s no -t exon -i gene_id "$bam_file" \
        "$ANNOTATION" > "${OUTPUT_DIR}/${sample_name}_counts.txt"
done

echo "Expression quantification completed."

#from the htseq document: htseq-count -i gene_id -i exon_number --additional-attr gene_name --additional-atr exon_number


