#!/bin/bash

#SBATCH --job-name=chr_pad_blocks_to_genes            # Job name
#SBATCH --partition=short                     # Partition or queue name
#SBATCH -N 1                                    # Number of nodes
#SBATCH -c 16                                    # Number of CPU cores
#SBATCH --mem=8G                               # Total memory
#SBATCH -t 4:00:00                              # Runtime (hh:mm:ss)
#SBATCH --mail-type=END,FAIL                    # Email notifications
#SBATCH --mail-user=patel.yashm@northeastern.edu            # Your email address
#SBATCH --output=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.log   # Standard output log
#SBATCH --error=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.err    # Standard error log

set -euo pipefail

###############################################
# USAGE:
#   sbatch chr_pad_blocks_to_genes.sh chr1
#   sbatch chr_pad_blocks_to_genes.sh chr2
#   sbatch chr_pad_blocks_to_genes.sh chrX
###############################################

# Usage check
if [ -z "$1" ]; then
    echo "Usage: sbatch run_chr_genes.sh <chromosome>"
    echo "Example: sbatch run_chr_genes.sh chr1"
    exit 1
fi

CHR="$1"

GTF="/home/patel.yashm/capstone_project/data/annotations/CCDS_hg38_one_id_per_gene.gtf"
BLOCKS_DIR="/home/patel.yashm/capstone_project/data/alignments/${CHR}/${CHR}_blocks_padded"
OUT_FASTA_DIR="/home/patel.yashm/capstone_project/data/alignments/${CHR}/${CHR}_genes_fa"
OUT_TSV=/home/patel.yashm/capstone_project/data/alignments/${CHR}/${CHR}_CCDSID_GENES.tsv


mkdir -p "${OUT_FASTA_DIR}" #if output gene level fasta output directory does not exist, it will create the directory



echo "[SBATCH] Chromosome       : ${CHR}"
echo "[SBATCH] GTF file         : ${GTF}"
echo "[SBATCH] Blocks directory : ${BLOCKS_DIR}"
echo "[SBATCH] FASTA output dir : ${OUT_FASTA_DIR}"
echo "[SBATCH] TSV summary path : ${OUT_TSV}"
echo "[SBATCH] Running Python..."

# Activate your environment (edit this to match Explorer setup)
# Example:
# source ~/.bashrc
# conda activate capstone_project

python chr_pad_blocks_to_genes.py \
    --gtf "${GTF}" \
    --blocks-dir "${BLOCKS_DIR}" \
    --out-fasta-dir "${OUT_FASTA_DIR}" \
    --chrom "${CHR}" \
    --out-tsv "${OUT_TSV}" \
    --hg-species "hg38"

echo "[INFO] Job finished at: $(date)"
