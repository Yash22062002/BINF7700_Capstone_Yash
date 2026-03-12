#!/bin/bash

#SBATCH --job-name=pad_all_blocks            # Job name
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
#   sbatch pad_all_blocks.sh chr1
#   sbatch pad_all_blocks.sh chr2
#   sbatch pad_all_blocks.sh chrX
###############################################

# === GENERIC USER INPUT CHROMOSOME FROM COMMAND LINE ===
CHR=$1

if [ -z "$CHR" ]; then
    echo "ERROR: No chromosome specified."
    echo "Usage: sbatch pad_blocks_chr.sbatch chr1"
    exit 1
fi

echo "[INFO] Chromosome: $CHR"

# === FIXED PATHS YOU CUSTOMIZE ONCE ===
BASE=/home/patel.yashm/capstone_project/data/alignments
SPECIES_LIST=$BASE/species_list.txt
SCRIPT=/home/patel.yashm/capstone_project/scripts/maf_to_fa/pad_all_blocks.py

# === AUTOMATIC PATHS BASED ON CHR VAR ===
INPUT_DIR=$BASE/${CHR}/${CHR}_blocks
OUTPUT_DIR=$BASE/${CHR}/${CHR}_blocks_padded
HEADER_DIR=$BASE/${CHR}
PREFIX=${CHR}_blocks_
HEADER_OUT=$HEADER_DIR/${CHR}_hg38_headers.txt

mkdir -p "$OUTPUT_DIR" "$HEADER_DIR" logs

echo "[INFO] Input dir     : $INPUT_DIR"
echo "[INFO] Output dir    : $OUTPUT_DIR"
echo "[INFO] Header output : $HEADER_OUT"

# === Run the Python script ===
module load python/  # or your conda env
# source activate capstone_project

python "$SCRIPT" \
    --species "$SPECIES_LIST" \
    --input-dir "$INPUT_DIR" \
    --prefix "$PREFIX" \
    --output-dir "$OUTPUT_DIR" \
    --header-out "$HEADER_OUT" \
    --ref-species hg38

echo "[INFO] Completed chromosome: $CHR"
