#!/bin/bash

#SBATCH --job-name=hiller_120_align.sh            # Job name
#SBATCH --partition=short                     # Partition or queue name
#SBATCH -N 1                                    # Number of nodes
#SBATCH -c 8                                    # Number of CPU cores
#SBATCH --mem=8G                               # Total memory
#SBATCH -t 8:00:00                              # Runtime (hh:mm:ss)
#SBATCH --mail-type=END,FAIL                    # Email notifications
#SBATCH --mail-user=patel.yashm@northeastern.edu            # Your email address
#SBATCH --output=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.log   # Standard output log
#SBATCH --error=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.err    # Standard error log


# -------------------------
# Step 1: Set up directories
# -------------------------
BASE_DIR=~/capstone_project/data/alignments
mkdir -p $BASE_DIR
cd $BASE_DIR

echo "Starting MAF download on $(date)"
echo "Files will be saved in: $BASE_DIR"

# -------------------------
# Step 2: Download files (chr1–22, X, Y)
# -------------------------
URL="https://bds.mpi-cbg.de/hillerlab/120MammalAlignment/Human120way/data/maf"

for chr in {1..22} X Y; do
    FILE="chr${chr}.maf.gz"
    echo "Downloading $FILE ..."
    wget -c ${URL}/${FILE}
done

echo "Download completed on $(date)"
