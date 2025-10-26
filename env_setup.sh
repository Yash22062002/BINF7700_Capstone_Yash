#!/bin/bash

#SBATCH --job-name=environment_setup            # Job name
#SBATCH --partition=short                     # Partition or queue name
#SBATCH -N 1                                    # Number of nodes
#SBATCH -c 8                                    # Number of CPU cores
#SBATCH --mem=8G                               # Total memory
#SBATCH -t 8:00:00                              # Runtime (hh:mm:ss)
#SBATCH --mail-type=END,FAIL                    # Email notifications
#SBATCH --mail-user=patel.yashm@northeastern.edu            # Your email address
#SBATCH --output=/home/patel.yashm/capstone_project/logs/%x_%j.log   # Standard output log
#SBATCH --error=/home/patel.yashm/capstone_project/logs/%x_%j.err    # Standard error log



# Setup environment
module load miniconda3/24.11.1 #version might differ as per the system
source activate
conda activate capstone_project

module load cmake/3.30.2 #version might differ as per the system
