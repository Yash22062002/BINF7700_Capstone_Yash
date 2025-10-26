#!/bin/bash

#SBATCH --job-name=test_env_setup            # Job name
#SBATCH --partition=short                     # Partition or queue name
#SBATCH -N 1                                    # Number of nodes
#SBATCH -c 8                                    # Number of CPU cores
#SBATCH --mem=8G                               # Total memory
#SBATCH -t 8:00:00                              # Runtime (hh:mm:ss)
#SBATCH --mail-type=END,FAIL                    # Email notifications
#SBATCH --mail-user=patel.yashm@northeastern.edu            # Your email address
#SBATCH --output=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.log   # Standard output log
#SBATCH --error=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.err    # Standard error log


echo "Testing environment setup..."
echo "Date: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"

# Test conda
conda --version

# Test HyPhy
HYPHYMP --version

# Test Python packages
python -c "import pandas; import Bio; print('Python packages OK')"

echo "Environment test complete!"
