#!/bin/bash
#SBATCH --job-name=chr_fresco
#SBATCH --partition=courses
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=patel.yashm@northeastern.edu
#SBATCH --output=/home/patel.yashm/capstone_project/scripts/logs/%x_%j_%a.log
#SBATCH --error=/home/patel.yashm/capstone_project/scripts/logs/%x_%j_%a.err
#SBATCH --array=1-300%100  #300 is  total job array runs, 300 individual job runs from this single sbatch file run.

# Usage:

# sbatch chr_fresco.sh chr1

CHR="$1"

FASTA_DIR="/home/patel.yashm/capstone_project/data/alignments/${CHR}/${CHR}_genes_fa_v2"
TREE_FILE="/home/patel.yashm/capstone_project/data/tree_hg38_assembly_120mammal.tree"
FRESCO_BF="/home/patel.yashm/capstone_project/software/fresco/fresco_code/FRESCO.bf"
OUTPUT_DIR="/home/patel.yashm/capstone_project/results/fresco_output/${CHR}_v2"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Record start time
START_TIME=$(date +%s)
echo "=========================================="
echo "Job started: $(date)"
echo "Chromosome: ${CHR}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "=========================================="

# Create file list directly from FASTA directory (each task does this independently)
# Convert to array for indexing
FASTA_FILES=($(ls ${FASTA_DIR}/*.fa ${FASTA_DIR}/*.fasta 2>/dev/null | sort))

# Get total number of files
TOTAL_FILES=${#FASTA_FILES[@]}

echo "[INFO] Found ${TOTAL_FILES} FASTA files total"

# Check if this task ID is within range
if [ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL_FILES} ]; then
    echo "[WARNING] Task ID ${SLURM_ARRAY_TASK_ID} exceeds total files (${TOTAL_FILES})"
    echo "[INFO] Nothing to process for this task. Exiting."
    exit 0
fi

# Get the specific file for this task (array index is 0-based, task ID is 1-based)
FASTA_FILE="${FASTA_FILES[$((SLURM_ARRAY_TASK_ID - 1))]}"

# Check if file exists
if [ ! -f "${FASTA_FILE}" ]; then
    echo "[ERROR] FASTA file not found: ${FASTA_FILE}"
    echo "[ERROR] Task ID: ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Get basename for output naming
BASENAME=$(basename "${FASTA_FILE}" .fa)
BASENAME=$(basename "${BASENAME}" .fasta)

echo "[INFO] Processing: ${BASENAME}"
echo "[INFO] FASTA file: ${FASTA_FILE}"

# Create temporary BF file
TEMP_BF="${OUTPUT_DIR}/${BASENAME}_runFRESCO.bf"

cat > "${TEMP_BF}" <<EOF
inputRedirect = {};
inputRedirect["01"] = "${FASTA_FILE}";
inputRedirect["02"] = "${TREE_FILE}";
inputRedirect["03"] = "9";
inputRedirect["04"] = "${FRESCO_BF}";
ExecuteAFile("${FRESCO_BF}", inputRedirect);
EOF

echo "[INFO] Created temp BF file: ${TEMP_BF}"

# Run HyPhy
cd "${OUTPUT_DIR}" || exit 1

HYPHY_OUTPUT="${BASENAME}_fresco.txt"

echo "[INFO] Running HYPHYMP..."
HYPHYMP "${TEMP_BF}" > "${HYPHY_OUTPUT}" 2>&1

# Calculate processing time
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINS=$(((DURATION % 3600) / 60))
SECS=$((DURATION % 60))

if [ $? -eq 0 ]; then
    echo "[SUCCESS] Completed ${BASENAME}"
    echo "[INFO] Job finished: $(date)"
    echo "[INFO] Processing time: ${HOURS}h ${MINS}m ${SECS}s"
else
    echo "[ERROR] Failed ${BASENAME}"
    echo "[INFO] Job finished: $(date)"
    echo "[INFO] Processing time: ${HOURS}h ${MINS}m ${SECS}s"
    exit 1
fi

echo "=========================================="
