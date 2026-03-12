#!/bin/bash
#SBATCH --job-name=remove_hg38_gaps
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=patel.yashm@northeastern.edu
#SBATCH --output=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.log
#SBATCH --error=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.err

# ==============================================================================
# Remove hg38 Gaps from Gene Alignments - Batch Script
# ==============================================================================
#
# Purpose: Remove all gap positions from hg38 sequence and apply same
#          position removal to all 120 species in alignment
#
# Usage:
#   Single file (test mode):
#     sbatch remove_hg38_gaps.sh single A3GALT2_CCDS60080.1.fa
#
#   All files (overwrite existing):
#     sbatch remove_hg38_gaps.sh all
#
#   All files (skip already processed):
#     sbatch remove_hg38_gaps.sh all skip
#
# Output:
#   - Cleaned FASTA files in chr1_genes_fa_v2 directory
#   - Summary of files processed, skipped, and failed
#
# ==============================================================================

echo "============================================================"
echo "Remove hg38 Gaps - Batch Processing Started"
echo "Start Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "============================================================"
echo

START_TIME=$(date +%s)

SUCCESS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
NO_GAPS_COUNT=0

# ==============================================================================
# CONFIGURATION - EDIT THESE PATHS IF NEEDED
# ==============================================================================

INPUT_DIR="/home/patel.yashm/capstone_project/data/alignments/chr1/chr1_genes_fa"
OUTPUT_DIR="/home/patel.yashm/capstone_project/data/alignments/chr1/chr1_genes_fa_v2"
SCRIPT_PATH="/home/patel.yashm/capstone_project/scripts/maf_to_fa/remove_hg38_gaps.py"

# ==============================================================================
# VALIDATE INPUTS
# ==============================================================================

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Validating inputs..."

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [ ! -f "$SCRIPT_PATH" ]; then
    echo "ERROR: Python script not found: $SCRIPT_PATH"
    exit 1
fi

echo "  ✓ Input directory: $INPUT_DIR"
echo "  ✓ Output directory: $OUTPUT_DIR"
echo "  ✓ Python script: $SCRIPT_PATH"
echo

# Create output directory if doesn't exist
mkdir -p "$OUTPUT_DIR"
echo "  ✓ Output directory ready"
echo

# ==============================================================================
# PARSE MODE AND OPTIONS
# ==============================================================================

MODE="${1:-all}"
OPTION="${2}"

if [ "$MODE" == "single" ]; then
    if [ -z "$OPTION" ]; then
        echo "ERROR: Single mode requires a filename"
        echo "Usage: sbatch remove_hg38_gaps.sh single GENE_CCDSID.fa"
        exit 1
    fi
    
    SINGLE_FILE="$OPTION"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running in SINGLE FILE mode"
    echo "  Target: $SINGLE_FILE"
    echo
    
    FILE_LIST=("${INPUT_DIR}/${SINGLE_FILE}")
    
    if [ ! -f "${FILE_LIST[0]}" ]; then
        echo "ERROR: File not found: ${FILE_LIST[0]}"
        exit 1
    fi
    
    RUN_MODE="overwrite"

elif [ "$MODE" == "all" ]; then
    
    if [ "$OPTION" == "skip" ]; then
        RUN_MODE="skip"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running in ALL FILES mode (SKIP EXISTING)"
    else
        RUN_MODE="overwrite"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running in ALL FILES mode (OVERWRITE EXISTING)"
    fi
    
    echo "  Scanning input directory..."
    
    # Find all .fa files
    FILE_LIST=(${INPUT_DIR}/*.fa)
    
    if [ ${#FILE_LIST[@]} -eq 0 ]; then
        echo "ERROR: No .fa files found in $INPUT_DIR"
        exit 1
    fi
    
    echo "  Found ${#FILE_LIST[@]} FASTA files"
    echo

else
    echo "ERROR: Invalid mode: $MODE"
    echo "Usage:"
    echo "  sbatch remove_hg38_gaps.sh single <filename>"
    echo "  sbatch remove_hg38_gaps.sh all [skip]"
    exit 1
fi

# ==============================================================================
# PROCESS FILES
# ==============================================================================

TOTAL_FILES=${#FILE_LIST[@]}

# Arrays to track issues
ALL_GAPS_FILES=()
WARNING_FILES=()

echo "============================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing $TOTAL_FILES files"
echo "Mode: $RUN_MODE"
echo "============================================================"
echo

for i in "${!FILE_LIST[@]}"; do
    input_file="${FILE_LIST[$i]}"
    file_num=$((i + 1))
    
    basename=$(basename "$input_file")
    
    echo "----------------------------------------"
    echo "[$file_num/$TOTAL_FILES] $basename"
    echo "  Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "----------------------------------------"
    
    output_file="${OUTPUT_DIR}/${basename}"
    
    # Check if output exists and we're in skip mode
    if [ -f "$output_file" ] && [ "$RUN_MODE" == "skip" ]; then
        echo "  SKIPPING: Output file already exists"
        ((SKIP_COUNT++))
        echo
        continue
    fi
    
    echo "  Input:  $input_file"
    echo "  Output: $output_file"
    echo
    
    # Run Python script and capture output
    OUTPUT=$(python "$SCRIPT_PATH" "$input_file" "$output_file" 2>&1)
    exit_code=$?
    
    # Print output
    echo "$OUTPUT"
    
    # Check for special cases
    if echo "$OUTPUT" | grep -q ""No gaps in hg38""; then
        ((NO_GAPS_COUNT++))
    fi
    
    if echo "$OUTPUT" | grep -q "ALL GAPS"; then
        ALL_GAPS_INFO=$(echo "$OUTPUT" | grep "ALL GAPS")
        ALL_GAPS_FILES+=("$basename: $ALL_GAPS_INFO")
    fi
    
    if echo "$OUTPUT" | grep -q "WARNING"; then
        WARNING_FILES+=("$basename")
    fi
    
    if [ $exit_code -eq 0 ]; then
        echo "  Result: ✓ SUCCESS"
        ((SUCCESS_COUNT++))
    else
        echo "  Result: ✗ FAILED (exit code $exit_code)"
        ((FAIL_COUNT++))
    fi
    
    echo "  Completed: $(date '+%Y-%m-%d %H:%M:%S')"
    echo
done

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINS=$(((DURATION % 3600) / 60))
SECS=$((DURATION % 60))

echo "============================================================"
echo "Processing Complete!"
echo "============================================================"
echo
echo "Summary:"
echo "  Total files processed:  $TOTAL_FILES"
echo "  Successful:             $SUCCESS_COUNT"
echo "  Failed:                 $FAIL_COUNT"
echo "  Skipped:                $SKIP_COUNT"
echo "  Files with no gaps:     $NO_GAPS_COUNT"

echo
echo "Output directory: $OUTPUT_DIR"

if [ -d "$OUTPUT_DIR" ]; then
    OUTPUT_COUNT=$(ls "$OUTPUT_DIR"/*.fa 2>/dev/null | wc -l)
    echo "  Total output files:     $OUTPUT_COUNT"
fi

echo

if [ ${#ALL_GAPS_FILES[@]} -gt 0 ]; then
    echo "Files with ALL GAPS in hg38 (${#ALL_GAPS_FILES[@]} total):"
    for file_info in "${ALL_GAPS_FILES[@]}"; do
        echo "  ⚠ $file_info"
    done
    echo
fi

if [ ${#WARNING_FILES[@]} -gt 0 ]; then
    echo "Files with warnings (${#WARNING_FILES[@]} total):"
    for file in "${WARNING_FILES[@]}"; do
        echo "  ⚠ $file"
    done
    echo
fi

echo "Timing:"
echo "  Start:          $(date -d @$START_TIME '+%Y-%m-%d %H:%M:%S' 2>/dev/null || date -r $START_TIME '+%Y-%m-%d %H:%M:%S')"
echo "  End:            $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Total runtime:  ${HOURS}h ${MINS}m ${SECS}s"
echo "============================================================"

if [ $FAIL_COUNT -gt 0 ]; then
    echo
    echo "WARNING: $FAIL_COUNT files failed. Check logs above."
    exit 1
fi

exit 0
