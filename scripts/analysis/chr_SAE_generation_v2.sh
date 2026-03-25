#!/bin/bash
#SBATCH --job-name=chr_SAE_generation_v2
#SBATCH --partition=courses
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=patel.yashm@northeastern.edu
#SBATCH --output=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.log
#SBATCH --error=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.err

# ==============================================================================
# SAE Region BED12 Generator - Batch Script
# ==============================================================================
#
# Purpose: Extract SAE regions from FRESCo output and generate BED12 format
#          Maps alignment coordinates to genomic positions
#          Excludes stop codons from analysis
#
# Usage:
#   Single file (test mode):
#     sbatch chr_SAE_generation_v2.sh single A3GALT2_CCDS60080.1_fresco.txt
#
#   All files (overwrite existing entries):
#     sbatch chr_SAE_generation_v2.sh all
#
#   All files (skip already processed genes):
#     sbatch chr_SAE_generation_v2.sh all skip
#
# Output:
#   - BED12 file: chr1_SAE.bed (SAE regions for all genes)
#   - Genes with no SAEs listed in log file summary
#
# ==============================================================================

echo "============================================================"
echo "SAE BED12 Generation Started"
echo "Start Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "============================================================"
echo

START_TIME=$(date +%s)

DUPLICATE_COUNT=0
SUCCESS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

# ==============================================================================
# CONFIGURATION - EDIT THESE PATHS IF NEEDED
# ==============================================================================

FRESCO_DIR="/home/patel.yashm/capstone_project/results/fresco_output/chr1_p_adj_v2"
MAPPING_DIR="/home/patel.yashm/capstone_project/results/fresco_output/chr1_mapping_table_v2"
OUTPUT_BED="/home/patel.yashm/capstone_project/results/fresco_output/chr1_SAE_v2.bed"
SCRIPT_PATH="/home/patel.yashm/capstone_project/scripts/analysis/chr_SAE_generation_v2.py"

# ==============================================================================
# VALIDATE INPUTS
# ==============================================================================

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Validating inputs..."

if [ ! -d "$FRESCO_DIR" ]; then
    echo "ERROR: FRESCo directory not found: $FRESCO_DIR"
    exit 1
fi

if [ ! -d "$MAPPING_DIR" ]; then
    echo "ERROR: Mapping directory not found: $MAPPING_DIR"
    exit 1
fi

if [ ! -f "$SCRIPT_PATH" ]; then
    echo "ERROR: Python script not found: $SCRIPT_PATH"
    exit 1
fi

echo "  ✓ FRESCo directory: $FRESCO_DIR"
echo "  ✓ Mapping directory: $MAPPING_DIR"
echo "  ✓ Python script: $SCRIPT_PATH"
echo "  ✓ Output BED: $OUTPUT_BED"
echo

OUTPUT_DIR=$(dirname "$OUTPUT_BED")
mkdir -p "$OUTPUT_DIR"

# ==============================================================================
# PARSE MODE AND OPTIONS
# ==============================================================================

MODE="${1:-all}"
OPTION="${2}"

RUN_MODE="overwrite"

if [ "$MODE" == "single" ]; then
    if [ -z "$OPTION" ]; then
        echo "ERROR: Single mode requires a filename"
        echo "Usage: sbatch chr_SAE_generation.sbatch single GENE_CCDSID_fresco.txt"
        exit 1
    fi
    
    SINGLE_FILE="$OPTION"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running in SINGLE FILE mode"
    echo "  Target: $SINGLE_FILE"
    echo "  Mode: OVERWRITE"
    echo
    
    FILE_LIST=("${FRESCO_DIR}/${SINGLE_FILE}")
    
    if [ ! -f "${FILE_LIST[0]}" ]; then
        echo "ERROR: File not found: ${FILE_LIST[0]}"
        exit 1
    fi

elif [ "$MODE" == "all" ]; then
    
    if [ "$OPTION" == "skip" ]; then
        RUN_MODE="skip"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running in ALL FILES mode (SKIP EXISTING)"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running in ALL FILES mode (OVERWRITE EXISTING)"
    fi
    
    echo "  Scanning FRESCo directory..."
    
    FRESCO_FILES=(${FRESCO_DIR}/*_fresco.txt)
    
    if [ ${#FRESCO_FILES[@]} -eq 0 ]; then
        echo "ERROR: No FRESCo output files found in $FRESCO_DIR"
        exit 1
    fi
    
    echo "  Found ${#FRESCO_FILES[@]} FRESCo files"
    
    echo
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking for duplicate files..."
    
    declare -A FILE_MAP
    
    for fresco_file in "${FRESCO_FILES[@]}"; do
        basename=$(basename "$fresco_file")
        key="${basename%_fresco.txt}"
        
        if [ -n "${FILE_MAP[$key]}" ]; then
            echo "  WARNING: DUPLICATE FOUND!"
            echo "      Gene/CCDS: $key"
            echo "      First: ${FILE_MAP[$key]}"
            echo "      Duplicate: $fresco_file"
            ((DUPLICATE_COUNT++))
        else
            FILE_MAP[$key]="$fresco_file"
        fi
    done
    
    if [ $DUPLICATE_COUNT -gt 0 ]; then
        echo "  Found $DUPLICATE_COUNT duplicates (will process unique only)"
    else
        echo "  No duplicates found"
    fi
    echo
    
    FILE_LIST=()
    for key in "${!FILE_MAP[@]}"; do
        FILE_LIST+=("${FILE_MAP[$key]}")
    done
    
    echo "  Will process ${#FILE_LIST[@]} unique files"
    echo

else
    echo "ERROR: Invalid mode: $MODE"
    echo "Usage:"
    echo "  sbatch chr_SAE_generation.sbatch single <filename>"
    echo "  sbatch chr_SAE_generation.sbatch all [skip]"
    exit 1
fi

# ==============================================================================
# PROCESS FILES
# ==============================================================================

TOTAL_FILES=${#FILE_LIST[@]}

# Array to track genes with no SAEs
NO_SAE_GENES=()

echo "============================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing $TOTAL_FILES files"
echo "Output mode: $RUN_MODE"
echo "============================================================"
echo

for i in "${!FILE_LIST[@]}"; do
    fresco_file="${FILE_LIST[$i]}"
    file_num=$((i + 1))
    
    basename=$(basename "$fresco_file")
    gene_ccds="${basename%_fresco.txt}"
    
    echo "----------------------------------------"
    echo "[$file_num/$TOTAL_FILES] $basename"
    echo "  Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "----------------------------------------"
    
    mapping_file="${MAPPING_DIR}/${gene_ccds}_map.txt"
    
    if [ ! -f "$mapping_file" ]; then
        echo "  ERROR: Mapping file not found: $mapping_file"
        ((FAIL_COUNT++))
        echo
        continue
    fi
    
    echo "  FRESCo: $fresco_file"
    echo "  Mapping: $mapping_file"
    echo
    
    # Capture output to check for NO_SAE_FOUND
    OUTPUT=$(python "$SCRIPT_PATH" "$fresco_file" "$mapping_file" "$OUTPUT_BED" "$RUN_MODE" 2>&1)
    exit_code=$?
    
    # Print the output
    echo "$OUTPUT"
    
    # Check if gene had no SAEs
    if echo "$OUTPUT" | grep -q "NO_SAE_FOUND"; then
        NO_SAE_INFO=$(echo "$OUTPUT" | grep "NO_SAE_FOUND" | sed 's/.*NO_SAE_FOUND: //')
        NO_SAE_GENES+=("$NO_SAE_INFO")
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
echo "  Total files:        $TOTAL_FILES"
echo "  Successful:         $SUCCESS_COUNT"
echo "  Failed:             $FAIL_COUNT"
echo "  Skipped:            $SKIP_COUNT"

if [ $DUPLICATE_COUNT -gt 0 ]; then
    echo "  Duplicates found:   $DUPLICATE_COUNT"
fi

echo
echo "Output BED file: $OUTPUT_BED"

if [ -f "$OUTPUT_BED" ]; then
    TOTAL_ENTRIES=$(wc -l < "$OUTPUT_BED")
    UNIQUE_CCDS=$(cut -f4 "$OUTPUT_BED" | sort | uniq | wc -l)
    echo "  Total SAE entries:  $TOTAL_ENTRIES"
    echo "  Unique CCDS IDs:    $UNIQUE_CCDS"
    echo "  Genes with no SAEs: ${#NO_SAE_GENES[@]}"
fi

echo
if [ ${#NO_SAE_GENES[@]} -eq 0 ]; then
    echo "Genes without SAE regions: None (all genes had SAEs)"
else
    echo "Genes without SAE regions (${#NO_SAE_GENES[@]} total):"
    for gene_info in "${NO_SAE_GENES[@]}"; do
        echo "  - $gene_info"
    done
fi

echo
echo "Timing:"
echo "  Start:          $(date -d @$START_TIME '+%Y-%m-%d %H:%M:%S' 2>/dev/null || date -r $START_TIME '+%Y-%m-%d %H:%M:%S' 2>/dev/null)"
echo "  End:            $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Total runtime:  ${HOURS}h ${MINS}m ${SECS}s"
echo "============================================================"

if [ $FAIL_COUNT -gt 0 ]; then
    echo
    echo "WARNING: Some files failed. Check logs above."
    exit 1
fi

exit 0
