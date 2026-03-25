#!/bin/bash
#SBATCH --job-name=chr_fresco_mapping_table
#SBATCH --partition=courses
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=patel.yashm@northeastern.edu
#SBATCH --output=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.log
#SBATCH --error=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.err

# ==============================================================================
# FRESCo Mapping Table Generator - Generate ALL Mapping Tables
# ==============================================================================
#
# Purpose: Generate mapping tables for ALL gene-level FASTA files
#          This creates a reusable resource for any future FRESCo runs
#
# Usage:
#   Single file:
#     sbatch chr_fresco_mapping_table.sh single A3GALT2_CCDS60080.1.fa
#
#   All files (overwrite existing):
#     sbatch chr_fresco_mapping_table.sh all
#
#   All files (skip existing):
#     sbatch chr_fresco_mapping_table.sh all skip
#
# ==============================================================================

echo "============================================================"
echo "FRESCo Mapping Table Generation Started"
echo "Start Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "============================================================"
echo

START_TIME=$(date +%s)


# Initialize counters to avoid unary operator errors
DUPLICATE_COUNT=0
SUCCESS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

# ==============================================================================
# CONFIGURATION - EDIT THESE PATHS
# ==============================================================================

GTF_FILE="/home/patel.yashm/capstone_project/data/annotations/CCDS_hg38_one_id_per_gene.gtf"
FASTA_DIR="/home/patel.yashm/capstone_project/data/alignments/chr1/chr1_genes_fa_v2"
OUTPUT_DIR="/home/patel.yashm/capstone_project/results/fresco_output/chr1_mapping_table_v2"
SCRIPT_PATH="/home/patel.yashm/capstone_project/scripts/analysis/chr_fresco_mapping_table.py"

# ==============================================================================
# VALIDATE INPUTS
# ==============================================================================

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Validating inputs..."

if [ ! -f "$GTF_FILE" ]; then
    echo "ERROR: GTF file not found: $GTF_FILE"
    exit 1
fi

if [ ! -d "$FASTA_DIR" ]; then
    echo "ERROR: FASTA directory not found: $FASTA_DIR"
    exit 1
fi

if [ ! -f "$SCRIPT_PATH" ]; then
    echo "ERROR: Python script not found: $SCRIPT_PATH"
    exit 1
fi

echo "  ✓ GTF file: $GTF_FILE"
echo "  ✓ FASTA directory: $FASTA_DIR"
echo "  ✓ Python script: $SCRIPT_PATH"
echo

# Create output directory
mkdir -p "$OUTPUT_DIR"
echo "  ✓ Output directory: $OUTPUT_DIR"
echo

# ==============================================================================
# PARSE MODE AND OPTIONS
# ==============================================================================

MODE="${1:-all}"  # Default to 'all' if not specified
OPTION="${2}"     # Can be filename (for single mode) or 'skip' (for all mode)

SKIP_EXISTING=false

if [ "$MODE" == "single" ]; then
    if [ -z "$OPTION" ]; then
        echo "ERROR: Single mode requires a filename"
        echo "Usage: sbatch chr_fresco_mapping_table.sbatch single GENENAME_CCDSID.fa"
        exit 1
    fi
    
    SINGLE_FILE="$OPTION"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running in SINGLE FILE mode"
    echo "  Target file: $SINGLE_FILE"
    echo
    
    # Build list with single file
    FILE_LIST=("${FASTA_DIR}/${SINGLE_FILE}")
    
    if [ ! -f "${FILE_LIST[0]}" ]; then
        echo "ERROR: File not found: ${FILE_LIST[0]}"
        exit 1
    fi

elif [ "$MODE" == "all" ]; then
    
    if [ "$OPTION" == "skip" ]; then
        SKIP_EXISTING=true
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running in ALL FILES mode (SKIP EXISTING)"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running in ALL FILES mode (OVERWRITE EXISTING)"
    fi
    
    echo "  Scanning FASTA directory for gene-level alignment files..."
    echo
    
    # Get ALL FASTA files from directory
    FASTA_FILES=(${FASTA_DIR}/*.fa ${FASTA_DIR}/*.fasta)
    
    if [ ${#FASTA_FILES[@]} -eq 0 ]; then
        echo "ERROR: No FASTA files found in $FASTA_DIR"
        exit 1
    fi
    
    echo "  Found ${#FASTA_FILES[@]} FASTA files"
    
    # ==============================================================================
    # DUPLICATE DETECTION IN FASTA FILES
    # ==============================================================================
    
    echo
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking for duplicate gene/CCDS combinations..."
    
    declare -A GENE_CCDS_MAP
    DUPLICATE_COUNT=0
    
    for fasta_file in "${FASTA_FILES[@]}"; do
        # Skip if file doesn't exist (glob didn't match)
        [ -e "$fasta_file" ] || continue
        
        basename=$(basename "$fasta_file")
        
        # Remove extension
        if [[ "$basename" == *.fa ]]; then
            name_part="${basename%.fa}"
        elif [[ "$basename" == *.fasta ]]; then
            name_part="${basename%.fasta}"
        else
            continue
        fi
        
        # Check if we've seen this gene/CCDS combo before
        if [ -n "${GENE_CCDS_MAP[$name_part]}" ]; then
            echo "  ⚠️  WARNING: DUPLICATE FOUND!"
            echo "      Gene/CCDS: $name_part"
            echo "      First occurrence: ${GENE_CCDS_MAP[$name_part]}"
            echo "      Duplicate: $fasta_file"
            ((DUPLICATE_COUNT++))
        else
            GENE_CCDS_MAP[$name_part]="$fasta_file"
        fi
    done
    
    if [ $DUPLICATE_COUNT -gt 0 ]; then
        echo
        echo "  ⚠️  WARNING: Found $DUPLICATE_COUNT duplicate gene/CCDS combinations!"
        echo "      Will process each unique combination only once."
        echo
    else
        echo "  ✓ No duplicates found in FASTA files"
        echo
    fi
    
    # Build list of unique FASTA files
    FILE_LIST=()
    for name_part in "${!GENE_CCDS_MAP[@]}"; do
        FILE_LIST+=("${GENE_CCDS_MAP[$name_part]}")
    done
    
    echo "  Will process ${#FILE_LIST[@]} unique FASTA files"
    echo

else
    echo "ERROR: Invalid mode: $MODE"
    echo "Usage:"
    echo "  sbatch chr_fresco_mapping_table.sbatch single <filename>"
    echo "  sbatch chr_fresco_mapping_table.sbatch all [skip]"
    exit 1
fi

# ==============================================================================
# CHECK FOR ALREADY PROCESSED FILES
# ==============================================================================

if [ "$SKIP_EXISTING" == "true" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking for already processed files..."
    
    ALREADY_PROCESSED=0
    NEW_FILE_LIST=()
    
    for fasta_file in "${FILE_LIST[@]}"; do
        basename=$(basename "$fasta_file")
        
        # Remove extension
        if [[ "$basename" == *.fa ]]; then
            name_part="${basename%.fa}"
        elif [[ "$basename" == *.fasta ]]; then
            name_part="${basename%.fasta}"
        else
            continue
        fi
        
        output_file="${OUTPUT_DIR}/${name_part}_map.txt"
        
        if [ -f "$output_file" ]; then
            ((ALREADY_PROCESSED++))
        else
            NEW_FILE_LIST+=("$fasta_file")
        fi
    done
    
    if [ $ALREADY_PROCESSED -gt 0 ]; then
        echo "  Found $ALREADY_PROCESSED already processed files (will skip)"
        echo "  Remaining files to process: ${#NEW_FILE_LIST[@]}"
        FILE_LIST=("${NEW_FILE_LIST[@]}")
    else
        echo "  No files have been processed yet"
    fi
    echo
fi

if [ ${#FILE_LIST[@]} -eq 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] All files already processed. Nothing to do."
    echo "============================================================"
    exit 0
fi

# ==============================================================================
# PROCESS FILES
# ==============================================================================

TOTAL_FILES=${#FILE_LIST[@]}
SUCCESS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

echo "============================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing $TOTAL_FILES files"
echo "  This will create mapping tables for all gene-level FASTAs"
echo "  These tables can be reused for any future FRESCo analysis"
echo "============================================================"
echo

for i in "${!FILE_LIST[@]}"; do
    fasta_file="${FILE_LIST[$i]}"
    file_num=$((i + 1))
    
    basename=$(basename "$fasta_file")
    
    echo "----------------------------------------"
    echo "[$file_num/$TOTAL_FILES] Processing: $basename"
    echo "  Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "----------------------------------------"
    
    # Run Python script (with skip option if needed)
    if [ "$SKIP_EXISTING" == "true" ]; then
        python "$SCRIPT_PATH" "$fasta_file" "$GTF_FILE" "$OUTPUT_DIR" skip
    else
        python "$SCRIPT_PATH" "$fasta_file" "$GTF_FILE" "$OUTPUT_DIR"
    fi
    
    exit_code=$?
    
    # Check result
    if [[ "$basename" == *.fa ]]; then
        name_part="${basename%.fa}"
    elif [[ "$basename" == *.fasta ]]; then
        name_part="${basename%.fasta}"
    fi
    
    output_file="${OUTPUT_DIR}/${name_part}_map.txt"
    
    if [ $exit_code -eq 0 ]; then
        if [ -f "$output_file" ]; then
            file_size=$(stat -f%z "$output_file" 2>/dev/null || stat -c%s "$output_file" 2>/dev/null)
            row_count=$(($(wc -l < "$output_file") - 1))  # Subtract header
            echo "  Result: ✓ SUCCESS"
            echo "  Output: $output_file"
            echo "  Size: $file_size bytes"
            echo "  Rows: $row_count"
            ((SUCCESS_COUNT++))
        else
            echo "  Result: ⊘ SKIPPED (already exists)"
            ((SKIP_COUNT++))
        fi
    else
        echo "  Result: ✗ FAILED - Exit code $exit_code"
        ((FAIL_COUNT++))
    fi
    
    echo "  Completed: $(date '+%Y-%m-%d %H:%M:%S')"
    echo
done

# ==============================================================================
# SUMMARY
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
echo "  Skipped (existing): $SKIP_COUNT"
echo "  Failed:             $FAIL_COUNT"

if [ $DUPLICATE_COUNT -gt 0 ]; then
    echo "  Duplicates found:   $DUPLICATE_COUNT"
fi

echo
echo "Timing:"
echo "  Start:          $(date -d @$START_TIME '+%Y-%m-%d %H:%M:%S' 2>/dev/null || date -r $START_TIME '+%Y-%m-%d %H:%M:%S' 2>/dev/null)"
echo "  End:            $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Total runtime:  ${HOURS}h ${MINS}m ${SECS}s"
echo
echo "Output directory: $OUTPUT_DIR"
echo "Mapping tables created: $SUCCESS_COUNT"
echo
echo "These mapping tables are now ready for use with any FRESCo analysis!"
echo "============================================================"

if [ $FAIL_COUNT -gt 0 ]; then
    echo
    echo "⚠️  WARNING: Some files failed to process. Check logs above."
    exit 1
fi

exit 0
