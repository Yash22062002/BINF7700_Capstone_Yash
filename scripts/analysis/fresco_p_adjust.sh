#!/bin/bash
#SBATCH --job-name=fresco_p_adjust
#SBATCH --partition=courses
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH --output=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.log
#SBATCH --error=/home/patel.yashm/capstone_project/scripts/logs/%x_%j.err

###########################
#usage
#sbatch fresco_p_adjust.sh
###########################

echo "==== FRESCo Batch Processing Started ===="
start_time=$(date +%s)
echo "Start Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo

# Input/output paths
INPUT_DIR="/home/patel.yashm/capstone_project/results/fresco_output/chr1_v2"
OUTPUT_DIR="/home/patel.yashm/capstone_project/results/fresco_output/chr1_p_adj_v2"
SCRIPT_PATH="/home/patel.yashm/capstone_project/scripts/analysis/fresco_p_adjust.py"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Enable extended globbing
shopt -s nullglob

# Count files first
echo "Checking input directory: $INPUT_DIR"
file_list=("$INPUT_DIR"/*.txt)
total_files=${#file_list[@]}

echo "Total .txt files found: $total_files"

if [ $total_files -eq 0 ]; then
    echo "ERROR: No .txt files found in $INPUT_DIR"
    echo "Files in directory:"
    ls -lh "$INPUT_DIR" | head -20
    exit 1
fi

echo "First 5 files to be processed:"
for i in "${!file_list[@]}"; do
    if [ $i -ge 5 ]; then break; fi
    echo "  $((i+1)). ${file_list[$i]}"
done
echo

file_count=0
success_count=0
fail_count=0

echo "Starting processing loop..."
echo "========================================="

for file in "$INPUT_DIR"/*.txt; do
    
    # Check if file exists (redundant but safe)
    if [[ ! -f "$file" ]]; then
        echo "WARNING: Not a file, skipping: $file"
        continue
    fi
    
    filename=$(basename "$file")
    outfile="$OUTPUT_DIR/$filename"
    
    ((file_count++))
    
    echo
    echo "[$file_count/$total_files] Processing: $filename"
    echo "  Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  Input:  $file"
    echo "  Output: $outfile"
    
    # Check file size
    filesize=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null)
    echo "  Size: $filesize bytes"
    
    # Check if output already exists
    if [[ -f "$outfile" ]]; then
        echo "  Status: Output exists, will replace"
        rm -f "$outfile"
    else
        echo "  Status: New file"
    fi
    
    # Run Python script
    echo "  Running: python $SCRIPT_PATH"
    
    python "$SCRIPT_PATH" "$file" "$OUTPUT_DIR"
    exit_code=$?
    
    # Check result
    if [[ $exit_code -eq 0 ]]; then
        if [[ -f "$outfile" ]]; then
            out_size=$(stat -f%z "$outfile" 2>/dev/null || stat -c%s "$outfile" 2>/dev/null)
            echo "  Result: SUCCESS (output size: $out_size bytes)"
            ((success_count++))
        else
            echo "  Result: FAILED - No output file created!"
            ((fail_count++))
        fi
    else
        echo "  Result: FAILED - Python exit code: $exit_code"
        ((fail_count++))
    fi
    
    echo "  Completed: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "----------------------------------------"
    
done

echo
echo "========================================="
echo "==== Processing Summary ===="
echo "Total files found:     $total_files"
echo "Files processed:       $file_count"
echo "Successful:            $success_count"
echo "Failed:                $fail_count"

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo
echo "Start Time: $(date -d @$start_time '+%Y-%m-%d %H:%M:%S' 2>/dev/null || date -r $start_time '+%Y-%m-%d %H:%M:%S' 2>/dev/null)"
echo "End Time:   $(date '+%Y-%m-%d %H:%M:%S')"
echo "Total Runtime: ${runtime} seconds ($((runtime/60)) minutes)"
echo "==== Completed ===="
