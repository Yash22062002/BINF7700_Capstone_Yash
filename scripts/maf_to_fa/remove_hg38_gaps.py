#!/usr/bin/env python3
"""
Remove gaps from hg38 sequence and apply same positions to all species
Ensures output has no gaps in hg38 while maintaining alignment
"""

import sys
import os
from datetime import datetime
from collections import OrderedDict

def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()

def read_fasta(fasta_file):
    """
    Read multi-species FASTA file into ordered dictionary
    Returns: OrderedDict {header: sequence}
    """
    sequences = OrderedDict()
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                # Accumulate sequence (handles multi-line FASTA)
                current_seq.append(line)
        
        # Don't forget last sequence
        if current_header is not None:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences

def write_fasta(sequences, output_file, width=80):
    """
    Write sequences to FASTA file with wrapped lines
    """
    with open(output_file, 'w') as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n")
            # Write sequence in lines of 'width' characters
            for i in range(0, len(seq), width):
                f.write(seq[i:i+width] + '\n')

def find_hg38_gap_positions(hg38_seq):
    """
    Find all gap positions in hg38 sequence
    Returns: list of 0-based positions where gaps occur
    """
    gap_positions = []
    for i, base in enumerate(hg38_seq):
        if base == '-':
            gap_positions.append(i)
    return gap_positions

def remove_positions(sequence, positions_to_remove):
    """
    Remove specified positions from sequence
    positions_to_remove: list of 0-based indices
    Returns: new sequence with positions removed
    """
    if not positions_to_remove:
        return sequence
    
    # Convert to list for easier manipulation
    seq_list = list(sequence)
    
    # Remove positions in reverse order to maintain indices
    for pos in sorted(positions_to_remove, reverse=True):
        if pos < len(seq_list):
            seq_list.pop(pos)
    
    return ''.join(seq_list)

def process_alignment(input_file, output_file):
    """
    Process one alignment file:
    1. Read all sequences
    2. Find gaps in hg38
    3. Remove those positions from all sequences
    4. Write output
    """
    log("="*80)
    log(f"Processing: {os.path.basename(input_file)}")
    log("="*80)
    
    # Read alignment
    log(f"Reading alignment: {input_file}")
    sequences = read_fasta(input_file)
    
    if not sequences:
        raise ValueError("No sequences found in file!")
    
    num_species = len(sequences)
    log(f"  Found {num_species} sequences")
    
    # Verify all sequences have same length
    seq_lengths = [len(seq) for seq in sequences.values()]
    unique_lengths = set(seq_lengths)
    
    if len(unique_lengths) > 1:
        log(f"  ERROR: Sequences have different lengths: {unique_lengths}")
        log(f"  Length distribution:")
        for header, seq in sequences.items():
            log(f"    {header}: {len(seq)} bp")
        raise ValueError("Sequences must all be same length in alignment!")
    
    original_length = seq_lengths[0]
    log(f"  Original alignment length: {original_length} bp")
    
    # Find hg38 sequence
    hg38_header = None
    hg38_seq = None
    
    for header in sequences.keys():
        if header.lower() == 'hg38':
            hg38_header = header
            hg38_seq = sequences[header]
            break
    
    if hg38_header is None:
        raise ValueError("hg38 sequence not found! Expected header: >hg38")
    
    log(f"  Found hg38 sequence: {hg38_header}")
    
    # Count gaps in hg38
    gap_positions = find_hg38_gap_positions(hg38_seq)
    num_gaps = len(gap_positions)
    
    log(f"  Gaps in hg38: {num_gaps} positions")
    
    if num_gaps == 0:
        log(f"  No gaps in hg38 - file is already clean!")
        log(f"  Copying to output unchanged")
        write_fasta(sequences, output_file)
        return {
            'original_length': original_length,
            'gaps_removed': 0,
            'final_length': original_length,
            'num_species': num_species,
            'status': 'no_gaps'
        }
    
    # Check if hg38 is ALL gaps
    non_gap_bases = sum(1 for base in hg38_seq if base != '-')
    if non_gap_bases == 0:
        log(f"  ERROR: hg38 sequence is ALL GAPS!")
        raise ValueError("Cannot process - hg38 has no bases, only gaps")
    
    log(f"  hg38 has {non_gap_bases} bases, {num_gaps} gaps")
    log(f"  Gap positions (0-indexed): {gap_positions[:20]}{'...' if len(gap_positions) > 20 else ''}")
    
    # Remove gap positions from ALL sequences
    log(f"  Removing {num_gaps} positions from all {num_species} sequences...")
    
    cleaned_sequences = OrderedDict()
    for header, seq in sequences.items():
        cleaned_seq = remove_positions(seq, gap_positions)
        cleaned_sequences[header] = cleaned_seq
    
    # Verify all cleaned sequences have same length
    cleaned_lengths = [len(seq) for seq in cleaned_sequences.values()]
    unique_cleaned = set(cleaned_lengths)
    
    if len(unique_cleaned) > 1:
        log(f"  ERROR: After cleaning, sequences have different lengths: {unique_cleaned}")
        raise ValueError("Cleaned sequences have mismatched lengths!")
    
    final_length = cleaned_lengths[0]
    expected_length = original_length - num_gaps
    
    if final_length != expected_length:
        log(f"  WARNING: Length mismatch!")
        log(f"    Expected: {expected_length} (original {original_length} - {num_gaps} gaps)")
        log(f"    Got: {final_length}")
    else:
        log(f"  ✓ All sequences now {final_length} bp (removed {num_gaps} positions)")
    
    # Verify hg38 has no gaps in cleaned version
    cleaned_hg38 = cleaned_sequences[hg38_header]
    remaining_gaps = cleaned_hg38.count('-')
    
    if remaining_gaps > 0:
        log(f"  ERROR: hg38 still has {remaining_gaps} gaps after cleaning!")
        raise ValueError("Gap removal failed - hg38 still contains gaps")
    
    log(f"  ✓ Verified: hg38 has 0 gaps in cleaned alignment")
    
    # Count gaps in other species
    other_species_gaps = {}
    for header, seq in cleaned_sequences.items():
        if header != hg38_header:
            gap_count = seq.count('-')
            if gap_count > 0:
                other_species_gaps[header] = gap_count
    
    if other_species_gaps:
        log(f"  Other species still have gaps (expected):")
        # Show first 5 examples
        for i, (header, gap_count) in enumerate(list(other_species_gaps.items())[:5]):
            log(f"    {header}: {gap_count} gaps")
        if len(other_species_gaps) > 5:
            log(f"    ... and {len(other_species_gaps) - 5} more species with gaps")
    
    # Write output
    log(f"Writing cleaned alignment: {output_file}")
    write_fasta(cleaned_sequences, output_file, width=80)
    
    log(f"✓ Successfully processed {os.path.basename(input_file)}")
    log(f"  Original length: {original_length} bp")
    log(f"  Gaps removed: {num_gaps} positions")
    log(f"  Final length: {final_length} bp")
    log(f"  Species: {num_species}")
    
    return {
        'original_length': original_length,
        'gaps_removed': num_gaps,
        'final_length': final_length,
        'num_species': num_species,
        'status': 'success'
    }

def main():
    start_time = datetime.now()
    
    log("="*80)
    log("Remove hg38 Gaps from Gene Alignments")
    log("="*80)
    
    if len(sys.argv) < 3:
        print("Usage: python remove_hg38_gaps.py <input_fasta> <output_fasta>")
        print()
        print("Arguments:")
        print("  input_fasta  : Input multi-species alignment FASTA file")
        print("  output_fasta : Output FASTA file (gaps removed)")
        print()
        print("Example:")
        print("  python remove_hg38_gaps.py AKIRIN1_CCDS433.1.fa output/AKIRIN1_CCDS433.1.fa")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    log(f"Input:  {input_file}")
    log(f"Output: {output_file}")
    log("")
    
    # Validate input
    if not os.path.exists(input_file):
        log(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        log(f"Created output directory: {output_dir}")
    
    try:
        result = process_alignment(input_file, output_file)
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        log("="*80)
        log(f"SUCCESS! Completed in {duration:.2f} seconds")
        log(f"Output: {output_file}")
        log("="*80)
        
    except Exception as e:
        log(f"ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
