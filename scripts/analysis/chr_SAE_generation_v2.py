#!/usr/bin/env python3
"""
Generate BED12 file of SAE regions from FRESCo output
Maps significant SAE windows to genomic coordinates
VERSION 2: Uses codon-overlap merging instead of consecutive-only
"""

import sys
import os
import pandas as pd
from datetime import datetime
from collections import defaultdict

def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()

def load_fresco_output(fresco_file):
    """Load FRESCo output with p_adj and classification"""
    log(f"Loading FRESCo output: {fresco_file}")
    df = pd.read_csv(fresco_file, sep='\t', comment=None)
    df.columns = df.columns.str.replace('^#', '', regex=True).str.strip()
    log(f"  Total windows: {len(df)}")
    if 'classification' not in df.columns:
        log(f"  ERROR: Available columns: {list(df.columns)}")
        raise ValueError("FRESCo file missing 'classification' column")
    if 'p_adj_bonf' not in df.columns:
        log(f"  ERROR: Available columns: {list(df.columns)}")
        raise ValueError("FRESCo file missing 'p_adj_bonf' column")
    sae_count = (df['classification'] == 'SAE').sum()
    sce_count = (df['classification'] == 'SCE').sum()
    normal_count = (df['classification'] == 'normal').sum()
    log(f"  SAE windows: {sae_count}")
    log(f"  SCE windows: {sce_count}")
    log(f"  Normal windows: {normal_count}")
    return df

def load_mapping_table(map_file):
    """Load mapping table for this gene"""
    log(f"Loading mapping table: {map_file}")
    df = pd.read_csv(map_file, sep='\t')
    log(f"  Total positions: {len(df)}")
    required = ['Chromosome', 'Alignment_Codon', 'Alignment_Index', 'Alignment_Base', 'CCDS_Index', 'CCDS_Base_Location', 'CCDS_Coordinates', 'CCDS_ID', 'Gene_Name', 'Strand']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Mapping table missing column: {col}")
    chrom = df['Chromosome'].iloc[0]
    strand = df['Strand'].iloc[0]
    ccds_id = df['CCDS_ID'].iloc[0]
    gene_name = df['Gene_Name'].iloc[0]
    log(f"  Gene: {gene_name}")
    log(f"  CCDS: {ccds_id}")
    log(f"  Chromosome: {chrom}")
    log(f"  Strand: {strand}")
    return df, chrom, strand, ccds_id, gene_name

def find_sae_windows(fresco_df):
    """Extract SAE windows from FRESCo output"""
    sae_df = fresco_df[fresco_df['classification'] == 'SAE'].copy()
    if len(sae_df) == 0:
        log("  No SAE windows found")
        return None
    sae_windows = sae_df.iloc[:, 0].tolist()
    sae_windows.sort()
    log(f"  Found {len(sae_windows)} SAE windows")
    log(f"  Window range: {min(sae_windows)}-{max(sae_windows)}")
    return sae_windows, sae_df

def merge_overlapping_windows(sae_windows, window_size=9):
    """
    Merge SAE windows based on overlapping codon coverage
    
    Window N covers codons N to N+window_size-1 (default: 9-codon window)
    Two windows merge if their codon ranges overlap
    
    Example with window_size=9:
      Window 0: covers codons 0-8
      Window 3: covers codons 3-11
      Overlap at codons 3-8 → MERGE into one region (windows 0-3)
    
    This is more biologically meaningful than merging only consecutive windows,
    as it accounts for the sliding window nature of the analysis.
    """
    if not sae_windows:
        return []
    
    # Sort windows
    sae_windows = sorted(sae_windows)
    
    merged_regions = []
    current_start = sae_windows[0]
    current_end = sae_windows[0]
    
    for i in range(1, len(sae_windows)):
        next_window = sae_windows[i]
        
        # Calculate codon coverage
        # Current region covers codons: current_start to (current_end + window_size - 1)
        current_codon_end = current_end + window_size - 1
        
        # Next window covers codons: next_window to (next_window + window_size - 1)
        next_codon_start = next_window
        
        # Check if codon ranges overlap
        if next_codon_start <= current_codon_end:
            # Overlapping codon coverage! Extend the current region
            current_end = next_window
        else:
            # No overlap - save current region and start new one
            merged_regions.append((current_start, current_end))
            current_start = next_window
            current_end = next_window
    
    # Don't forget the last region
    merged_regions.append((current_start, current_end))
    
    # Log the merging results
    log(f"  Merged {len(sae_windows)} windows into {len(merged_regions)} SAE regions (codon-overlap method):")
    for i, (start, end) in enumerate(merged_regions, 1):
        num_windows = end - start + 1
        codon_start = start
        codon_end = end + window_size - 1
        log(f"    Region {i}: windows {start}-{end} ({num_windows} windows) → codons {codon_start}-{codon_end}")
    
    return merged_regions

# [REST OF THE FUNCTIONS REMAIN THE SAME - verify_stop_codon, map_sae_region_to_genomic, etc.]

def verify_stop_codon(bases):
    """Check if 3 bases form a stop codon"""
    codon = ''.join(bases).upper()
    stop_codons = ['TAA', 'TAG', 'TGA']
    return codon in stop_codons

def map_sae_region_to_genomic(start_window, end_window, mapping_df, strand, window_size=9, max_window_in_gene=None):
    """Map SAE region to genomic coordinates, excluding stop codon if needed"""
    start_codon = start_window
    end_codon = end_window + window_size - 1
    log(f"    Mapping windows {start_window}-{end_window} (codons {start_codon}-{end_codon})")
    region_rows = mapping_df[(mapping_df['Alignment_Codon'] >= start_codon) & (mapping_df['Alignment_Codon'] <= end_codon)].copy()
    if len(region_rows) == 0:
        raise ValueError(f"No mapping found for codons {start_codon}-{end_codon}")
    non_gap_rows = region_rows[region_rows['Alignment_Base'] != '-'].copy()
    total_positions = len(region_rows)
    gap_positions = len(region_rows) - len(non_gap_rows)
    log(f"      Total positions: {total_positions} ({len(non_gap_rows)} bases, {gap_positions} gaps)")
    if len(non_gap_rows) == 0:
        log(f"      WARNING: Entire region is gaps! Skipping this SAE.")
        return None
    if max_window_in_gene is not None and end_window == max_window_in_gene:
        log(f"      SAE includes last window ({max_window_in_gene}) - checking for stop codon")
        valid_ccds = mapping_df[mapping_df['CCDS_Index'] != '-'].copy()
        valid_ccds['CCDS_Index_int'] = valid_ccds['CCDS_Index'].astype(int)
        max_ccds_index = valid_ccds['CCDS_Index_int'].max()
        stop_codon_indices = [max_ccds_index - 2, max_ccds_index - 1, max_ccds_index]
        log(f"      Last 3 CCDS indices: {stop_codon_indices}")
        stop_codon_rows = mapping_df[mapping_df['CCDS_Index'].astype(str).str.replace('-', '0').astype(int).isin(stop_codon_indices)]
        if len(stop_codon_rows) > 0:
            stop_bases = stop_codon_rows['Alignment_Base'].tolist()
            stop_bases_clean = [b for b in stop_bases if b != '-']
            log(f"      Bases at last 3 positions: {stop_bases_clean}")
            if len(stop_bases_clean) == 3:
                is_stop = verify_stop_codon(stop_bases_clean)
                codon_sequence = ''.join(stop_bases_clean).upper()
                if is_stop:
                    log(f"      Confirmed stop codon: {codon_sequence}")
                    log(f"      Excluding last 3 CCDS positions from SAE region")
                    before_count = len(non_gap_rows)
                    non_gap_rows = non_gap_rows[~non_gap_rows['CCDS_Index'].astype(str).str.replace('-', '0').astype(int).isin(stop_codon_indices)]
                    after_count = len(non_gap_rows)
                    excluded = before_count - after_count
                    log(f"      Excluded {excluded} bases ({codon_sequence})")
                    if len(non_gap_rows) == 0:
                        log(f"      WARNING: Only stop codon in this SAE region! Skipping.")
                        return None
                else:
                    log(f"      Last 3 bases ({codon_sequence}) are NOT a stop codon - keeping all positions")
            else:
                log(f"      Could not verify stop codon (found {len(stop_bases_clean)} bases, expected 3)")
        else:
            log(f"      No bases found at last 3 CCDS positions")
    genomic_positions = non_gap_rows['CCDS_Base_Location'].astype(int).tolist()
    exon_coords = non_gap_rows['CCDS_Coordinates'].unique()
    exon_coords = [e for e in exon_coords if e != '-']
    log(f"      Spans {len(exon_coords)} exon(s): {exon_coords}")
    blocks = []
    current_block_start = genomic_positions[0]
    current_block_end = genomic_positions[0]
    current_exon = non_gap_rows.iloc[0]['CCDS_Coordinates']
    for i in range(1, len(non_gap_rows)):
        pos = int(non_gap_rows.iloc[i]['CCDS_Base_Location'])
        exon = non_gap_rows.iloc[i]['CCDS_Coordinates']
        if exon == current_exon:
            current_block_end = pos
        else:
            blocks.append({'start': int(current_block_start), 'end': int(current_block_end), 'exon': current_exon})
            current_block_start = pos
            current_block_end = pos
            current_exon = exon
    blocks.append({'start': int(current_block_start), 'end': int(current_block_end), 'exon': current_exon})
    log(f"      Genomic blocks: {len(blocks)}")
    for i, block in enumerate(blocks, 1):
        log(f"        Block {i}: {block['start']}-{block['end']} (exon {block['exon']})")
    return {'blocks': blocks, 'strand': strand, 'num_windows': end_window - start_window + 1}

def blocks_to_bed12(blocks, chrom, ccds_id, strand, score=0):
    """Convert genomic blocks to BED12 format"""
    if not blocks:
        return None
    if strand == '-':
        for block in blocks:
            if block['start'] > block['end']:
                block['start'], block['end'] = block['end'], block['start']
        blocks_sorted = sorted(blocks, key=lambda b: b['start'])
    else:
        blocks_sorted = sorted(blocks, key=lambda b: b['start'])
    all_starts = [int(b['start']) for b in blocks_sorted]
    all_ends = [int(b['end']) for b in blocks_sorted]
    overall_start_1based = min(all_starts)
    overall_end_1based = max(all_ends)
    chrom_start = overall_start_1based - 1
    chrom_end = overall_end_1based
    block_count = len(blocks_sorted)
    block_sizes = []
    chrom_starts = []
    for block in blocks_sorted:
        block_start_0based = int(block['start']) - 1
        block_end_0based = int(block['end'])
        size = block_end_0based - block_start_0based
        block_sizes.append(size)
        offset = block_start_0based - chrom_start
        chrom_starts.append(offset)
    block_sizes_str = ','.join(map(str, block_sizes))
    chrom_starts_str = ','.join(map(str, chrom_starts))
    bed_line = [chrom, str(chrom_start), str(chrom_end), ccds_id, str(score), strand, str(chrom_start), str(chrom_end), "255,0,0", str(block_count), block_sizes_str, chrom_starts_str]
    return '\t'.join(bed_line)

def process_gene(fresco_file, mapping_file):
    """Process one gene: extract SAEs and convert to BED12"""
    log("="*80)
    log(f"Processing gene: {os.path.basename(fresco_file)}")
    log("="*80)
    fresco_df = load_fresco_output(fresco_file)
    mapping_df, chrom, strand, ccds_id, gene_name = load_mapping_table(mapping_file)
    fresco_basename = os.path.basename(fresco_file).replace('_fresco.txt', '')
    mapping_basename = os.path.basename(mapping_file).replace('_map.txt', '')
    if fresco_basename != mapping_basename:
        log(f"WARNING: Filename mismatch!")
        log(f"   FRESCo: {fresco_basename}")
        log(f"   Mapping: {mapping_basename}")
    result = find_sae_windows(fresco_df)
    if not result:
        log("  No SAE windows - skipping gene")
        return []
    sae_windows, sae_df = result
    
    # USE NEW MERGING FUNCTION HERE!
    sae_regions = merge_overlapping_windows(sae_windows, window_size=9)
    
    if not sae_regions:
        log("  No SAE regions after merging - skipping gene")
        return []
    max_window_in_gene = fresco_df.iloc[:, 0].max()
    log(f"  Max window in gene: {max_window_in_gene}")
    bed_lines = []
    for region_num, (start_win, end_win) in enumerate(sae_regions, 1):
        log(f"  Processing SAE region {region_num}/{len(sae_regions)}:")
        genomic_info = map_sae_region_to_genomic(start_win, end_win, mapping_df, strand, window_size=9, max_window_in_gene=max_window_in_gene)
        if genomic_info is None:
            log(f"    WARNING: Could not map region - skipping")
            continue
        bed_line = blocks_to_bed12(genomic_info['blocks'], chrom, ccds_id, strand, score=0)
        if bed_line:
            bed_lines.append(bed_line)
            log(f"    Generated BED entry")
    log(f"Completed: {len(bed_lines)} SAE regions for {gene_name} ({ccds_id})")
    return bed_lines

def read_existing_bed(bed_file):
    """Read existing BED file into dictionary keyed by CCDS ID"""
    if not os.path.exists(bed_file):
        log(f"BED file doesn't exist - will create new: {bed_file}")
        return {}
    log(f"Reading existing BED file: {bed_file}")
    bed_entries = defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                ccds_id = parts[3]
                bed_entries[ccds_id].append(line.strip())
    total_entries = sum(len(v) for v in bed_entries.values())
    log(f"  Found {total_entries} existing entries for {len(bed_entries)} CCDS IDs")
    return bed_entries

def write_bed_file(bed_file, all_entries):
    """Write all BED entries to file (sorted)"""
    log(f"Writing BED file: {bed_file}")
    all_lines = []
    for ccds_id in all_entries.keys():
        all_lines.extend(all_entries[ccds_id])
    def sort_key(line):
        parts = line.split('\t')
        chrom = parts[0]
        start = int(parts[1])
        return (chrom, start)
    all_lines.sort(key=sort_key)
    with open(bed_file, 'w') as out:
        for line in all_lines:
            out.write(line + '\n')
    log(f"  Wrote {len(all_lines)} BED entries (sorted)")

def main():
    start_time = datetime.now()
    log("="*80)
    log("SAE Region BED12 Generator - VERSION 2 (Codon-Overlap Merging)")
    log("="*80)
    if len(sys.argv) < 5:
        print("Usage: python chr_SAE_generation_v2.py <fresco_file> <mapping_file> <output_bed> <mode>")
        print("Arguments:")
        print("  fresco_file  : FRESCo output with p_adj")
        print("  mapping_file : Mapping table")
        print("  output_bed   : Output BED12 file path")
        print("  mode         : 'skip' or 'overwrite'")
        sys.exit(1)
    fresco_file = sys.argv[1]
    mapping_file = sys.argv[2]
    output_bed = sys.argv[3]
    mode = sys.argv[4]
    if mode not in ['skip', 'overwrite']:
        log(f"ERROR: Invalid mode '{mode}'")
        sys.exit(1)
    log(f"Input FRESCo: {fresco_file}")
    log(f"Input Mapping: {mapping_file}")
    log(f"Output BED: {output_bed}")
    log(f"Mode: {mode.upper()}")
    log("")
    if not os.path.exists(fresco_file):
        log(f"ERROR: FRESCo file not found: {fresco_file}")
        sys.exit(1)
    if not os.path.exists(mapping_file):
        log(f"ERROR: Mapping file not found: {mapping_file}")
        sys.exit(1)
    try:
        new_bed_lines = process_gene(fresco_file, mapping_file)
        _, _, _, ccds_id, gene_name = load_mapping_table(mapping_file)
        if not new_bed_lines:
            log(f"NO_SAE_FOUND: {gene_name} ({ccds_id}) - No SAE regions in this gene")
            log("="*80)
            sys.exit(0)
        existing_entries = read_existing_bed(output_bed)
        if ccds_id in existing_entries:
            if mode == 'skip':
                log(f"SKIPPING: {ccds_id} already in BED file")
                log("="*80)
                sys.exit(0)
            else:
                log(f"OVERWRITING: Replacing {len(existing_entries[ccds_id])} existing entries for {ccds_id}")
                existing_entries[ccds_id] = new_bed_lines
        else:
            log(f"ADDING: New CCDS ID {ccds_id} ({len(new_bed_lines)} regions)")
            existing_entries[ccds_id] = new_bed_lines
        write_bed_file(output_bed, existing_entries)
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        log("="*80)
        log(f"SUCCESS! Completed in {duration:.2f} seconds")
        log(f"Output: {output_bed}")
        log("="*80)
    except Exception as e:
        log(f"ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()

