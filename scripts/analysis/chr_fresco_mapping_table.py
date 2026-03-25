#!/usr/bin/env python3
"""
Create mapping table for gene-level FASTA to genomic coordinates
Maps alignment positions to genomic coordinates accounting for gaps
"""

import sys
import os
from datetime import datetime
from collections import defaultdict

def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into dictionary"""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        if ' ' in item:
            key, value = item.split(' ', 1)
            attrs[key] = value.strip('"')
    return attrs

def load_gtf_cds(gtf_file, gene_name, ccds_id):
    """
    Load CDS coordinates for a specific gene and CCDS ID from GTF
    Returns: list of (chrom, start, end, strand) tuples, sorted appropriately
    """
    log(f"Loading GTF file: {gtf_file}")
    log(f"Searching for gene: {gene_name}, CCDS: {ccds_id}")
    
    cds_list = []
    found = False
    
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, frame, attributes = parts
            
            if feature != 'CDS':
                continue
            
            attrs = parse_gtf_attributes(attributes)
            
            # Check if this matches our gene and CCDS
            gtf_gene = attrs.get('gene_name', '')
            gtf_ccds = attrs.get('transcript_id', '')
            
            if gtf_gene == gene_name and gtf_ccds == ccds_id:
                cds_list.append({
                    'chrom': chrom,
                    'start': int(start),  # GTF is 1-based
                    'end': int(end),
                    'strand': strand
                })
                found = True
    
    if not found:
        raise ValueError(f"No CDS found for gene '{gene_name}' with CCDS ID '{ccds_id}' in GTF")
    
    if not cds_list:
        raise ValueError(f"No CDS coordinates found")
    
    # Get strand (should be same for all)
    strand = cds_list[0]['strand']
    chrom = cds_list[0]['chrom']
    
    # Sort by genomic position (always by start position first)
    cds_list.sort(key=lambda x: x['start'])
    
    log(f"Found {len(cds_list)} CDS exons on {chrom} ({strand} strand)")
    for i, cds in enumerate(cds_list, 1):
        log(f"  Exon {i}: {cds['start']}-{cds['end']} ({cds['end']-cds['start']+1} bp)")
    
    return cds_list, chrom, strand

def load_hg38_sequence(fasta_file):
    """
    Load hg38 sequence from gene-level FASTA file
    Returns: sequence string (with gaps)
    """
    log(f"Loading FASTA file: {fasta_file}")
    
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Save previous sequence
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                
                current_header = line[1:].strip()  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_header is not None:
            sequences[current_header] = ''.join(current_seq)
    
    # Find hg38 sequence
    if 'hg38' not in sequences:
        raise ValueError(f"No 'hg38' sequence found in FASTA. Available: {list(sequences.keys())}")
    
    hg38_seq = sequences['hg38']
    
    total_bases = len(hg38_seq)
    gap_count = hg38_seq.count('-')
    real_bases = total_bases - gap_count
    
    log(f"hg38 sequence loaded: {total_bases} positions ({real_bases} bases, {gap_count} gaps)")
    log(f"Total codons in alignment: {total_bases // 3}")
    
    return hg38_seq

def create_mapping_table(gene_name, ccds_id, hg38_seq, cds_list, chrom, strand, output_file):
    """
    Create the mapping table
    """
    log("Creating mapping table...")
    
    # Calculate total CDS length from GTF
    total_cds_length = sum(cds['end'] - cds['start'] + 1 for cds in cds_list)
    log(f"Total CDS length from GTF: {total_cds_length} bp")
    
    # Count real bases (no gaps) in hg38 sequence
    real_bases_in_seq = len([b for b in hg38_seq if b != '-'])
    log(f"Real bases in hg38 sequence: {real_bases_in_seq} bp")
    
    if real_bases_in_seq != total_cds_length:
        log(f"WARNING: CDS length mismatch! GTF: {total_cds_length}, FASTA: {real_bases_in_seq}")
    
    # Prepare coordinate mapping based on strand
    if strand == '+':
        # Forward strand: use coordinates as-is, ascending order
        coord_list = []
        for cds in cds_list:
            for pos in range(cds['start'], cds['end'] + 1):
                coord_list.append((pos, f"{cds['start']}-{cds['end']}"))
    else:
        # Reverse strand: reverse exon order, descending within each
        cds_list_reversed = list(reversed(cds_list))
        coord_list = []
        for cds in cds_list_reversed:
            for pos in range(cds['end'], cds['start'] - 1, -1):  # Descending
                coord_list.append((pos, f"{cds['start']}-{cds['end']}"))
    
    log(f"Generated {len(coord_list)} coordinate positions")
    
    # Verify we have enough coordinates
    if len(coord_list) != real_bases_in_seq:
        raise ValueError(f"Coordinate count ({len(coord_list)}) doesn't match real bases ({real_bases_in_seq})")
    
    # Create table rows
    rows = []
    ccds_index = 1
    coord_index = 0
    
    for aln_index in range(len(hg38_seq)):
        base = hg38_seq[aln_index]
        alignment_codon = aln_index // 3
        alignment_index_1based = aln_index + 1
        
        if base == '-':
            # Gap
            row = {
                'Chromosome': chrom,
                'Alignment_Codon': alignment_codon,
                'Alignment_Index': alignment_index_1based,
                'Alignment_Base': '-',
                'CCDS_Index': '-',
                'CCDS_Base_Location': '-',
                'CCDS_Coordinates': '-',
                'CCDS_ID': ccds_id,
                'Gene_Name': gene_name,
                'Strand': strand
            }
        else:
            # Real base
            if coord_index >= len(coord_list):
                raise ValueError(f"Ran out of coordinates at alignment position {aln_index}")
            
            genomic_pos, coord_range = coord_list[coord_index]
            
            row = {
                'Chromosome': chrom,
                'Alignment_Codon': alignment_codon,
                'Alignment_Index': alignment_index_1based,
                'Alignment_Base': base,
                'CCDS_Index': ccds_index,
                'CCDS_Base_Location': genomic_pos,
                'CCDS_Coordinates': coord_range,
                'CCDS_ID': ccds_id,
                'Gene_Name': gene_name,
                'Strand': strand
            }
            
            ccds_index += 1
            coord_index += 1
        
        rows.append(row)
    
    log(f"Created {len(rows)} rows in mapping table")
    
    # Write output
    log(f"Writing output to: {output_file}")
    
    header = [
        'Chromosome',
        'Alignment_Codon',
        'Alignment_Index',
        'Alignment_Base',
        'CCDS_Index',
        'CCDS_Base_Location',
        'CCDS_Coordinates',
        'CCDS_ID',
        'Gene_Name',
        'Strand'
    ]
    
    with open(output_file, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for row in rows:
            out.write('\t'.join(str(row[h]) for h in header) + '\n')
    
    log(f"Successfully wrote mapping table with {len(rows)} rows")
    
    # Summary statistics
    gap_rows = sum(1 for r in rows if r['Alignment_Base'] == '-')
    base_rows = len(rows) - gap_rows
    log(f"Summary: {base_rows} bases, {gap_rows} gaps, {len(rows)} total positions")

def check_output_exists(output_file, skip_existing=False):
    """
    Check if output file already exists and handle accordingly
    Returns: True if should process, False if should skip
    """
    if os.path.exists(output_file):
        log(f"⚠️  WARNING: Output file already exists!")
        log(f"   File: {output_file}")
        
        if skip_existing:
            log(f"   Action: SKIPPING (skip_existing=True)")
            return False
        else:
            log(f"   Action: OVERWRITING")
            return True
    
    return True

def main():
    start_time = datetime.now()
    log("="*80)
    log("FRESCo Mapping Table Generator")
    log("="*80)
    
    # Parse arguments
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python chr_fresco_mapping_table.py <fasta_file> <gtf_file> <output_dir> [skip_existing]")
        print()
        print("Arguments:")
        print("  fasta_file     : Path to gene-level FASTA file (GENENAME_CCDSID.fa)")
        print("  gtf_file       : Path to GTF file with CDS annotations")
        print("  output_dir     : Directory to save mapping table")
        print("  skip_existing  : Optional. Set to 'skip' to skip existing files")
        print()
        print("Example:")
        print("  python chr_fresco_mapping_table.py A3GALT2_CCDS60080.1.fa CCDS.gtf output/")
        print("  python chr_fresco_mapping_table.py A3GALT2_CCDS60080.1.fa CCDS.gtf output/ skip")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    gtf_file = sys.argv[2]
    output_dir = sys.argv[3]
    skip_existing = (len(sys.argv) == 5 and sys.argv[4] == 'skip')
    
    log(f"Input FASTA: {fasta_file}")
    log(f"Input GTF: {gtf_file}")
    log(f"Output directory: {output_dir}")
    log(f"Skip existing: {skip_existing}")
    
    # Parse gene name and CCDS ID from filename
    basename = os.path.basename(fasta_file)
    if not basename.endswith('.fa'):
        raise ValueError(f"FASTA file must end with .fa: {basename}")
    
    name_part = basename[:-3]  # Remove .fa
    
    # Split on underscore: GENENAME_CCDSID
    if '_' not in name_part:
        raise ValueError(f"Filename must be GENENAME_CCDSID.fa format: {basename}")
    
    parts = name_part.split('_', 1)
    gene_name = parts[0]
    ccds_id = parts[1]
    
    log(f"Parsed gene name: {gene_name}")
    log(f"Parsed CCDS ID: {ccds_id}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Output file name
    output_file = os.path.join(output_dir, f"{gene_name}_{ccds_id}_map.txt")
    
    # Check if output exists
    if not check_output_exists(output_file, skip_existing):
        log("="*80)
        log("SKIPPED - File already processed")
        log("="*80)
        sys.exit(0)
    
    try:
        # Load data
        cds_list, chrom, strand = load_gtf_cds(gtf_file, gene_name, ccds_id)
        hg38_seq = load_hg38_sequence(fasta_file)
        
        # Create mapping table
        create_mapping_table(gene_name, ccds_id, hg38_seq, cds_list, chrom, strand, output_file)
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        log("="*80)
        log(f"✓ SUCCESS! Completed in {duration:.2f} seconds")
        log(f"Output: {output_file}")
        log("="*80)
        
    except Exception as e:
        log(f"✗ ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
