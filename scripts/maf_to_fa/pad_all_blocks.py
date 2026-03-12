#!/usr/bin/env python3
"""
Pad MafFilter block FASTA files so that each block has all species,
and write a summary file of hg38 headers for all blocks.

Example usage (from sbatch):

    python pad_all_blocks.py \
        --species /path/to/species_120.txt \
        --input-dir /path/to/chr1_blocks \
        --prefix chr1_blocks_ \
        --output-dir /path/to/chr1_blocks_padded \
        --header-out /path/to/chr1_hg38_headers.txt \
        --ref-species hg38
"""

import os
import sys
import argparse
import re
from collections import OrderedDict


def parse_args():
    p = argparse.ArgumentParser(
        description="Pad MafFilter block FASTAs and log hg38 headers."
    )
    p.add_argument(
        "--species", required=True,
        help="Text file with 1 species name per line (must include ref species, e.g. hg38)."
    )
    p.add_argument(
        "--input-dir", required=True,
        help="Directory containing block FASTA files."
    )
    p.add_argument(
        "--prefix", required=True,
        help="Filename prefix for blocks (e.g. chr1_blocks_). Only files starting with this prefix "
             "and ending in .fa will be processed."
    )
    p.add_argument(
        "--output-dir", required=True,
        help="Directory to write padded FASTA files."
    )
    p.add_argument(
        "--header-out", required=True,
        help="Path to write hg38 header summary TXT file."
    )
    p.add_argument(
        "--ref-species", default="hg38",
        help="Reference species name (default: hg38)."
    )
    return p.parse_args()


def load_species_list(path):
    species = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            species.append(line)
    if not species:
        sys.stderr.write(f"[ERROR] Species list {path} is empty.\n")
        sys.exit(1)
    sys.stderr.write(f"[INFO] Loaded {len(species)} species from {path}\n")
    return species


def natural_block_sort_key(filename):
    """
    Extract a numeric index from a filename like 'chr1_blocks_123.fa'
    so sorting is 1,2,3,...,10,11,... not 1,10,100,11,...
    """
    m = re.search(r"(\d+)(?=\.fa$)", filename)
    if m:
        return int(m.group(1))
    return 0


def list_block_files(input_dir, prefix):
    files = [
        f for f in os.listdir(input_dir)
        if f.startswith(prefix) and f.endswith(".fa")
    ]
    if not files:
        sys.stderr.write(f"[ERROR] No files found in {input_dir} with prefix '{prefix}' and suffix '.fa'\n")
        sys.exit(1)

    files.sort(key=natural_block_sort_key)
    sys.stderr.write(f"[INFO] Found {len(files)} block files to process.\n")
    return [os.path.join(input_dir, f) for f in files]


def read_block_fasta(path, ref_species):
    """
    Read a single multi-FASTA block file and return:
      - headers: OrderedDict[species] = header_string (without '>')
      - seqs   : OrderedDict[species] = full_sequence_string
      - ref_header: str or None (hg38 header)
    We assume there is exactly ONE block per file.
    """
    headers = OrderedDict()
    seq_chunks = OrderedDict()
    current_species = None
    ref_header = None

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                # flush previous sequence if any
                if current_species is not None:
                    seq_chunks[current_species] = "".join(seq_chunks[current_species])

                header = line[1:]
                species = header.split("-", 1)[0]

                headers[species] = header
                seq_chunks[species] = []

                if species == ref_species and ref_header is None:
                    ref_header = header

                current_species = species
            else:
                if current_species is None:
                    continue
                seq_chunks[current_species].append(line)

    # flush last
    if current_species is not None:
        seq_chunks[current_species] = "".join(seq_chunks[current_species])

    return headers, seq_chunks, ref_header


def pad_block(headers, seqs, species_list, ref_species):
    """
    Given a single block (headers + seqs for some species),
    return new ordered headers and seqs for ALL species in species_list.

    Missing species get:
      - header = hg38 header with species name replaced
      - seq    = all gaps of length len(hg38_seq)
    """
    if ref_species not in seqs:
        sys.stderr.write(f"[WARN] Block missing ref species '{ref_species}', skipping padding.\n")
        return headers, seqs  # nothing we can do

    ref_header = headers[ref_species]
    ref_seq = seqs[ref_species]
    L = len(ref_seq)

    new_headers = OrderedDict()
    new_seqs = OrderedDict()

    # First keep existing species in original order
    for sp, h in headers.items():
        new_headers[sp] = h
        new_seqs[sp] = seqs[sp]

    # Then add missing species in species_list order
    for sp in species_list:
        if sp in new_headers:
            continue
        # build new header by replacing the species name at the start
        # (before the first '-')
        # ref_header like: "hg38-chr1(+)/65565-65574"
        parts = ref_header.split("-", 1)
        if len(parts) == 2:
            new_header = sp + "-" + parts[1]
        else:
            # fallback, just prefix species
            new_header = sp + "-" + ref_header
        new_headers[sp] = new_header
        new_seqs[sp] = "-" * L

    return new_headers, new_seqs


def write_block_fasta(path, headers, seqs):
    """
    Write a single multi-FASTA block file given headers + seqs (OrderedDict).
    """
    with open(path, "w") as out:
        for sp, header in headers.items():
            seq = seqs[sp]
            out.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")


def pad_file(input_fasta, output_fasta, species_list, ref_species, header_log_handle):
    """
    Pad all blocks in a single MafFilter FASTA file (one block per file),
    write padded FASTA, and log the hg38 header.
    """
    sys.stderr.write(f"[INFO] Padding file: {os.path.basename(input_fasta)}\n")

    headers, seqs, ref_header = read_block_fasta(input_fasta, ref_species)

    if ref_header is None:
        sys.stderr.write(f"[WARN] No {ref_species} header found in {input_fasta}; "
                         "writing file unchanged and logging NA.\n")
        # still write unchanged file
        write_block_fasta(output_fasta, headers, seqs)
        header_log_handle.write(f"{os.path.basename(input_fasta)}\tNA\n")
        return

    # Log hg38 header
    header_log_handle.write(f"{os.path.basename(input_fasta)}\t{ref_header}\n")

    # Pad missing species
    new_headers, new_seqs = pad_block(headers, seqs, species_list, ref_species)

    # Write padded FASTA
    write_block_fasta(output_fasta, new_headers, new_seqs)


def main():
    args = parse_args()

    species_list = load_species_list(args.species)
    if args.ref_species not in species_list:
        sys.stderr.write(
            f"[WARN] ref-species '{args.ref_species}' not found in species list. "
            f"It should be included for consistent padding.\n"
        )

    os.makedirs(args.output_dir, exist_ok=True)
    header_out_dir = os.path.dirname(args.header_out)
    if header_out_dir:
        os.makedirs(header_out_dir, exist_ok=True)

    files = list_block_files(args.input_dir, args.prefix)

    # Open header summary file and write header line
    with open(args.header_out, "w") as header_log:
        header_log.write("block_file\thg38_header\n")

        for i, infile in enumerate(files, start=1):
            base = os.path.basename(infile)
            name, ext = os.path.splitext(base)
            outfile = os.path.join(args.output_dir, f"{name}_padded{ext}")

            sys.stderr.write(
                f"[INFO] ({i}/{len(files)}) Processing {base} -> {os.path.basename(outfile)}\n"
            )
            pad_file(infile, outfile, species_list, args.ref_species, header_log)

    sys.stderr.write("[INFO] All files padded successfully.\n")


if __name__ == "__main__":
    main()
