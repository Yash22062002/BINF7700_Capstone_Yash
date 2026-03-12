#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict

# -----------------------------
# 1. Parse arguments
# -----------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Aggregate padded CDS blocks into (gene,CCDS) multi-species FASTA."
    )
    p.add_argument(
        "--gtf", required=True,
        help="CCDS CDS-only GTF (e.g. CCDS_public_hg38.CDS.gtf)"
    )
    p.add_argument(
        "--blocks-dir", required=True,
        help="Directory with padded block FASTAs for one chromosome"
    )
    p.add_argument(
        "--out-fasta-dir", required=True,
        help="Directory to write per-(gene,CCDS) FASTA files"
    )
    p.add_argument(
        "--out-tsv", required=True,
        help="Path for summary TSV describing each written (gene,CCDS) FASTA"
    )
    p.add_argument(
        "--chrom", required=True,
        help="Chromosome name to restrict to (e.g. 'chr1')"
    )
    p.add_argument(
        "--hg-species", default="hg38",
        help="Reference species prefix in headers (default: hg38)"
    )
    return p.parse_args()

# -----------------------------
# 2. GTF loading
# -----------------------------

def parse_attributes(attr_str):
    """
    Parse 9th column of GTF into a dict: key -> value (no quotes).
    """
    attrs = {}
    for field in attr_str.split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            continue
        key, val = field.split(" ", 1)
        val = val.strip().strip('"')
        attrs[key] = val
    return attrs

def load_ccds_gtf(gtf_path, chrom):
    """
    Load CCDS CDS GTF rows for a single chromosome.

    Returns:
      transcripts: dict keyed by (gene_name, transcript_id, chrom, strand) with:
          {
            'gene_name': ...,
            'transcript_id': ...,
            'chrom': ...,
            'strand': '+' or '-',
            'cds_intervals': [(start,end), ...]
          }
      cds_intervals_by_chrom: dict chrom -> list of
          (start,end,gene_name,transcript_id,strand)
    """
    transcripts = {}
    cds_intervals_by_chrom = defaultdict(list)
    cds_count = 0

    print(f"[INFO] Loading GTF: {gtf_path}")
    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            g_chrom, source, feature, start, end, score, strand, frame, attrs = parts

            # Restrict to requested chromosome
            if g_chrom != chrom:
                continue

            # Only CDS
            if feature != "CDS":
                continue

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            attr_dict = parse_attributes(attrs)
            gene_name = attr_dict.get("gene_name") or attr_dict.get("gene_id")
            transcript_id = attr_dict.get("transcript_id")
            if gene_name is None or transcript_id is None:
                continue

            cds_count += 1

            key = (gene_name, transcript_id, g_chrom, strand)
            if key not in transcripts:
                transcripts[key] = {
                    "gene_name": gene_name,
                    "transcript_id": transcript_id,
                    "chrom": g_chrom,
                    "strand": strand,
                    "cds_intervals": []
                }

            transcripts[key]["cds_intervals"].append((start_i, end_i))
            cds_intervals_by_chrom[g_chrom].append(
                (start_i, end_i, gene_name, transcript_id, strand)
            )

    # Sort CDS intervals per transcript and per chromosome
    for key in transcripts:
        transcripts[key]["cds_intervals"].sort(key=lambda x: x[0])
    for c in cds_intervals_by_chrom:
        cds_intervals_by_chrom[c].sort(key=lambda x: x[0])

    print(
        f"[INFO] Loaded {cds_count} CDS rows "
        f"for {len(transcripts)} (gene,CCDS) transcripts on {chrom}."
    )
    return transcripts, cds_intervals_by_chrom

# -----------------------------
# 3. FASTA helpers
# -----------------------------

def read_fasta_entries(path):
    """
    Simple FASTA reader: returns list of (header_without_>, sequence_string).
    """
    entries = []
    header = None
    seq_chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_chunks)))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        entries.append((header, "".join(seq_chunks)))

    return entries

def parse_block_header(header):
    """
    Parse header like: hg38-chr1(+)/65565-65574
    Returns: (species, chrom, strand, start, end)
    """
    if "-" not in header:
        return None, None, None, None, None
    species, rest = header.split("-", 1)

    if "/" not in rest:
        return species, None, None, None, None

    coord_part, pos_part = rest.split("/", 1)

    if "(" in coord_part and ")" in coord_part:
        chrom = coord_part.split("(", 1)[0]
        strand = coord_part.split("(", 1)[1].split(")", 1)[0]
    else:
        chrom = coord_part
        strand = "."

    if "-" not in pos_part:
        return species, chrom, strand, None, None

    try:
        s_str, e_str = pos_part.split("-", 1)
        start = int(s_str)
        end = int(e_str)
    except ValueError:
        start = end = None

    return species, chrom, strand, start, end

def find_transcripts_for_block(chrom, start, end, cds_intervals_by_chrom):
    """
    For a block on (chrom, start, end), return list of (gene_name, transcript_id, strand)
    whose CDS intervals contain this block's start coordinate:
       cds_start <= start <= cds_end
    """
    hits = []
    if chrom not in cds_intervals_by_chrom:
        return hits

    intervals = cds_intervals_by_chrom[chrom]  # sorted by cds_start

    for cds_start, cds_end, gene_name, tid, strand in intervals:
        if cds_start > start:
            # later CDS all start after block_start -> no more hits
            break
        if cds_end < start:
            # this CDS ends before block_start
            continue
        # cds_start <= start <= cds_end
        hits.append((gene_name, tid, strand))

    return hits

# -----------------------------
# 4. Main aggregation
# -----------------------------

def aggregate_blocks(args):
    # Load GTF
    transcripts, cds_intervals_by_chrom = load_ccds_gtf(args.gtf, args.chrom)

    # Accumulator for per-transcript concatenation
    # key = (gene_name, transcript_id, chrom, strand)
    # value = { "segments": [ {start,end,seqs}, ... ], "seen_coords": set((chrom,start,end)) }
    transcript_seqs = {}

    species_order = None

    # Prepare output dirs
    os.makedirs(args.out_fasta_dir, exist_ok=True)
    tsv_dir = os.path.dirname(args.out_tsv)
    if tsv_dir:
        os.makedirs(tsv_dir, exist_ok=True)

    # Find all block FASTAs
    block_files = sorted(
        f for f in os.listdir(args.blocks_dir)
        if f.endswith(".fa") or f.endswith(".fasta")
    )
    if not block_files:
        print(f"[ERROR] No FASTA files found in {args.blocks_dir}")
        sys.exit(1)

    print(f"[INFO] Found {len(block_files)} padded block FASTAs in {args.blocks_dir}")

    total_blocks = 0          # hg38 blocks on this chrom processed
    mapped_blocks = 0
    unmapped_blocks = 0

    # Scan blocks
    for idx, fname in enumerate(block_files, start=1):
        fpath = os.path.join(args.blocks_dir, fname)
        entries = read_fasta_entries(fpath)
        if not entries:
            continue

        # Infer species order from the first non-empty block
        if species_order is None:
            seen_sp = set()
            species_order = []
            for h, _ in entries:
                sp = h.split("-", 1)[0]
                if sp not in seen_sp:
                    seen_sp.add(sp)
                    species_order.append(sp)
            print(f"[INFO] Species order inferred from first block: {len(species_order)} species.")
            print(f"[INFO] Example species list: {species_order[:5]} ...")

        # Build mapping species -> seq for this block
        block_seqs = {}
        hg38_info = None
        for h, seq in entries:
            sp, b_chrom, b_strand, b_start, b_end = parse_block_header(h)
            block_seqs[sp] = seq
            if sp == args.hg_species:
                hg38_info = (b_chrom, b_start, b_end)

        if hg38_info is None:
            # No hg38 sequence found in this block
            unmapped_blocks += 1
            continue

        chrom, start, end = hg38_info

        # Restrict to requested chromosome
        if chrom != args.chrom:
            continue

        total_blocks += 1

        # Which transcripts use this block?
        hits = find_transcripts_for_block(chrom, start, end, cds_intervals_by_chrom)
        if not hits:
            unmapped_blocks += 1
        else:
            mapped_blocks += 1
            for gene_name, tid, strand in hits:
                t_key = (gene_name, tid, chrom, strand)
                if t_key not in transcript_seqs:
                    transcript_seqs[t_key] = {
                        "segments": [],
                        "seen_coords": set()
                    }

                seg_key = (chrom, start, end)
                if seg_key in transcript_seqs[t_key]["seen_coords"]:
                    # duplicate block for this transcript; skip
                    continue

                transcript_seqs[t_key]["seen_coords"].add(seg_key)
                transcript_seqs[t_key]["segments"].append(
                    {"start": start, "end": end, "seqs": block_seqs}
                )

        if idx % 1000 == 0:
            print(
                f"[INFO] Processed {idx} block files "
                f"(hg38 blocks considered on {args.chrom}: {total_blocks}, "
                f"mapped: {mapped_blocks}, unmapped: {unmapped_blocks})"
            )

    print("[INFO] Finished scanning padded blocks.")
    print(f"[INFO] Total hg38 blocks considered on {args.chrom}: {total_blocks}")
    print(f"[INFO] Mapped blocks: {mapped_blocks}")
    print(f"[INFO] Unmapped blocks: {unmapped_blocks}")
    print(f"[INFO] Transcripts with ≥1 block: {len(transcript_seqs)}")

    # -------------------------
    # 5. Write per-transcript FASTA + TSV summary
    # -------------------------

    tsv_lines = []
    header_cols = [
        "gene_name",
        "transcript_id",
        "chrom",
        "strand",
        "coord_order",
        "n_cds_exons",
        "n_blocks_used",
        "cds_start_min",
        "cds_end_max",
        "fasta_path",
    ]
    tsv_lines.append("\t".join(header_cols))

    written = 0

    for t_key, t_data in transcript_seqs.items():
        gene_name, tid, chrom, strand = t_key
        segments = t_data["segments"]
        if not segments:
            continue

        # Sort segments according to transcript strand
        if strand == "+":
            coord_order = "ascending"
            segments_sorted = sorted(segments, key=lambda s: s["start"])
        else:
            coord_order = "descending"
            segments_sorted = sorted(segments, key=lambda s: s["start"], reverse=True)

        # Build full sequences per species
        full_seqs = {sp: [] for sp in species_order}
        for seg in segments_sorted:
            seg_seqs = seg["seqs"]
            for sp in species_order:
                seq = seg_seqs.get(sp)
                if seq is None:
                    # If somehow missing, pad with gaps of hg38 segment length
                    hg_seq = seg_seqs.get(args.hg_species, "")
                    seq = "-" * len(hg_seq)
                full_seqs[sp].append(seq)

        # Concatenate chunks
        for sp in full_seqs:
            full_seqs[sp] = "".join(full_seqs[sp])

        # Sanity check: all species same length
        lengths = {len(s) for s in full_seqs.values()}
        if len(lengths) > 1:
            print(
                f"[WARN] Transcript {gene_name}|{tid} has unequal sequence lengths across species: {lengths}. Skipping."
            )
            continue

        # Output FASTA name: gene_transcript.fa
        safe_gene = gene_name.replace(" ", "_")
        safe_tid = tid.replace(" ", "_").replace(":", "_")
        fasta_name = f"{safe_gene}_{safe_tid}.fa"
        fasta_path = os.path.join(args.out_fasta_dir, fasta_name)

        with open(fasta_path, "w") as out_fh:
            for sp in species_order:
                out_fh.write(f">{sp}\n")
                seq = full_seqs[sp]
                for i in range(0, len(seq), 80):
                    out_fh.write(seq[i:i+80] + "\n")

        written += 1
        if written % 100 == 0:
            print(f"[INFO] Written {written} transcript FASTAs so far...")

        # Summary row: link back to CDS info from GTF
        cds_info = transcripts[(gene_name, tid, chrom, strand)]
        cds_starts = [s for (s, e) in cds_info["cds_intervals"]]
        cds_ends = [e for (s, e) in cds_info["cds_intervals"]]

        row = [
            gene_name,
            tid,
            chrom,
            strand,
            coord_order,
            str(len(cds_info["cds_intervals"])),
            str(len(segments_sorted)),
            str(min(cds_starts)),
            str(max(cds_ends)),
            fasta_path,
        ]
        tsv_lines.append("\t".join(row))

    with open(args.out_tsv, "w") as tsv_fh:
        tsv_fh.write("\n".join(tsv_lines))

    print(f"[INFO] Done. Wrote {written} transcript FASTAs.")
    print(f"[INFO] Summary TSV: {args.out_tsv}")

def main():
    args = parse_args()
    aggregate_blocks(args)

if __name__ == "__main__":
    main()
