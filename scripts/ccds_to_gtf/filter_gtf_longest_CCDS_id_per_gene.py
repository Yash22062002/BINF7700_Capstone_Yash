#!/usr/bin/env python3

import argparse
import re
import gzip
from collections import defaultdict

def smart_open(filename, mode="rt"):
    """
    Open plain text or gzipped files transparently.
    mode: "rt" for reading text, "wt" for writing text.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Keep only the longest CDS transcript_id per gene_id in a GTF file."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input GTF file (.gtf or .gtf.gz)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output filtered GTF file"
    )
    parser.add_argument(
        "-s", "--summary",
        required=False,
        help="Optional summary text file (default: <output>.summary.txt)"
    )
    return parser.parse_args()

def main():
    args = parse_args()
    input_gtf = args.input
    output_gtf = args.output
    summary_file = args.summary if args.summary else output_gtf + ".summary.txt"

    # Regex to extract gene_id and transcript_id
    gene_re = re.compile(r'gene_id "([^"]+)"')
    tx_re   = re.compile(r'transcript_id "([^"]+)"')

    # First pass: compute total CDS length per transcript per gene
    cds_lengths = defaultdict(lambda: defaultdict(int))  # gene_id -> tx_id -> length

    print(f"[PASS 1] Scanning CDS to compute transcript lengths...")
    with smart_open(input_gtf, "rt") as fin:
        for line in fin:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            feature_type = cols[2]
            if feature_type != "CDS":
                # Only CDS lines define CDS length
                continue

            attr_field = cols[8]
            gene_match = gene_re.search(attr_field)
            tx_match   = tx_re.search(attr_field)

            if not gene_match or not tx_match:
                continue

            gene_id = gene_match.group(1)
            tx_id   = tx_match.group(1)

            try:
                start = int(cols[3])
                end   = int(cols[4])
            except ValueError:
                continue

            length = end - start + 1
            if length < 0:
                continue

            cds_lengths[gene_id][tx_id] += length

    # Decide which transcript to keep per gene: longest CDS
    chosen_tx_for_gene = {}
    for gene_id, tx_dict in cds_lengths.items():
        # Choose max CDS length; if tie, choose lexicographically smallest tx_id
        best_tx = sorted(tx_dict.items(), key=lambda x: (-x[1], x[0]))[0][0]
        chosen_tx_for_gene[gene_id] = best_tx

    print(f"[INFO] Number of genes with CDS data: {len(chosen_tx_for_gene)}")

    # Second pass: write only the chosen transcript per gene
    ignored_tx_for_gene = defaultdict(set)

    total_lines = 0
    kept_lines = 0
    skipped_lines = 0

    print(f"[PASS 2] Filtering GTF using longest CDS isoform per gene...")
    print(f"Input GTF : {input_gtf}")
    print(f"Output GTF: {output_gtf}")
    print(f"Summary   : {summary_file}")

    with smart_open(input_gtf, "rt") as fin, open(output_gtf, "w") as fout:
        for line in fin:
            total_lines += 1

            # Keep header lines as-is
            if line.startswith("#"):
                fout.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                fout.write(line)
                kept_lines += 1
                continue

            attr_field = cols[8]
            gene_match = gene_re.search(attr_field)
            tx_match   = tx_re.search(attr_field)

            if not gene_match or not tx_match:
                # No gene or transcript info -> just keep
                fout.write(line)
                kept_lines += 1
                continue

            gene_id = gene_match.group(1)
            tx_id   = tx_match.group(1)

            # If we have CDS-based choice for this gene, use it.
            chosen_tx = chosen_tx_for_gene.get(gene_id)

            # If gene_id had no CDS lines in pass 1 (unlikely in your CCDS GTF),
            # we don't know which is longest, so keep everything for this gene.
            if chosen_tx is None:
                fout.write(line)
                kept_lines += 1
                continue

            if tx_id == chosen_tx:
                fout.write(line)
                kept_lines += 1
            else:
                ignored_tx_for_gene[gene_id].add(tx_id)
                skipped_lines += 1

    # Summary file
    with open(summary_file, "w") as sf:
        sf.write("Summary: longest CDS transcript per gene\n")
        sf.write("---------------------------------------\n\n")

        for gene_id in sorted(chosen_tx_for_gene.keys()):
            chosen_tx = chosen_tx_for_gene[gene_id]
            chosen_len = cds_lengths[gene_id][chosen_tx]

            # all transcripts for this gene + their lengths
            all_tx = cds_lengths[gene_id]
            ignored_entries = []
            for tx_id, length in all_tx.items():
                if tx_id != chosen_tx:
                    ignored_entries.append(f"{tx_id} ({length} bp)")

            if ignored_entries:
                sf.write(
                    f"Gene {gene_id}: kept transcript {chosen_tx} ({chosen_len} bp); "
                    f"ignored: {', '.join(sorted(ignored_entries))}\n"
                )
            else:
                sf.write(
                    f"Gene {gene_id}: only transcript {chosen_tx} ({chosen_len} bp) observed (no others to ignore).\n"
                )

        sf.write("\nOverall stats:\n")
        sf.write(f"  Total lines read: {total_lines}\n")
        sf.write(f"  Lines kept:       {kept_lines}\n")
        sf.write(f"  Lines skipped:    {skipped_lines}\n")

    # Quick stdout summary
    print("\n[DONE] Finished filtering GTF.")
    print(f"  Total lines read: {total_lines}")
    print(f"  Lines kept:       {kept_lines}")
    print(f"  Lines skipped:    {skipped_lines}")

    print("\nExample gene summary (first few genes):")
    for i, gene_id in enumerate(sorted(chosen_tx_for_gene.keys())):
        if i >= 10:
            break
        chosen_tx = chosen_tx_for_gene[gene_id]
        chosen_len = cds_lengths[gene_id][chosen_tx]
        all_tx = cds_lengths[gene_id]
        ignored_entries = []
        for tx_id, length in all_tx.items():
            if tx_id != chosen_tx:
                ignored_entries.append(f"{tx_id} ({length} bp)")
        if ignored_entries:
            print(
                f"  Gene {gene_id}: kept {chosen_tx} ({chosen_len} bp); "
                f"ignored {', '.join(sorted(ignored_entries))}"
            )
        else:
            print(
                f"  Gene {gene_id}: only {chosen_tx} ({chosen_len} bp) (no other transcripts)"
            )

if __name__ == "__main__":
    main()
