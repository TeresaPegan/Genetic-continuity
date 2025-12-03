#!/usr/bin/env python3

"""
Human note: I had ChatGPT make this script, starting below:

Convert a per-chromosome MAF alignment into an ancestral states file.

For each reference position with at least one informative outgroup base,
we compute a consensus ancestral base from 3 outgroups and a support value
(1, 2, or 3) = number of outgroups supporting that consensus base.

Output format (one-based coordinates, tab-separated):

    <pos>    <anc_base>    <support>

This matches what the Anchor Continuity PART1 script expects
(aside from the fact that PART1 only uses lines with support == '3').

Usage:
    python maf_to_ancestral_states.py \
        <maf_file> \
        <ref_genome_name> \
        <chrom_name> \
        <outgroup1_name> \
        <outgroup2_name> \
        <outgroup3_name> \
        <output_txt>

Example:
    python maf_to_ancestral_states.py \
        chr1.maf \
        Human \
        chr1 \
        Chimp \
        Gorilla \
        Orangutan \
        chr1_ancestral.txt

Notes / assumptions:
- The MAF file should contain an alignment with:
    - one row per block for the reference genome (ref_genome_name.chrom_name),
    - zero or one row per block for each of the 3 outgroups.
- We ignore blocks where the reference row is not on <chrom_name>.
- Coordinates are derived from the reference row only (MAF is 0-based,
  we output 1-based).
- We ignore reference bases that are not A/C/G/T.
- Outgroup bases that are gap ('-') or not A/C/G/T (e.g. N) are treated
  as missing.
- If there is a tie in the outgroup base counts at a site, we skip that site.
"""

import sys
from collections import Counter

VALID_BASES = set("ACGT")


def parse_src(src):
    """
    Parse the 'src' field from a MAF 's' line.

    Typical format: "<genome>.<chromosome>" or just "<genome>".

    Returns:
        genome_name, chrom_name (chrom_name may be None).
    """
    parts = src.split(".", 1)
    genome = parts[0]
    chrom = parts[1] if len(parts) > 1 else None
    return genome, chrom


def process_block(ref_seq, ref_start, out_seqs, out_names, out_handle):
    """
    Given a completed alignment block, compute ancestral states and write them.

    Args:
        ref_seq  : reference sequence string (aligned, including gaps) in uppercase.
        ref_start: 0-based start position on the reference chromosome.
        out_seqs : dict mapping outgroup name -> sequence string (aligned).
                   Missing outgroups may have empty strings.
        out_names: list of outgroup names in desired order.
        out_handle: open file handle for writing output.
    """
    if ref_seq is None:
        return

    block_len = len(ref_seq)
    if block_len == 0:
        return

    # Ensure all outgroup seqs are the same length; if missing, fill with gaps.
    filled_out_seqs = []
    for name in out_names:
        seq = out_seqs.get(name, "")
        if not seq:
            seq = "-" * block_len
        if len(seq) != block_len:
            # Alignment block is inconsistent; skip it.
            return
        filled_out_seqs.append(seq)

    offset = 0  # number of non-gap bases we've seen in the reference so far

    for i in range(block_len):
        r = ref_seq[i]

        if r == "-":
            # Gap in reference: no coordinate, no ancestor call
            continue

        if r not in VALID_BASES:
            # Non-ACGT ref base; still increment coordinate
            offset += 1
            continue

        # Compute 1-based genomic coordinate from ref_start (0-based in MAF)
        pos = ref_start + offset + 1
        offset += 1

        # Collect outgroup bases at this column
        valid_out_bases = []
        for out_seq in filled_out_seqs:
            b = out_seq[i]
            if b in VALID_BASES:
                valid_out_bases.append(b)

        if not valid_out_bases:
            # No informative outgroups at this site
            continue

        counts = Counter(valid_out_bases)
        most_common = counts.most_common()

        # If there's a tie for the top count between different bases, skip site
        if len(most_common) > 1 and most_common[0][1] == most_common[1][1]:
            continue

        anc_base = most_common[0][0]
        support = most_common[0][1]  # 1, 2, or 3 (max number of outgroups)

        # Write: <pos>\t<anc_base>\t<support>\n
        out_handle.write(f"{pos}\t{anc_base}\t{support}\n")


def main():
    if len(sys.argv) != 8:
        sys.stderr.write(
            "Usage:\n"
            "  python maf_to_ancestral_states.py "
            "<maf_file> <ref_genome> <chrom_name> "
            "<outgroup1> <outgroup2> <outgroup3> <output_txt>\n"
        )
        sys.exit(1)

    maf_file = sys.argv[1]
    ref_genome = sys.argv[2]
    chrom_name = sys.argv[3]
    out1 = sys.argv[4]
    out2 = sys.argv[5]
    out3 = sys.argv[6]
    out_file = sys.argv[7]

    out_names = [out1, out2, out3]

    ref_seq = None
    ref_start = None
    out_seqs = {name: "" for name in out_names}
    in_block = False

    with open(maf_file, "r") as mf, open(out_file, "w") as outf:
        for line in mf:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if line.startswith("a"):
                # Start of a new block: process previous one
                if in_block and ref_seq is not None:
                    process_block(ref_seq, ref_start, out_seqs, out_names, outf)

                # Reset for new block
                in_block = True
                ref_seq = None
                ref_start = None
                out_seqs = {name: "" for name in out_names}
                continue

            if not in_block:
                continue

            if line.startswith("s "):
                parts = line.split()
                # s src start size strand srcSize seq
                if len(parts) < 7:
                    continue
                _, src, start, size, strand, src_size, seq = parts[:7]
                genome, chrom = parse_src(src)
                start = int(start)
                seq = seq.upper()

                if genome == ref_genome:
                    # Only use rows on the specified chromosome
                    if chrom_name is not None and chrom is not None and chrom != chrom_name:
                        # This block's reference is on some other chromosome; ignore it
                        continue
                    ref_seq = seq
                    ref_start = start
                elif genome in out_seqs:
                    # Take the sequence for this outgroup
                    out_seqs[genome] = seq

        # Process the last block if needed
        if in_block and ref_seq is not None:
            process_block(ref_seq, ref_start, out_seqs, out_names, outf)


if __name__ == "__main__":
    main()
