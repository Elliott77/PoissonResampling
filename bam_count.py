#!/usr/bin/env python3
"""
bam_count.py
Elliott Ferris, 2019

Count reads for each allele at known strain-distinguishing SNP positions.

This script avoids indel-induced alignment bias by directly inspecting read
sequences at SNP positions, making it suitable for short-read (<50 bp) data
where standard allele-splitting tools may misassign reads near indels.

Requires: pysam

Usage:
    python bam_count.py <BAM> <SNP_file> <output_file>

Arguments:
    BAM          Sorted, indexed BAM file
    SNP_file     Tab-delimited SNP file with columns:
                   [0] chromosome, [1] position, [3] ref_allele (e.g. C57),
                   [4] alt_allele (e.g. CAST), [6] gene/exon ID,
                   [8] SNP quality, [9] near_indel flag ("True"/"False")
    output_file  Path for tab-delimited output (columns: ID, C57_count, CAST_count)
"""

import argparse
import csv
import os
import re
import sys
import time

import pysam


###############################################################################
# CONFIGURATION
###############################################################################

SNP_QUALITY_MINIMUM = 0  # Minimum SNP quality score to include


###############################################################################
# HELPER FUNCTIONS
###############################################################################

def remove_trailing_comma(allele):
    """Strip anything after the first comma in an allele string."""
    match = re.search(r",", allele)
    if match:
        return allele[:match.start()]
    return allele


def matches_allele(sequence, allele, index):
    """Check whether the read sequence matches the given allele at the index."""
    allele_len = len(allele)
    return sequence[index: index + allele_len] == allele.upper()


def excise_insertion(query_seq, cigar_tuples, snp_index):
    """
    Correct the query sequence for upstream insertions so that the SNP index
    into the reference remains valid.

    Returns the (possibly corrected) query sequence.
    """
    match_len = None
    insert_len = None

    for op, length in cigar_tuples:
        if op == 0:  # CIGAR M (match/mismatch)
            match_len = length
        elif op == 1:  # CIGAR I (insertion)
            insert_len = length

        if match_len is not None and insert_len is not None:
            if snp_index > match_len:
                # Insertion is upstream of the SNP — excise it
                return query_seq[:match_len] + query_seq[match_len + insert_len:]
            else:
                return query_seq

    return query_seq


###############################################################################
# MAIN
###############################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Count allele-specific reads at known SNP positions."
    )
    parser.add_argument("bam", help="Sorted, indexed BAM file")
    parser.add_argument("snp_file", help="Tab-delimited SNP reference file")
    parser.add_argument("out_file", help="Output file path")
    parser.add_argument(
        "--min-quality", type=float, default=SNP_QUALITY_MINIMUM,
        help="Minimum SNP quality score (default: %(default)s)"
    )
    args = parser.parse_args()

    # Open input files
    try:
        snp_reader = csv.reader(open(args.snp_file, "r"), delimiter="\t")
    except IOError:
        sys.exit(f"Error: cannot open SNP file '{args.snp_file}'")

    try:
        samfile = pysam.AlignmentFile(args.bam, "rb")
    except IOError:
        sys.exit(f"Error: cannot open BAM file '{args.bam}'")

    # Open output file and write header
    column_name = os.path.basename(args.bam).replace(".bam", "")
    f_out = open(args.out_file, "w")
    f_out.write(f"ChrPos\t{column_name}_C57\t{column_name}_CAST\n")

    # Skip header line
    next(snp_reader)

    # Tracking variables
    start_time = time.time()
    last_chr = None
    last_pos = 0
    snp_count = 0
    total_read_count = 0
    already_written = []
    aw_count = 400

    # Process each SNP
    for line in snp_reader:
        snp_quality = float(line[8])
        near_indel = line[9]

        # Skip low-quality SNPs
        if snp_quality < args.min_quality:
            continue

        # Skip SNPs near indels
        if near_indel == "True":
            continue

        chromosome = line[0]
        pos = int(line[1])
        c57_allele = remove_trailing_comma(line[3])
        cast_allele = remove_trailing_comma(line[4])
        feature_id = line[6]

        snp_count += 1
        count_c57 = 0
        count_cast = 0
        count_neither = 0
        total_depth = 0

        if chromosome != last_chr:
            elapsed = time.time() - start_time
            last_chr = chromosome

        # Cull the already-written list to limit memory usage
        last_snp_depth = aw_count
        if len(already_written) > (2 * last_snp_depth + 20):
            already_written = already_written[-(2 * last_snp_depth + 1):]
        aw_count = 0

        # Iterate over reads overlapping this SNP position
        for read in samfile.fetch(chromosome, pos - 1, pos):
            aw_count += 1

            # Skip reads already counted at a nearby SNP (avoids double-counting)
            if abs(last_pos - pos) < 60 and read.query_name in already_written:
                already_written.append(read.query_name)
                continue

            total_read_count += 1
            ref_positions = read.get_reference_positions()

            if (pos - 1) not in ref_positions:
                continue

            index = ref_positions.index(pos - 1)
            query_seq = excise_insertion(
                read.query_sequence, read.cigartuples or [], index
            )

            total_depth += 1
            alleles_equal_length = (len(c57_allele) == len(cast_allele))

            if alleles_equal_length:
                # Simple SNP — alleles are the same length
                if matches_allele(query_seq, c57_allele, index):
                    count_c57 += 1
                elif matches_allele(query_seq, cast_allele, index):
                    count_cast += 1
                else:
                    count_neither += 1
                already_written.append(read.query_name)

            elif len(c57_allele) > len(cast_allele):
                # C57 has the longer allele (insertion relative to CAST)
                read_end = read.reference_end + 1 if read.reference_end else 0
                if read_end <= pos + len(c57_allele):
                    continue  # Read doesn't span the full insertion
                if matches_allele(query_seq, c57_allele, index):
                    count_c57 += 1
                elif matches_allele(query_seq, cast_allele, index):
                    count_cast += 1
                else:
                    count_neither += 1
                already_written.append(read.query_name)

            else:
                # CAST has the longer allele (insertion relative to C57)
                read_end = read.reference_end + 1 if read.reference_end else 0
                if read_end <= pos + len(cast_allele):
                    count_neither += 1
                    continue
                if matches_allele(query_seq, cast_allele, index):
                    count_cast += 1
                elif matches_allele(query_seq, c57_allele, index):
                    count_c57 += 1
                else:
                    count_neither += 1
                already_written.append(read.query_name)

        last_pos = pos
        f_out.write(f"{feature_id}\t{count_c57}\t{count_cast}\n")

    # Cleanup
    elapsed_time = time.time() - start_time
    print(f"Processed {snp_count} SNPs, {total_read_count} reads")
    print(f"Elapsed time: {elapsed_time / 3600:.2f} hours")

    samfile.close()
    f_out.close()


if __name__ == "__main__":
    main()

