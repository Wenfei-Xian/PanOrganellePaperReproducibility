import sys
import pysam
from collections import defaultdict

def is_primary_alignment(flag):
    return (flag & 0x100) == 0 and (flag & 0x800) == 0

def calculate_coverage(reference_length, covered_positions):
    return (len(covered_positions) / reference_length) * 100

def filter_alignments(input_stream, output_stream, min_coverage=90):
    samfile = pysam.AlignmentFile(input_stream, "r")
    filtered_samfile = pysam.AlignmentFile(output_stream, "w", header=samfile.header)
    
    reference_lengths = {seq: samfile.get_reference_length(seq) for seq in samfile.references}
    coverage_positions = defaultdict(set)
    primary_alignments = []

    for read in samfile.fetch(until_eof=True):
        if is_primary_alignment(read.flag):
            if read.reference_start is not None and read.reference_end is not None:
                primary_alignments.append(read)
                for pos in range(read.reference_start, read.reference_end):
                    coverage_positions[read.reference_name].add(pos)

    for read in primary_alignments:
        ref_name = read.reference_name
        ref_length = reference_lengths[ref_name]
        covered_positions = coverage_positions[ref_name]
        coverage = calculate_coverage(ref_length, covered_positions)

        if coverage >= min_coverage:
            filtered_samfile.write(read)

    samfile.close()
    filtered_samfile.close()

if __name__ == "__main__":
    filter_alignments(sys.stdin, sys.stdout)
