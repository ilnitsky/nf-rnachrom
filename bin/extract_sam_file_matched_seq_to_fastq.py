import pysam
import sys
from Bio.Seq import Seq

def remove_soft_clipped(alignment):
    cigar_tuples = alignment.cigartuples
    seq = alignment.query_sequence
    qual = alignment.qual
    
    # Initialize start and end indices for slicing
    start, end = 0, len(seq)
    
    if cigar_tuples[0][0] == 4:  # Soft clipping at the start
        start = cigar_tuples[0][1]
    if cigar_tuples[-1][0] == 4:  # Soft clipping at the end
        end -= cigar_tuples[-1][1]
    
    # Adjust for reverse complemented reads
    if alignment.is_reverse:
        seq = seq[::-1]  # Reverse the sequence
        qual = qual[::-1]  # Reverse the quality scores
        # Adjust start and end accordingly
        start, end = len(seq) - end, len(seq) - start

                # Create a Seq object and get the reverse complement
        seq_obj = Seq(seq)
        seq = str(seq_obj.complement())
    
    # Return the sliced sequence and quality scores
    return seq[start:end], qual[start:end]

# def convert_to_fastq(input_sam, output_fastq):
#     samfile = pysam.AlignmentFile(input_sam, "r")
#     with open(output_fastq, 'w') as fastq:
#         for read in samfile:
#             if read.is_unmapped or len(read.query_sequence) < 20:
#                 continue
#             seq, qual = remove_soft_clipped(read)
#             fastq.write("@%s\n%s\n+\n%s\n" % (read.query_name, seq, qual))
#     samfile.close()


def convert_to_fastq(input_sam):
    samfile = pysam.AlignmentFile(input_sam, "r")
    for read in samfile:
        if read.is_unmapped or len(read.query_sequence) < 20:
#  or read.is_secondary or read.is_supplementary
            continue
        seq, qual = remove_soft_clipped(read)
        sys.stdout.write("@%s\n%s\n+\n%s\n" % (read.query_name, seq, qual))
    samfile.close()
#
convert_to_fastq(sys.stdin)