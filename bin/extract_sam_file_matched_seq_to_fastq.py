import pysam
import sys
from Bio.Seq import Seq

def remove_soft_clipped(alignment):
    cigar_tuples = alignment.cigartuples
    original_read_sequence = get_forward_sequence(alignment)
    #The alignment is returned as a list of tuples of (operation, length). BAM_CSOFT_CLIP = 4
    if cigar_tuples[0][0] == 4:
        #If aligned fragment starts with softclip start coordinate equals to length of soft-clip
        start = cigar_tuples[0][1]
    else:
        start = 0
    if cigar_tuples[-1][0] == 4:
        #If aligned fragment ends with softclip end coordinate equals to length of soft-clip
        end = -cigar_tuples[-1][1]
    else:
        end = None
    return alignment.seq[start:end], alignment.qual[start:end]


def convert_to_fastq(input_sam):
    samfile = pysam.AlignmentFile(input_sam, "r")
    for read in samfile:
        if read.is_unmapped or len(read.seq) < 20:
#  or read.is_secondary or read.is_supplementary
            continue
        seq, qual = remove_soft_clipped(read)
        sys.stdout.write("@%s\n%s\n+\n%s\n" % (read.qname, seq, qual))
    samfile.close()
#
convert_to_fastq(sys.stdin)