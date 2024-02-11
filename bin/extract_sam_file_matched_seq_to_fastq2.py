import pysam
import sys


def remove_soft_clipped(alignment):
    cigar_tuples = alignment.cigartuples
    original_read_sequence = alignment.get_forward_sequence()
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
        end = 99999

    if not alignment.is_reverse:
        #Cut out original read sequence based on orientation
        return original_read_sequence[start:end], original_read_sequence[start:end]
    else:
        return original_read_sequence[-end:-start], original_read_sequence[-end:-start]


    


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