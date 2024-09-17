import argparse
import pysam

def main(infile, outfile, chim_segment_min):
    bamfile = pysam.AlignmentFile(infile, "rb")
    processed_bam = pysam.AlignmentFile(outfile, "wb", template=bamfile)
    chim_segment_min = int(chim_segment_min)

    for read in bamfile.fetch():
        if read.has_tag("SA"):
            if read.is_reverse:
                # For reverse complemented reads, the 5' end is towards the end of the alignment
                if read.cigartuples[-1][0] in (4, 5) and read.cigartuples[-1][1] > chim_segment_min: # 4=S (soft clipping), 5=H (hard clipping) Check if soft-clipping (S) is longer then 13 nts 
                    read.is_supplementary = True
                    read.is_secondary = True
                    processed_bam.write(read)
                    print(f"Read {read.query_name} {read.is_reverse} {read.reference_name} {read.cigarstring} Start: {read.reference_start},End: {read.reference_end}, is 3\' end.")
                else:
                    # Mark the 5' alignment as primary
                    read.is_supplementary = False
                    read.is_secondary = False   # Ensure the read is marked as not 
                    processed_bam.write(read)
                    print(f"Read {read.query_name} {read.is_reverse} {read.reference_name} {read.cigarstring} Start: {read.reference_start},End: {read.reference_end}, is 5\' end ")
            else:
                # For forward reads, the 5' end is towards the beginning of the alignment
                if read.cigartuples[0][0] in (4, 5) and read.cigartuples[0][1] > chim_segment_min:
                    read.is_supplementary = True
                    read.is_secondary = True
                    processed_bam.write(read)
                    print(f"Read {read.query_name} {read.is_reverse} {read.reference_name} {read.cigarstring} Start: {read.reference_start},End: {read.reference_end}, is 3\' end")
                else:
                    read.is_supplementary = False
                    read.is_secondary = False
                    processed_bam.write(read)
                    print(f"Read {read.query_name} {read.is_reverse} {read.reference_name} {read.cigarstring} Start: {read.reference_start},End: {read.reference_end}, is 5\' end ")

    bamfile.close()
    processed_bam.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM files.")
    parser.add_argument("infile", type=str, help="Input BAM file path")
    parser.add_argument("outfile", type=str, help="Output BAM file path")
    parser.add_argument("chim_segment_min", type=int, help="Minimum length of chimeric segment")

    args = parser.parse_args()

    main(args.infile, args.outfile, args.chim_segment_min)