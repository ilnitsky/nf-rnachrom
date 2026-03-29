#!/usr/bin/env python
### last update: 29.03.2026
import pysam
import argparse
import sys

def process_bam(filename):
    return pysam.AlignmentFile(filename, "rb").fetch(until_eof=True)

def extract_nm_tag(read):
    if not read:
        return '*'
    try:
        return read.get_tag('NM')
    except (KeyError, AttributeError):
        return '*'

def extract_other_tags(read, other_tags):
    if not read:
        return '*'
    tags = []
    for tag in other_tags:
        if read.has_tag(tag):
            value = read.get_tag(tag)
            tags.append(f'{tag}:"{value}"')
    return ';'.join(tags) if tags else '*'

def format_alignment(read):
    if read.is_unmapped:
        return ['*'] * 7  # chr, start, end, strand, cigar, NM, mapq
    return [
        read.reference_name,
        str(read.reference_start + 1), #1-based
        str(read.reference_end),       #1-based
        '-' if read.is_reverse else '+',
        read.cigarstring,
        str(extract_nm_tag(read)),
        str(read.mapping_quality)
    ]

def format_secondary_alignment(read):
    strand = '-' if read.is_reverse else '+'
    nm = read.get_tag('NM') if read.has_tag('NM') else '*'
    return f"({read.reference_name},{read.reference_start},{strand}{read.cigarstring},{nm})"

def get_mapping_status(alignments, secondary_alignments):
    # Defensive: treat missing or empty alignments as Not mapped ('N').
    # If the caller uses ['*'] to mark explicit unmapped, handle that too.
    if not alignments:
        return 'N'
    try:
        if alignments[0] == '*':
            return 'N'
    except Exception:
        return 'N'
    return 'M' if secondary_alignments else 'U'

def get_base_name(query_name):
    return query_name.split('.')[-1] if query_name else None

def process_reads(r1_iter, r2_iter, other_tags): #, rna_mode='STAR', dna_mode='STAR'
    # 1. SUPPLEMENTARY READS FILTRATION
    r1_iter_filtered = (r for r in r1_iter if not r.is_supplementary)  # Generator for RNA
    r2_iter_filtered = (r for r in r2_iter if not r.is_supplementary)  # Generator for DNA
    # Local bindings to speed up hot loop
    fmt_secondary = format_secondary_alignment

    # 2. INITIALIZING THE FIRST READS
    r1, r2 = next(r1_iter_filtered, None), next(r2_iter_filtered, None)

    while r1 is not None or r2 is not None:
        query_name = r1.query_name if r1 else r2.query_name
        base_name = get_base_name(query_name)
        r1_data = {'alignments': [], 'secondary_alignments': [], 'other_tags': '*'}
        r2_data = {'alignments': [], 'secondary_alignments': [], 'other_tags': '*'}

        # Process R1 reads (RNA)
        while r1:
            # compare query_name suffix quickly (avoid repeated split/rfind)
            r1_qname = r1.query_name
            if not (r1_qname == base_name or r1_qname.endswith('.' + base_name)):
                break
            # if rna_mode == 'BWA':
            #     if r1.is_unmapped:
            #         r1_data['alignments'] = ['*'] * 7
            #         r1_data['other_tags'] = '*'
            #     else:
            #         # Use existing helpers: format_alignment and extract_other_tags
            #         r1_data['alignments'] = format_alignment(r1)
            #         r1_data['other_tags'] = extract_other_tags(r1, other_tags)
            #         if r1.has_tag('XA'):
            #             r1_data['secondary_alignments'].append(r1.get_tag('XA'))
            # else:
            if r1.is_secondary:
                r1_data['secondary_alignments'].append(fmt_secondary(r1))
            else:
                if r1.is_unmapped:
                    r1_data['alignments'] = ['*'] * 7
                    r1_data['other_tags'] = '*'
                else:
                    r1_data['alignments'] = format_alignment(r1)
                    r1_data['other_tags'] = extract_other_tags(r1, other_tags)
            r1 = next(r1_iter_filtered, None)

        # Process R2 reads (DNA)
        while r2:
            r2_qname = r2.query_name
            if not (r2_qname == base_name or r2_qname.endswith('.' + base_name)):
                break
            # if dna_mode == 'BWA':
            #     if r2.is_unmapped:
            #         r2_data['alignments'] = ['*'] * 7
            #         r2_data['other_tags'] = '*'
            #     else:
            #         # Use existing helpers for BWA primary alignments
            #         r2_data['alignments'] = format_alignment(r2)
            #         r2_data['other_tags'] = extract_other_tags(r2, other_tags)
            #         if r2.has_tag('XA'):
            #             r2_data['secondary_alignments'].append(r2.get_tag('XA'))
            # else:
            if r2.is_secondary:
                r2_data['secondary_alignments'].append(fmt_secondary(r2))
            else:
                if r2.is_unmapped:
                    r2_data['alignments'] = ['*'] * 7
                    r2_data['other_tags'] = '*'
                else:
                    r2_data['alignments'] = format_alignment(r2)
                    r2_data['other_tags'] = extract_other_tags(r2, other_tags)
            r2 = next(r2_iter_filtered, None)

        yield query_name, r1_data, r2_data

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process BAM files and output alignment information.')
    parser.add_argument('-r1', '--rna_bam', required=True, help='Input RNA BAM file (R1)')
    parser.add_argument('-r2', '--dna_bam', help='Input DNA BAM file (R2)')
    # parser.add_argument('-mr', '--rna_mode', choices=['BWA', 'STAR', 'HISAT'], required=True, 
    #                     help='RNA alignment mode: BWA, STAR, or HISAT')
    # parser.add_argument('-md', '--dna_mode', choices=['BWA', 'STAR', 'HISAT'], 
    #                     help='DNA alignment mode: BWA, STAR, or HISAT (defaults to RNA mode if not specified)')
    parser.add_argument('-e', '--exp_type', choices=['OTA_PE', 'OTA_SE', 'ATA', 'RNA_SEQ_PE', 'RNA_SEQ_SE'], required=True,
                        help='Experiment type: OTA_PE, OTA_SE, ATA, RNA_SEQ_PE, or RNA_SEQ_SE')
    parser.add_argument('-t', '--other_tags', nargs='+', default=["NH"],
                        help='List of additional SAM tags to extract (default: NH)')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Output file prefix')
    return parser.parse_args()

def write_header(file, exp_type):
    if exp_type == 'ATA':
        header = ["read_id", "ATA_pairtype", 
                "rna_chr", "rna_start", "rna_end", "rna_strand", "rna_cigar", "rna_NM", "rna_mapq",
                "dna_chr", "dna_start", "dna_end", "dna_strand", "dna_cigar", "dna_NM", "dna_mapq",
                "rna_secondary_alignments", "dna_secondary_alignments",
                "rna_other_tags", "dna_other_tags"]
    elif exp_type == 'OTA_PE':
        header = ["read_id", "OTA_PE_pairtype", 
                "dna1_chr", "dna1_start", "dna1_end", "dna1_strand", "dna1_cigar", "dna1_NM", "dna1_mapq",
                "dna2_chr", "dna2_start", "dna2_end", "dna2_strand", "dna2_cigar", "dna2_NM", "dna2_mapq",
                "dna1_secondary_alignments", "dna2_secondary_alignments",
                "dna1_other_tags", "dna2_other_tags"]
    elif exp_type == 'OTA_SE':          
        header = ["read_id", "OTA_SE_pairtype", 
                "dna1_chr", "dna1_start", "dna1_end", "dna1_strand", "dna1_cigar", "dna1_NM", "dna1_mapq",
                "dna2_chr", "dna2_start", "dna2_end", "dna2_strand", "dna2_cigar", "dna2_NM", "dna2_mapq",
                "dna1_secondary_alignments", "dna2_secondary_alignments",
                "dna1_other_tags", "dna2_other_tags"]
    elif exp_type == 'RNA_SEQ_PE':          
        header = ["read_id", "RNAseq_PE_pairtype", 
                "rna1_chr", "rna1_start", "rna1_end", "rna1_strand", "rna1_cigar", "rna1_NM", "rna1_mapq",
                "rna2_chr", "rna2_start", "rna2_end", "rna2_strand", "rna2_cigar", "rna2_NM", "rna2_mapq",
                "rna1_secondary_alignments", "rna2_secondary_alignments",
                "rna1_other_tags", "rna2_other_tags"]
    elif exp_type == 'RNA_SEQ_SE':          
        header = ["read_id", "RNAseq_SE_pairtype", 
                "rna1_chr", "rna1_start", "rna1_end", "rna1_strand", "rna1_cigar", "rna1_NM", "rna1_mapq",
                "rna2_chr", "rna2_start", "rna2_end", "rna2_strand", "rna2_cigar", "rna2_NM", "rna2_mapq",
                "rna1_secondary_alignments", "rna2_secondary_alignments",
                "rna1_other_tags", "rna2_other_tags"]
              
    file.write('\t'.join(header) + '\n')

def are_bams_queryname_sorted(file1, file2):
    """Checks the sorting of two BAM files by queryname"""
    for filename in [file1, file2]:
        try:
            with pysam.AlignmentFile(filename, "rb") as bam:
                header = bam.header
                if not ('HD' in header and 'SO' in header['HD'] and header['HD']['SO'] == 'queryname'):
                    return False
        except:
            return False
    return True


def count_unique_read_ids(filename):
    """Counts unique read IDs in a BAM file. Assumes the file is sorted by queryname."""
    unique_count = 0
    prev_read_id = None
    try:
        with pysam.AlignmentFile(filename, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                current_read_id = read.query_name
                if prev_read_id is None or current_read_id != prev_read_id:
                    unique_count += 1
                    prev_read_id = current_read_id
        return unique_count
    except Exception as e:
        return -1  # Return -1 if an error occurs.



def main():
    args = parse_arguments()
    
    # Set default DNA mode to RNA mode if not specified
    # if not args.dna_mode:
    #     args.dna_mode = args.rna_mode
    try:
        unique_file = f"{args.prefix}_Unique_RNA.tab.rc"
        other_file = f"{args.prefix}_Other.tab.rc"

        # For RNA-seq SE / OTA_SE we do NOT want to process the same BAM twice.
        # Avoid setting args.dna_bam = args.rna_bam which would cause redundant work.
        if args.exp_type in ['RNA_SEQ_SE','OTA_SE']:
            processing_exp_type = args.exp_type
        elif args.exp_type == 'RNA_SEQ_PE':
            if not args.dna_bam:
                # args.dna_bam = args.rna_bam  # Use the same file for PE RNA-seq
                raise ValueError("BAM file (-r2) is required for PE RNA-seq experiments")
            processing_exp_type = 'RNA_SEQ_PE'
        else:
            processing_exp_type = args.exp_type
            if not args.dna_bam:
                raise ValueError("DNA BAM file (-r2) is required for RNA-DNA experiments")
        
        if processing_exp_type in ['ATA', 'OTA_PE', 'RNA_SEQ_PE']:
            bam_files_sorted = are_bams_queryname_sorted(args.rna_bam, args.dna_bam)
            if not bam_files_sorted:
                raise ValueError("Both BAM files must be sorted by query name for paired-end experiments.")
            else:
                count1 = count_unique_read_ids(args.rna_bam)
                count2 = count_unique_read_ids(args.dna_bam)
                if (count1 != count2) or count1 == -1 or count2 == -1:
                    raise ValueError("Two BAM files have different numbers of unique read identifiers, but should have the same number.")
                
        unique_set = {'UU'} if args.exp_type in ['OTA_PE', 'RNA_SEQ_PE'] else {'UU', 'UM', 'U'}
                
        # use a larger buffering size for output files to reduce syscalls
        with open(unique_file, 'w', buffering=262144) as f_unique, open(other_file, 'w', buffering=262144) as f_other:
            write_header(f_unique, processing_exp_type)
            write_header(f_other, processing_exp_type)

            # Prepare iterators: for SE experiments we avoid opening the DNA BAM
            # and pass an empty iterator for r2 to prevent redundant work.
            r1_iter = process_bam(args.rna_bam)
            if processing_exp_type in ['RNA_SEQ_SE', 'OTA_SE']:
                r2_iter = iter(())
            else:
                r2_iter = process_bam(args.dna_bam)

            # batched writing to reduce Python->syscall overhead
            # increase to 100k lines per flush to further reduce syscalls
            FLUSH_LINES = 100000
            buf_unique = []
            buf_other = []

            for query_name, r1_data, r2_data in process_reads(
                r1_iter,
                r2_iter,
                args.other_tags,
                # args.rna_mode,
                # args.dna_mode
            ):
                r1_status = get_mapping_status(r1_data['alignments'], r1_data['secondary_alignments'])
                if processing_exp_type in ['RNA_SEQ_SE', 'OTA_SE']:
                    r2_status = '' # No r2 data for SE experiments, so status is empty string
                else:
                    r2_status = get_mapping_status(r2_data['alignments'], r2_data['secondary_alignments'])
                pairtype = f"{r1_status}{r2_status}"

                r1_alignments = r1_data['alignments'] if r1_data['alignments'] else ['*'] * 7
                r2_alignments = r2_data['alignments'] if r2_data['alignments'] else ['*'] * 7

                # format_alignment already returns strings, use them directly
                r1_alignments_str = r1_alignments
                r2_alignments_str = r2_alignments

                r1_secondary = ';'.join(r1_data['secondary_alignments']) if r1_data['secondary_alignments'] else '*'
                r2_secondary = ';'.join(r2_data['secondary_alignments']) if r2_data['secondary_alignments'] else '*'

                r1_other = r1_data['other_tags'] if r1_data['other_tags'] else '*'
                r2_other = r2_data['other_tags'] if r2_data['other_tags'] else '*'

                line = [
                    query_name, pairtype,
                    *r1_alignments_str,
                    *r2_alignments_str,
                    r1_secondary, r2_secondary,
                    r1_other, r2_other
                ]

                output_line = '\t'.join(line) + '\n'

                if pairtype in unique_set:
                    buf_unique.append(output_line)
                    if len(buf_unique) >= FLUSH_LINES:
                        f_unique.write(''.join(buf_unique))
                        buf_unique.clear()
                else:
                    buf_other.append(output_line)
                    if len(buf_other) >= FLUSH_LINES:
                        f_other.write(''.join(buf_other))
                        buf_other.clear()

            # flush remaining buffers
            if buf_unique:
                f_unique.write(''.join(buf_unique))
            if buf_other:
                f_other.write(''.join(buf_other))

    except IOError as e:
        print(f"Error processing files: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
    # # Optional line-by-line profiling with line_profiler:
    # # Run as: python Bam_to_Contacts_grisha_fix_old.py --profile --profile-output profile.txt [other args...]
    # if '--profile' in sys.argv:
    #     try:
    #         from line_profiler import LineProfiler
    #     except ImportError:
    #         print("line_profiler is not installed. Install with: pip install line_profiler", file=sys.stderr)
    #         sys.exit(1)
    #     # --profile requires --profile-output <file>
    #     if '--profile-output' not in sys.argv:
    #         print("When using --profile you must provide --profile-output <file>", file=sys.stderr)
    #         sys.exit(1)
    #     # extract profile output path
    #     po_idx = sys.argv.index('--profile-output')
    #     try:
    #         profile_out = sys.argv[po_idx + 1]
    #     except Exception:
    #         print("Missing value for --profile-output <file>", file=sys.stderr)
    #         sys.exit(1)
    #     # remove the flags/values so argparse in main() won't see them
    #     for token in ['--profile', '--profile-output', profile_out]:
    #         while token in sys.argv:
    #             sys.argv.remove(token)

    #     lp = LineProfiler()
    #     # register functions to profile
    #     for fn in (process_bam, process_reads, format_alignment,
    #                format_secondary_alignment, extract_other_tags,
    #                extract_nm_tag, get_mapping_status):
    #         try:
    #             lp.add_function(fn)
    #         except Exception:
    #             pass

    #     # run main under the profiler and write results to file
    #     lp.runcall(main)
    #     try:
    #         # write profile output with buffering to reduce I/O overhead
    #         with open(profile_out, 'w', buffering=65536) as f:
    #             lp.print_stats(stream=f)
    #     except Exception as e:
    #         print(f"Failed to write profile output to {profile_out}: {e}", file=sys.stderr)
    #         sys.exit(1)
    # else:
    #     main()
