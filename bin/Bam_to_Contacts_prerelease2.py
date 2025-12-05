#!/usr/bin/env python
import pysam
import argparse
import sys

def process_bam(filename):
    return pysam.AlignmentFile(filename, "rb").fetch(until_eof=True)

def extract_nm_tag(read):
    return str(read.get_tag('NM')) if read and read.has_tag('NM') else '*'

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
        str(read.reference_start),
        str(read.reference_end),
        '-' if read.is_reverse else '+',
        read.cigarstring,
        extract_nm_tag(read),
        str(read.mapping_quality)
    ]

def format_secondary_alignment(read):
    strand = '-' if read.is_reverse else '+'
    nm = read.get_tag('NM') if read.has_tag('NM') else '0'
    return f"({read.reference_name},{read.reference_start},{strand}{read.cigarstring},{nm})"

def get_mapping_status(alignments, secondary_alignments):
    if not alignments or alignments[0] == '*':
        return 'N'
    return 'M' if secondary_alignments else 'U'

def process_reads(r1_iter, r2_iter, other_tags, rna_mode='STAR', dna_mode='STAR'):
    def get_base_name(query_name):
        return query_name.split('.')[-1] if query_name else None

    r1, r2 = next(r1_iter, None), next(r2_iter, None)
    while r1 is not None or r2 is not None:
        query_name = r1.query_name if r1 else r2.query_name
        base_name = get_base_name(query_name)
        r1_data = {'alignments': [], 'secondary_alignments': [], 'other_tags': '*'}
        r2_data = {'alignments': [], 'secondary_alignments': [], 'other_tags': '*'}

        # Process R1 reads (RNA)
        while r1 and get_base_name(r1.query_name) == base_name:
            if rna_mode == 'BWA':
                r1_data['alignments'] = format_alignment(r1)
                r1_data['other_tags'] = extract_other_tags(r1, other_tags)
                if r1.has_tag('XA'):
                    r1_data['secondary_alignments'].append(r1.get_tag('XA'))
            else:
                # HISAT2/STAR handling for RNA
                if r1.is_secondary or r1.is_supplementary:
                    r1_data['secondary_alignments'].append(format_secondary_alignment(r1))
                else:
                    r1_data['alignments'] = format_alignment(r1)
                    r1_data['other_tags'] = extract_other_tags(r1, other_tags)
            r1 = next(r1_iter, None)

        # Process R2 reads (DNA)
        while r2 and get_base_name(r2.query_name) == base_name:
            if dna_mode == 'BWA':
                r2_data['alignments'] = format_alignment(r2)
                r2_data['other_tags'] = extract_other_tags(r2, other_tags)
                if r2.has_tag('XA'):
                    r2_data['secondary_alignments'].append(r2.get_tag('XA'))
            else:
                # HISAT2/STAR handling for DNA
                if r2.is_secondary or r2.is_supplementary:
                    r2_data['secondary_alignments'].append(format_secondary_alignment(r2))
                else:
                    r2_data['alignments'] = format_alignment(r2)
                    r2_data['other_tags'] = extract_other_tags(r2, other_tags)
            r2 = next(r2_iter, None)

        yield query_name, r1_data, r2_data

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process BAM files and output alignment information.')
    parser.add_argument('-r1', '--rna_bam', required=True, help='Input RNA BAM file (R1)')
    parser.add_argument('-r2', '--dna_bam', help='Input DNA BAM file (R2)')
    parser.add_argument('-mr', '--rna_mode', choices=['BWA', 'STAR', 'HISAT', 'BOWTIE'], required=True, 
                        help='RNA alignment mode: BWA, STAR, HISAT, or BOWTIE')
    parser.add_argument('-md', '--dna_mode', choices=['BWA', 'STAR', 'HISAT', 'BOWTIE'], 
                        help='DNA alignment mode: BWA, STAR, HISAT, or BOWTIE (defaults to RNA mode if not specified)')
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

def main():
    args = parse_arguments()
    
    # Set default DNA mode to RNA mode if not specified
    if not args.dna_mode:
        args.dna_mode = args.rna_mode
    
    try:
        unique_file = f"{args.prefix}_Unique_RNA.tab.rc"
        other_file = f"{args.prefix}_Other.tab.rc"

        # For RNA-seq SE, use dummy values for r2
        if args.exp_type == 'RNA_SEQ_SE':
            args.dna_bam = args.rna_bam  
            processing_exp_type = 'RNA_SEQ_SE'
        elif args.exp_type == 'RNA_SEQ_PE':
            if not args.dna_bam:
                args.dna_bam = args.rna_bam  # Use the same file for PE RNA-seq
            processing_exp_type = 'RNA_SEQ_PE'
        else:
            processing_exp_type = args.exp_type
            if not args.dna_bam:
                raise ValueError("DNA BAM file (-r2) is required for RNA-DNA experiments")

        with open(unique_file, 'w') as f_unique, open(other_file, 'w') as f_other:
            write_header(f_unique, processing_exp_type)
            write_header(f_other, processing_exp_type)

            for query_name, r1_data, r2_data in process_reads(
                process_bam(args.rna_bam), 
                process_bam(args.dna_bam), 
                args.other_tags, 
                args.rna_mode,
                args.dna_mode
            ):
                r1_status = get_mapping_status(r1_data['alignments'], r1_data['secondary_alignments'])
                r2_status = get_mapping_status(r2_data['alignments'], r2_data['secondary_alignments'])
                pairtype = f"{r1_status}{r2_status}"
                print(r1_data, r2_data, pairtype)
                r1_alignments = r1_data['alignments'] if r1_data['alignments'] else ['*'] * 7
                r2_alignments = r2_data['alignments'] if r2_data['alignments'] else ['*'] * 7

                r1_secondary = ';'.join(r1_data['secondary_alignments']) if r1_data['secondary_alignments'] else '*'
                r2_secondary = ';'.join(r2_data['secondary_alignments']) if r2_data['secondary_alignments'] else '*'

                line = [
                    query_name, pairtype,
                    *r1_alignments,
                    *r2_alignments,
                    r1_secondary, r2_secondary,
                    r1_data['other_tags'], r2_data['other_tags']
                ]

                output_line = '\t'.join(line) + '\n'

                if pairtype in ['UU', 'UM']:
                    f_unique.write(output_line)
                else:
                    f_other.write(output_line)

    except IOError as e:
        print(f"Error processing files: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()