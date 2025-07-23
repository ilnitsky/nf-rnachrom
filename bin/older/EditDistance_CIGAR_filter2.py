#!/usr/bin/env python
import pandas as pd
from matplotlib import pyplot as plt
from collections import Counter
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description='Filter and process contact data based on edit distance and CIGAR string characteristics.')
    parser.add_argument('edit_dist_type', type=str, help='Type of edit distance calculation')
    parser.add_argument('r1_final_threshold', type=int, help='Final threshold for read 1')
    parser.add_argument('r2_final_threshold', type=int, help='Final threshold for read 2')
    parser.add_argument('r1_mapq_threshold', type=int, help='Mapping quality threshold for read 1')
    parser.add_argument('r2_mapq_threshold', type=int, help='Mapping quality threshold for read 2')
    parser.add_argument('distance_threshold', type=int, help='Distance threshold between reads')
    parser.add_argument('assembly_of_ucaRNAs', type=str, help='Assemble ucaRNAs (yes/no)')
    parser.add_argument('mode', type=str, help='Processing mode (basic/explorer)')
    parser.add_argument('experiment_type', type=str, help='Type of experiment')
    parser.add_argument('input_file_name', type=str, help='Name of the input file')
    parser.add_argument('input_file_path', type=str, help='Path to the input file')
    parser.add_argument('output_file_path', type=str, help='Path to the output files')
    return parser.parse_args()


def header_parser(experiment_type):
    """Parses the header information based on the experiment type."""
    srr_id = 'read_id'

    if experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI']:
        return parse_experiment(srr_id, 'ATA', ['rna', 'dna'])

    elif experiment_type in ['OTA_SE', 'OTA_PE']:
        return parse_experiment(srr_id, experiment_type, ['dna1', 'dna2'])

    return parse_experiment(srr_id, experiment_type, ['rna1', 'rna2'])


def parse_experiment(srr_id, experiment_type, field_prefixes):
    """Generates a tuple of header fields for a generic experiment."""
    pairtype = f"{experiment_type}_pairtype"
    fields = [srr_id, pairtype]

    for prefix in field_prefixes:
        fields.extend([
            f'{prefix}_chr', f'{prefix}_start', f'{prefix}_end', f'{prefix}_strand', 
            f'{prefix}_cigar', f'{prefix}_secondary_alignments', f'{prefix}_other_tags', 
            f'{prefix}_NM', f'{prefix}_mapq', f'{prefix}_cigar_type', 
            f'{prefix}_N_softClipp_bp', f'{prefix}_softClipp_type'
        ])
        
    return tuple(fields)

def calculate_edit_distance(nm, softClipp_bp, edit_dist_type):
    """Calculates the final edit distance based on NM field and soft clipping."""
    if edit_dist_type == 'NM + N_softClipp_bp':
        return int(nm) + sum(softClipp_bp)
    return int(nm)


def check_distance_and_chromosomes(r1_chr, r1_start, r1_end, r2_chr, r2_start, r2_end, distance_threshold):
    """Checks if the reads are on the same chromosome and within the specified distance."""
    distance = max(int(r1_start) - int(r2_end) - 1, int(r2_start) - int(r1_end) - 1)
    return (distance <= distance_threshold) and (r1_chr == r2_chr)


def process_contacts(raw_contacts, filtered_contacts, filtered_out_contacts, id_reads, edit_dist_type,
                     r1_final_threshold, r2_final_threshold, r1_mapq_threshold, r2_mapq_threshold, distance_threshold,
                     mode, experiment_type):
    """Processes each contact line for filtering."""
    header_dict = {}
    count = 0
    
    for line in raw_contacts:
        count += 1
        line = line.strip()
        
        if count == 1:
            header_info = header_parser(experiment_type)
            # Update header dict and write output headers
            header_dict.update({item: i for i, item in enumerate(line.split("\t"))})
            # Define output headers based on mode
            output_headers = define_output_headers(mode, header_info)
            filtered_contacts.write(output_headers + "\n")
            filtered_out_contacts.write(output_headers.replace("filtered", "filtered_out") + "\n")
            continue
            
        contact = line.split('\t')
        r1_info = extract_read_info(contact, header_dict, "1")
        r2_info = extract_read_info(contact, header_dict, "2") if r1_info['pairtype'] == 'UU' else None
        
        if r2_info and validate_pairs(r1_info, r2_info, edit_dist_type, distance_threshold, mode):
            output = build_output_line(r1_info, r2_info, mode, header_dict)
            filtered_contacts.write(output + "\n")
            if id_reads:
                id_reads.write(f"{r1_info['SRR_ID']}\t{r1_info['pairtype']}\n")
        else:
            output = build_output_line(r1_info, r2_info, mode, header_dict, filtered=True)
            filtered_out_contacts.write(output + "\n")


def define_output_headers(mode, header_info):
    """Defines output headers based on mode and experiment type."""
    if mode != "explorer":
        return "\t".join(header_info[:-3])  # Exclude last three headers for basic mode
    else:
        return "\t".join(header_info + ['extra_field1', 'extra_field2'])  # Add additional fields for explorer


def extract_read_info(contact, header_dict, read_suffix):
    """Extracts and returns the information for a specific read."""
    return { 
        'SRR_ID': contact[header_dict['read_id']],
        'pairtype': contact[header_dict[f'pairtype']],
        'chr': contact[header_dict[f'rna{read_suffix}_chr']],
        'start': int(contact[header_dict[f'rna{read_suffix}_start']]),
        'end': int(contact[header_dict[f'rna{read_suffix}_end']]),
        'strand': contact[header_dict[f'rna{read_suffix}_strand']],
        'cigar': contact[header_dict[f'rna{read_suffix}_cigar']],
        'NM': int(contact[header_dict[f'rna{read_suffix}_NM']]),
        'mapQ': int(contact[header_dict[f'rna{read_suffix}_mapQ']]),
        'secondary_alignments': contact[header_dict[f'rna{read_suffix}_secondary_alignments']],
        'other_tags': contact[header_dict[f'rna{read_suffix}_other_tags']]
    }


def validate_pairs(r1_info, r2_info, edit_dist_type, distance_threshold, mode):
    """Validates if the pairs meet certain criteria for filtering."""
    # Calculate the edit distances
    r1_final_edit_dist = calculate_edit_distance(r1_info['NM'], [], edit_dist_type)
    r2_final_edit_dist = calculate_edit_distance(r2_info['NM'], [], edit_dist_type)

    if r1_info['pairtype'] == 'UU':
        return (r1_final_edit_dist <= r1_final_threshold and
                r2_final_edit_dist <= r2_final_threshold and
                check_distance_and_chromosomes(r1_info['chr'], r1_info['start'], r1_info['end'], 
                                                r2_info['chr'], r2_info['start'], r2_info['end'], 
                                                distance_threshold))

    return r1_final_edit_dist <= r1_final_threshold


def build_output_line(r1_info, r2_info, mode, header_dict, filtered=False):
    """Builds the output line based on read information. Can be for filtered or filtered_out."""
    line_parts = [
        r1_info['SRR_ID'], r1_info['pairtype'],
        r1_info['chr'], str(r1_info['start']), str(r1_info['end']), r1_info['strand'],
        r1_info['cigar'], str(r1_info['NM']), str(r1_info['mapQ']),
    ]

    if r2_info:
        line_parts += [
            r2_info['chr'], str(r2_info['start']), str(r2_info['end']), r2_info['strand'],
            r2_info['cigar'], str(r2_info['NM']), str(r2_info['mapQ']),
        ]
    else:
        line_parts += ["*"] * 7  # Placeholder if r2_info is not available

    if mode == "explorer":
        extra_fields = [r1_info['cigar_type'], str(sum(r1_info['N_softClipp_bp'])), r1_info['softClipp_type']]
        line_parts += extra_fields

    return "\t".join(line_parts)


def edit_distance_and_cigar_filter(args):
    """Main function to filter based on edit distance and CIGAR string characteristics."""
    with open(f"{args.input_file_path}{args.input_file_name}", 'r') as raw_contacts, \
            open(f"{args.output_file_path}filtered_{args.input_file_name}", 'w') as filtered_contacts, \
            open(f"{args.output_file_path}filtered_out_{args.input_file_name}", 'w') as filtered_out_contacts:
        
        id_reads = None
        if args.assembly_of_ucaRNAs == "yes":
            id_reads = open(f"{args.output_file_path}id_reads_for_ucaRNAs_{args.input_file_name}", 'w')

        process_contacts(raw_contacts, filtered_contacts, filtered_out_contacts, id_reads, args.edit_dist_type,
                         args.r1_final_threshold, args.r2_final_threshold, args.r1_mapq_threshold, args.r2_mapq_threshold,
                         args.distance_threshold, args.mode, args.experiment_type)

        if id_reads:
            id_reads.close()


if __name__ == "__main__":
    args = parse_arguments()
    edit_distance_and_cigar_filter(args)