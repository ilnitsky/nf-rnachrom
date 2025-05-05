import argparse
from pathlib import Path
import sys

def process_batki_files(pos_file, neg_file, output_file):
    """Process batki files and write combined results."""
    with open(output_file, 'w') as out:
        # Write header
        out.write("gene_name\tchr\tstart\tend\tgene_type\tfrom_source\tN_counts\n")
        
        # Process both files
        for file_path in [pos_file, neg_file]:
            with open(file_path) as f:
                for line_num, line in enumerate(f, 1):
                    # Remove outer brackets and split by "], ["
                    line = line.strip()[2:-2]  # Remove [[ from start and ]] from end
                    gene_lists = line.split("], [")
                    
                    for gene_info_str in gene_lists:
                        # Split by comma and clean up values
                        gene_info = [item.strip().strip('"') for item in gene_info_str.split(',')]
                        
                        gene, chr_, start, end, rel_start, rel_end, gene_type, from_source, coverage = gene_info
                        
                        # Convert numeric values
                        start = int(start)
                        end = int(end)
                        rel_start = int(rel_start)
                        rel_end = int(rel_end)
                        coverage = float(coverage)
                        
                        # Calculate new value
                        new_value = int(coverage * (rel_end - rel_start + 1))
                        # Write as tab-separated values
                        out.write(f"{gene}\t{chr_}\t{start}\t{end}\t{gene_type}\t{from_source}\t{new_value}\n")
                        

def main():
    parser = argparse.ArgumentParser(description='Combine batki files and calculate new values')
    parser.add_argument('--output-dir', type=str, required=True, help='Output directory')
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    pos_file = output_dir / 'batki.pos_strand.bed'
    neg_file = output_dir / 'batki.neg_strand.bed'
    output_file = output_dir / 'counts.tsv'
    
    try:
        process_batki_files(pos_file, neg_file, output_file)
    except Exception as e:
        print(f"Fatal error: {e}", file=sys.stderr)
        sys.exit(1)

main() 