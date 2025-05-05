import argparse
from pathlib import Path

def read_genes_file(genes_file):
    """Read genes from bedrc file."""
    genes = set()
    with open(genes_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            chr_, start, end, gene_name, strand, gene_type, from_source = fields
            genes.add((gene_name, chr_, int(start), int(end), gene_type, from_source))
    return genes

def read_batki_genes(batki_file):
    """Read genes that are already in counts.tsv."""
    genes = set()
    with open(batki_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] != 'gene_name':
                gene_name, chr_, start, end, gene_type, from_source, _ = fields
                genes.add((gene_name, chr_, int(start), int(end), gene_type, from_source))
    return genes

def add_missing_genes(genes_file, batki_file):
    """Add missing genes to batki file with zero coverage."""
    all_genes = read_genes_file(genes_file)
    existing_genes = read_batki_genes(batki_file)
    
    # Find missing genes
    missing_genes = all_genes - existing_genes
    
    # Append missing genes to batki file
    with open(batki_file, 'a') as f:
        for gene_name, chr_, start, end, gene_type, from_source in missing_genes:
            # Write with zero coverage
            f.write(f"{gene_name}\t{chr_}\t{start}\t{end}\t{gene_type}\t{from_source}\t0\n")

def main():
    parser = argparse.ArgumentParser(description='Add missing genes to counts.tsv')
    parser.add_argument('--genes-file', type=str, required=True, help='Path to genes.bedrc file')
    parser.add_argument('--output-dir', type=str, required=True, help='Output directory')
    args = parser.parse_args()
    
    genes_file = args.genes_file
    batki_file = Path(args.output_dir) / 'counts.tsv'
    
    add_missing_genes(genes_file, batki_file)

if __name__ == '__main__':
    main() 