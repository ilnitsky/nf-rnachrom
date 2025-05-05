#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -d <distance> -i <input_path> -o <output_path> -u <uu_file> -m <um_file> -n <N_contacts_min> -f <fdr_threshold> -r <input_path_RNAseq>"
    echo "Example: $0 -d 500000 -i /path/to/input -o /path/to/output -u contacts.voting.UU.bed -m contacts.voting.UM.bed -n 100 -f 0.05 -r /path/to/rnaseq"
    exit 1
}

# Parse command line arguments
while getopts "d:i:o:u:m:n:f:r:" opt; do
    case $opt in
        d) d="$OPTARG";;
        i) input_path="$OPTARG";;
        o) output_path="$OPTARG";;
        u) uu_file="$OPTARG";;
        m) um_file="$OPTARG";;
        n) n_contacts_min="$OPTARG";;
        f) fdr_threshold="$OPTARG";;
        r) input_path_RNAseq="$OPTARG";;
        *) usage;;
    esac
done

# Check if all required parameters are provided
if [ -z "$d" ] || [ -z "$input_path" ] || [ -z "$output_path" ] || [ -z "$uu_file" ] || [ -z "$um_file" ] || [ -z "$n_contacts_min" ] || [ -z "$fdr_threshold" ] || [ -z "$input_path_RNAseq" ]; then
    usage
fi

# Check if input directory and files exist
if [ ! -d "$input_path" ]; then
    echo "Error: Input directory does not exist: $input_path"
    exit 1
fi

if [ ! -f "$input_path/$uu_file" ]; then
    echo "Error: UU contacts file does not exist: $input_path/$uu_file"
    exit 1
fi

if [ ! -f "$input_path/$um_file" ]; then
    echo "Error: UM contacts file does not exist: $input_path/$um_file"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_path"

# Define output filenames
output_file_uu_all="$output_path/counts_contacts_UU_all.tsv"
output_file_uu_dist="$output_path/counts_contacts_UU_filter_dist_${d}.tsv"
output_file_um_dist="$output_path/counts_contacts_UM_filter_dist_${d}.tsv"

# Create temporary files
temp_file_uu_all=$(mktemp)
temp_file_uu_dist=$(mktemp)

# Process UU contacts
awk -F'\t' -v d="$d" '
function abs(x) { return x < 0 ? -x : x }

NR > 1 {
    if (!($7 in gene_info)) {
        gene_info[$7] = $8 "\t" $9
    }
    # Count all contacts
    genes_all[$7]++
    
    # Check distance for the second file
    rna_chr = $3
    rna_start = $4
    dna_chr = $13
    dna_start = $14
    
    if (rna_chr == dna_chr) {  # If chromosomes match
        if (abs(rna_start - dna_start) > d) {  # If distance is greater than d
            genes_dist[$7]++
        }
    } else {  # If chromosomes are different, always count
        genes_dist[$7]++
    }
}

END {
    # Output results for all contacts
    for (gene in genes_all) {
        print gene "\t" gene_info[gene] "\t" genes_all[gene]
    }
    
    print "---SPLIT---"
    
    # Output results for distance-filtered contacts
    for (gene in genes_dist) {
        print gene "\t" gene_info[gene] "\t" genes_dist[gene]
    }
}' "$input_path/$uu_file" | \
awk -v file1="$temp_file_uu_all" -v file2="$temp_file_uu_dist" '
    /---SPLIT---/ { split_found=1; next }
    !split_found { print > file1 }
    split_found { print > file2 }
'

# Process UM contacts
awk -F'\t' -v d="$d" '
function abs(x) { return x < 0 ? -x : x }

NR > 1 {
    if (!($7 in gene_info)) {
        gene_info[$7] = $8 "\t" $9
    }
}

NR > 1 {
    rna_chr = $3
    rna_start = $4
    dna_chr = $13
    dna_start = $14
    sec_aligns = $21
    
    # Check primary contact
    found_match = 0
    if (rna_chr == dna_chr) {
        if (abs(rna_start - dna_start) < d) {
            found_match = 1
        }
    }
    
    # If primary contact doesnt match, check secondary alignments
    if (!found_match && sec_aligns != "*") {
        gsub(/[()]/, "", sec_aligns)
        split(sec_aligns, aligns, ";")
        
        for (i in aligns) {
            split(aligns[i], parts, ",")
            sec_chr = parts[1]
            sec_pos = parts[2]
            
            if (sec_chr == rna_chr && abs(rna_start - sec_pos) < d) {
                found_match = 1
                break
            }
        }
    }
    
    if (!found_match) {
        genes[$7]++
    }
}

END {
    for (gene in genes) {
        print gene "\t" gene_info[gene] "\t" genes[gene]
    }
}' "$input_path/$um_file" > "$output_file_um_dist"

# Create final UU files with headers
echo "gene_name\tgene_type\tfrom_source\tN_counts" > "$output_file_uu_all"
cat "$temp_file_uu_all" >> "$output_file_uu_all"

echo "gene_name\tgene_type\tfrom_source\tN_counts" > "$output_file_uu_dist"
cat "$temp_file_uu_dist" >> "$output_file_uu_dist"

# Add header to UM file
sed -i '1i gene_name\tgene_type\tfrom_source\tN_counts' "$output_file_um_dist"

# Remove temporary files
rm "$temp_file_uu_all" "$temp_file_uu_dist"

# Create combined output file
output_file_combined="$output_path/counts_contacts_UU_UM_filter_dist_${d}.tsv"

# Merge files and sum N_counts for identical gene_names
awk -F'\t' '
    # Skip header for second file
    FNR == 1 { next }
    
    {
        # Key is gene_name, gene_type and from_source
        key = $1 "\t" $2 "\t" $3
        # Sum N_counts
        counts[key] += $4
        # Preserve order of appearance
        if (!(key in seen)) {
            seen[key] = 1
            keys[++n] = key
        }
    }
    
    END {
        # Output header
        print "gene_name\tgene_type\tfrom_source\tN_counts"
        # Output results in order of appearance
        for (i = 1; i <= n; i++) {
            key = keys[i]
            print key "\t" counts[key]
        }
    }
' "$output_file_uu_dist" "$output_file_um_dist" > "$output_file_combined"


rm "$output_file_um_dist"
cp "$input_path/counts.tsv" "$output_path"

# Run Python script for different datasets
python3 RD_chP.py \
    --input_path "$output_path" \
    --output_path "$output_path" \
    --input_path_RNAseq "$input_path_RNAseq" \
    --counts_RNAseq "counts.tsv" \
    --counts_contacts "counts_contacts_UU_all.tsv" \
    --N_contacts_min "$n_contacts_min" \
    --fdr_threshold "$fdr_threshold" \
    --type "UU_all"

python3 RD_chP.py \
    --input_path "$output_path" \
    --output_path "$output_path" \
    --input_path_RNAseq "$input_path_RNAseq" \
    --counts_RNAseq "counts.tsv" \
    --counts_contacts "counts_contacts_UU_filter_dist_${d}.tsv" \
    --N_contacts_min "$n_contacts_min" \
    --fdr_threshold "$fdr_threshold" \
    --type "UU_filter_dist_${d}"

python3 RD_chP.py \
    --input_path "$output_path" \
    --output_path "$output_path" \
    --input_path_RNAseq "$input_path_RNAseq" \
    --counts_RNAseq "counts.tsv" \
    --counts_contacts "counts.tsv" \
    --N_contacts_min "$n_contacts_min" \
    --fdr_threshold "$fdr_threshold" \
    --type "UU_UM_all"

python3 RD_chP.py \
    --input_path "$output_path" \
    --output_path "$output_path" \
    --input_path_RNAseq "$input_path_RNAseq" \
    --counts_RNAseq "counts.tsv" \
    --counts_contacts "counts_contacts_UU_UM_filter_dist_${d}.tsv" \
    --N_contacts_min "$n_contacts_min" \
    --fdr_threshold "$fdr_threshold" \
    --type "UU_UM_filter_dist_${d}"

# /usr/bin/time -v sh count_contacts_all.sh \
#                        -d 500000 \
#                        -i /home/snap/projects/lncRNA_app/voting/output_SRR17331253_UU_UM \
#                        -o /home/snap/projects/lncRNA_app/voting/output_SRR17331253_UU_UM/chP \
#                        -u contacts.voting.UU.bed \
#                        -m contacts.voting.UM.bed \
#                        -n 100 \
#                        -f 0.05 \
#                        -r /home/snap/Downloads