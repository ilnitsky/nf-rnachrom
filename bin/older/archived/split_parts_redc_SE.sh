#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 bridge_len assembled assembled_positions min_seq_len prefix not_found f_br r_br"
    exit 1
fi

# Assign arguments to variables
bridge_len=$1
assembled=$2
assembled_positions=$3
min_seq_len=$4
prefix=$5
not_found=$6 
f_br=$7 
r_br=$8 

# Validate if input files exist
if [ ! -f "$assembled" ]; then
    echo "Error: File $assembled not found."
    exit 1
fi

if [ ! -f "$assembled_positions" ]; then
    echo "Error: File $assembled_positions not found."
    exit 1
fi

# $1  bridge forward or reverse
# $2  bridge coordinate/position
# $3  read ID 
# $4  sequence 
# $5  quality 


passed_reads_counter=0

# Execute the bioawk command with the provided and calculated parameters
bioawk -c fastx '{print}' "$assembled" | cut -f1-3 | \
paste <(cat "$assembled_positions") - | \
bioawk -v bridge_len="$bridge_len" -v min_seq_len="$min_seq_len" -v prefix="$prefix" -v not_found="$not_found" -v f_br="$f_br" -v r_br="$r_br" 'BEGIN{OFS="\n";}
{
    rna_seq=substr($4, $2+bridge_len, 999);
    rna_qual=substr($5, $2+bridge_len, 999);
    dna_seq=substr($4, 1, $2-1);
    dna_qual=substr($5, 1, $2-1);

    # Remove trailing "CCC" from sequences
    if (substr(rna_seq, length(rna_seq)-2) == "CCC") {
        rna_seq = substr(rna_seq, 1, length(rna_seq)-3);
        rna_qual = substr(rna_qual, 1, length(rna_qual)-3);
    }

    if  ( !($1==not_found) && !($1==r_br) && !($1==4) && (length(rna_seq)>=min_seq_len) &&  (length(dna_seq)>=min_seq_len) )  
    {
        print "@"$3" length="length(rna_seq), rna_seq, "+", rna_qual >> (prefix"_1.fastq");
        print "@"$3" length="length(dna_seq), dna_seq"CATG", "+", dna_qual"AAAA" >> (prefix"_2.fastq");
        passed_reads_counter++;
    }  
}'

echo "Passed Reads Counter: $passed_reads_counter" > $prefix.bridge_SE.stats.txt
