#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 bridge_len assembled assembled_positions min_seq_len prefix not_found f_br r_br"
    exit 1
fi

bridge_len=$1
assembled=$2
assembled_positions=$3
min_seq_len=$4
prefix=$5
not_found=$6 
f_br=$7 
r_br=$8 

passed_reads_counter=0
reads_for=0
reads_rev=0

# $1  bridge forward or reverse
# $2  bridge coordinate/position
# $3  read ID 
# $4  sequence 
# $5  quality 

bioawk -c fastx '{print}' ${assembled} | cut -f1-3 | \
paste <(cat ${assembled_positions} | cut -f1,3) - | \
bioawk -v bridge_len=${bridge_len} -v min_seq_len=${min_seq_len} -v prefix=${prefix} -v not_found="$not_found" -v f_br="$f_br" -v r_br="$r_br" 'BEGIN{OFS="\n";} 
{   
  if  ( $1==f_br ) 
  {
    rna_seq=substr($4, 1, $2-1);
    rna_qual=substr($5, 1, $2-1);
    dna_seq=substr($4, $2+bridge_len, 999);
    dna_qual=substr($5, $2+bridge_len, 999);

    if (  (length(rna_seq)>=min_seq_len) &&  (length(dna_seq)>=min_seq_len) ) 
    {
      print "@"$3" length="length(rna_seq), rna_seq, "+", rna_qual >> (prefix"_1.fastq");
      print "@"$3" length="length(dna_seq), dna_seq, "+", dna_qual >> (prefix"_2.fastq");
      passed_reads_counter++;
      reads_for++;
    }
  } 

  if  ( $1==r_br ) 
  {
    rna_seq=substr($4, $2+bridge_len, 999);
    rna_qual=substr($5, $2+bridge_len, 999);
    dna_seq=substr($4, 1, $2-1);
    dna_qual=substr($5, 1, $2-1);

    if (  (length(rna_seq)>=min_seq_len) &&  (length(dna_seq)>=min_seq_len) ) 
    {
      print "@"$3" length="length(rna_seq), revcomp(rna_seq), "+", revcomp(rna_qual) >> (prefix"_1.fastq");
      print "@"$3" length="length(dna_seq), revcomp(dna_seq), "+", revcomp(dna_qual) >> (prefix"_2.fastq");
      passed_reads_counter++;
      reads_rev++;
    }
  }
}' 

echo "Passed Reads Counter: $passed_reads_counter" > $prefix.stats.txt
echo "Reads Forward: $reads_for" >> $prefix.stats.txt
echo "Reads Reverse: $reads_rev" >> $prefix.stats.txt