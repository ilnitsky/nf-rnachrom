#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 10 ]; then
    echo "Usage: $0 bridge_len assembled_F assembled_positions_F assembled_R assembled_positions_R min_seq_len prefix not_found f_br r_br"
    exit 1
fi


bridge_len=$1
unassembled_F=$2
unassembled_positions_F=$3
unassembled_R=$4
unassembled_positions_R=$5
min_seq_len=$6
prefix=$7
not_found=$8
f_br=$9
r_br=${10}


passed_reads_counter=0
reads_for=0
reads_rev=0

# $1  unassembled_positions_F bridge forward or reverse
# $2  unassembled_positions_F bridge coordinate/position
# $3  read ID (unassembled_F)
# $4  sequence (unassembled_F)
# $5  quality (unassembled_F)
# $6  unassembled_positions_R bridge forward or reverse
# $7  unassembled_positions_R bridge coordinate/position
# $8  read ID (unassembled_R)
# $9  sequence (unassembled_R)
# $10 quality (unassembled_R)

paste \
  <(cat ${unassembled_positions_F} | cut -f1,3) \
  <(bioawk -c fastx '{print}' ${unassembled_F} | cut -f1-3) \
  <(cat ${unassembled_positions_R} | cut -f1,3) \
  <(bioawk -c fastx '{print}' ${unassembled_R} | cut -f1-3) | \
bioawk -v bridge_len=${bridge_len} -v prefix=${prefix} -v min_seq_len=${min_seq_len} -v not_found="$not_found" -v f_br="$f_br" -v r_br="$r_br" 'BEGIN{OFS="\n";} 
{
  if  ( ($1==f_br) && ($6==not_found)  ) 
  {
    rna_seq=substr($4, $2+bridge_len, 999);
    rna_qual=substr($5, $2+bridge_len, 999);
    dna_seq=substr($4, 1, $2-1);
    dna_qual=substr($5, 1, $2-1);

    if (  (length(rna_seq)>=min_seq_len) &&  (length(dna_seq)>=min_seq_len) ) 
    {
      print "@"$3" length="length(rna_seq), rna_seq, "+", rna_qual >> (prefix"_1.fastq");
      print "@"$3" length="length(dna_seq), dna_seq, "+", dna_qual >> (prefix"_2.fastq");
      passed_reads_counter++;
      reads_for++;
    }
  }

  if  ( ($1==not_found) && ($6==r_br)  ) 
  {
    dna_seq=substr($9, $7+bridge_len, 999);
    dna_qual=substr($10, $7+bridge_len, 999);
    rna_seq=substr($9, 1, $7-1);
    rna_qual=substr($10, 1, $7-1);

    print "@"$8" length="length(rna_seq), revcomp(rna_seq), "+", revcomp(rna_qual) >> (prefix"_1.fastq");
    print "@"$8" length="length(dna_seq), revcomp(dna_seq), "+", revcomp(dna_qual) >> (prefix"_2.fastq");
    passed_reads_counter++;
    reads_rev++;
  }
}'



echo "Passed Reads Counter: $passed_reads_counter" > $prefix.bridge_PE.stats.txt
echo "Reads Forward: $reads_for" >> $prefix.bridge_PE.stats.txt
echo "Reads Reverse: $reads_rev" >> $prefix.bridge_PE.stats.txt