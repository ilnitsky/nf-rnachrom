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
min_seq_len=$6                # RNA/DNA parts length filter
prefix=$7                     # SRR prefix
not_found=$8                  # bridge not found identifier: chartools = 1, bitap = 0
f_br=$9                       # forward bridge identifier: chartools = 2, bitap = 1
r_br=${10}                    # reverse bridge identifier: chartools = 3, bitap = 2

passed_reads_counter=0


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
  <(cat ${unassembled_positions_F}) \
  <(bioawk -c fastx '{print}' ${unassembled_F} | cut -f1-3) \
  <(cat ${unassembled_positions_R}) \
  <(bioawk -c fastx '{print}' ${unassembled_R} | cut -f1-3) | \
  
  bioawk -v bridge_len=${bridge_len} -v min_seq_len=${min_seq_len} -v prefix=${prefix} -v OFS="\n" '
  function process_seq(seq, qual, pos, is_rev, file_suffix) {
    rna_seq  = substr(seq, 1, pos-1);
    rna_qual = substr(qual, 1, pos-1);
    dna_seq  = substr(seq, pos+bridge_len);
    dna_qual = substr(qual, pos+bridge_len);
    if (is_rev) {
      rna_seq  = revcomp(rna_seq);
      rna_qual = revcomp(rna_qual);
      dna_seq  = revcomp(dna_seq);
      dna_qual = revcomp(dna_qual);
    }
    if (length(rna_seq) >= min_seq_len && length(dna_seq) >= min_seq_len) {
      print "@"$3" length="length(rna_seq), rna_seq, "+", rna_qual >> (prefix"_1"file_suffix".fastq");
      print "@"$3" length="length(dna_seq), dna_seq, "+", dna_qual >> (prefix"_2"file_suffix".fastq");
    }
  }
  {
    if (($1 == f_br && $6 == not_found) || ($1 == not_found && $6 == r_br)) {
      process_seq($4, $5, $2, 0, "");
      passed_reads_counter++;
    } else if (($1 == r_br && $6 == not_found) || ($1 == not_found && $6 == r_br)) {
      process_seq($9, $10, $7, 1, "_revcomp");
      passed_reads_counter++;
    }
  }'

echo "Passed Reads Counter: $passed_reads_counter" > $prefix.stats.txt
