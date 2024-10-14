#!/bin/bash

bwa_index=/home/ryabykhgrigory/nf-rnachrom/data/genome_index/bwa/GRCh38.p13
hisat_index=/home/ryabykhgrigory/nf-rnachrom/clean_genomes/hg38/GRCh38.p13
hisat_SS=/home/ryabykhgrigory/nf-rnachrom/data/genes/gencode.v43.ss

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 input_dna_file input_rna_file output_dna_file output_rna_file"
    exit 1
fi

input_dna_file=$1
input_rna_file=$2
output_dna_file=$3
output_rna_file=$4


if [ ! -f "$input_dna_file" ]; then
    echo "Error: $input_dna_file does not exist."
    exit 1
fi

if [ ! -f "$input_rna_file" ]; then
    echo "Error: $input_rna_file does not exist."
    exit 1
fi

input_dna_prefix="${input_dna_file%%.*}"
input_rna_prefix="${input_rna_file%%.*}"
output_dna_prefix="${output_dna_file%%.*}"
output_rna_prefix="${output_rna_file%%.*}"


bwa mem -t 12 -5M $bwa_index  \
  $input_rna_file > \
  $input_rna_prefix.sam 2> $input_rna_prefix.bwa.log

bwa mem -t 12 -5M $bwa_index \
  $input_dna_file > \
  $input_dna_prefix.sam 2> $input_dna_prefix.bwa.log

samtools view -h -F 256 $input_rna_prefix.sam | python ~/gpfs/nf-core-rnachrom/bin/extract_sam_file_matched_seq_to_fastq2.py > $input_rna_prefix.removed_softclip.fastq
samtools view -h -F 256 $input_dna_prefix.sam | python ~/gpfs/nf-core-rnachrom/bin/extract_sam_file_matched_seq_to_fastq2.py > $input_dna_prefix.removed_softclip.fastq

hisat2 -p 12 -x $hisat_index \
  -U $input_rna_prefix.removed_softclip.fastq -k 100 --no-softclip --known-splicesite-infile $hisat_SS --dta-cufflinks \
  -S $output_rna_file 2> $input_rna_prefix.hisat.log

hisat2 -p 12 -x $hisat_index \
  -U $input_dna_prefix.removed_softclip.fastq --no-spliced-alignment -k 100 --no-softclip -S $output_dna_file 2> $input_dna_prefix.hisat.log

