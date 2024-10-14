#!/bin/bash
base=$(basename -s .fastq pear.assembled.fastq)
bridge_len=37

#get max length of merged reads and add Ns to reads shorter 
#then max length for rklib to work

max_len=$(awk '{if(NR%4==2) {print length($0)}}' $base.fastq | sort -nr | head -1)
awk -v max_len=$max_len '
BEGIN {OFS="\n"}
{
    if(NR%4==1) {header=$0}    # header line
    if(NR%4==2) {              # sequence line
        seq=$0
        len=length(seq)
        for(i=len+1;i<=max_len;i++) {
            seq=seq"N"         # append Ns to the sequence
        }
    }
    if(NR%4==3) {plus=$0}     # plus line
    if(NR%4==0) {              # quality line
        qual=$0
        len=length(qual)
        for(i=len+1;i<=max_len;i++) {
            qual=qual"I"       # append Is to the quality
        }
        print header, seq, plus, qual
    }
}' $base.fastq > $base.N.fastq
echo $max_len


# fasta2hash redc_bridge_for.fasta redc_bridge_for.bin
fastq2hash $base.N.fastq $base.bin


# Extract the first three fields from the fastq file
bioawk -c fastx '{print}' $base.N.fastq | cut -f1-3 | \

# Append the output of two rk_querysearch commands
paste - <(rk_querysearch redc_bridge_for.bin $base.bin 37 $max_len 2 0 200 1 | tail -n+2 | cut -f4) | \

# Process the combined output with bioawk
# Bridge design:
# (5') AANNN-AAACCGGCGTCCAAG (3')
# (3') TTNNN-TTTGGCCGCAGGTTC (5')
# Bridge sequence used for search is AAACCGGCGTCCAAG. Thus, we need to substract 5 nts from found bridge start coordinate.

bioawk -v bridge_len=37 'BEGIN{OFS="\n";} 
{   
    rna=substr($2, 1, $4-1);
    dna=substr($2, $4+bridge_len+1, 999);

    if  ( !($4==99999) )  
    {
        print $1, substr($2, 1, $4), "+", substr($3, 1, $4) > "dna.fastq"; 
        print $1, substr($2, $4+bridge_len+1, 999), "+", substr($3, $4+bridge_len+1, 999) > "rna.fastq"
    }  
}'
