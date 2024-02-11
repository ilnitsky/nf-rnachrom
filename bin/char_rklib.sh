bridge_len=19

base=$(basename -s .fastq SRR25582110_1.400.fastq)

fasta2hash char_bridge_for.fasta char_bridge_for.bin
fasta2hash char_bridge_rev.fasta char_bridge_rev.bin
fastq2hash $base.fastq $base.bin


# Extract the first three fields from the fastq file
bioawk -c fastx '{print}' $base.fastq | cut -f1-3 | \

# Append the output of two rk_querysearch commands
paste - <(rk_querysearch char_bridge_for.bin $base.bin $bridge_len 150 2 0 140 1 | tail -n+2 | cut -f4) \
<(rk_querysearch char_bridge_rev.bin $base.bin $bridge_len 150 2 0 140 1 | tail -n+2 | cut -f5) | \

# Process the combined output with bioawk
# Bridge design:
# (5') AANNN-AAACCGGCGTCCAAG (3')
# (3') TTNNN-TTTGGCCGCAGGTTC (5')
# Bridge sequence used for search is AAACCGGCGTCCAAG. Thus, we need to substract 5 nts from found bridge start coordinate.

bioawk -v bridge_len="$bridge_len" 'BEGIN{OFS="\n";} 
{   
    rna=substr($2, 1, $4-5);
    dna=substr($2, $4+bridge_len+1, 999);
    rna_rev=substr($2, $5+6, 999);
    dna_rev=substr($2, 1, $5-bridge_len);

    if ( (!($5<99999 && $4<99999)) || (!($4==99999 && $5==99999)) )
    {
        if  ( ($5==99999) && (length(rna)>=21) &&  (length(dna)>=21) )  
        {
            print $1, substr($2, 1, $4-5), "+", substr($3, 1, $4-5) > "rna.fastq"; 
            print $1, "GATC"substr($2, $4+bridge_len+1, 999), "+", "FFFF"substr($3, $4+bridge_len+1, 999) > "dna.fastq"
        }  
        if ( ($4==99999) && (length(dna_rev)>=21) &&  (length(rna_rev)>=21) )  
        {
            print $1, "GATC"revcomp(substr($2, 1, $5-bridge_len)), "+", "FFFF"revcomp(substr($3, 1, $5-bridge_len)) > "dna.fastq"; 
            print $1, revcomp(substr($2, $5+6, 999)), "+", revcomp(substr($3, $5+6, 999)) > "rna.fastq"
        } 
    }
}'
