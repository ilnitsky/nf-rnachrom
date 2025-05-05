#!/bin/bash

export LC_NUMERIC="C"
set -e

if [ $# -lt 2 ]; then
    echo "No options found!"
    exit 1
fi


while getopts "p:g:d:o:" opt; do
    case $opt in
    p)
        GENESPATH=$OPTARG
        # echo "GENES PATH: $GENESPATH"
        ;;
    g)
        GENES=$OPTARG
        # echo "GENES FILE: $GENES"
        ;;
    d)
        DISTANCE=$OPTARG
        # echo "DISTANCE: $DISTANCE, it should be greater, than read size"
        ;;
    o)
        OUTPUT_DIR=$OPTARG
        # echo "OUTPUT DIR: $OUTPUT_DIR"
    esac
done

if [[ ! -f $GENESPATH/${GENES} ]] ; then
    echo 'File $GENESPATH/$GENES is not there, aborting.'
    exit
fi

mkdir -p $OUTPUT_DIR

awk '($5=="+")' $GENESPATH/$GENES | awk 'BEGIN { OFS = "\t"}; {print $1, $2, $3, $4, 0, $5, $6, $7}' >$OUTPUT_DIR/$GENES.genes_pos_strand.bed
bedtools merge -i $OUTPUT_DIR/$GENES.genes_pos_strand.bed -d $DISTANCE | awk -F '\t' 'BEGIN { OFS = "\t"}; {print $0, NR}' > $OUTPUT_DIR/$GENES.clusters_dist_${DISTANCE}_pos_strand.bed
bedmap --echo --echo-map-id --delim '\t' $OUTPUT_DIR/$GENES.genes_pos_strand.bed $OUTPUT_DIR/${GENES}.clusters_dist_${DISTANCE}_pos_strand.bed >$OUTPUT_DIR/${GENES}.genes_pos_strand.clusters_dist_${DISTANCE}.bed

awk '($5=="-")' $GENESPATH/$GENES | awk 'BEGIN { OFS = "\t"}; {print $1, $2, $3, $4, 0, $5, $6, $7}' >$OUTPUT_DIR/$GENES.genes_neg_strand.bed
bedtools merge -i $OUTPUT_DIR/$GENES.genes_neg_strand.bed -d $DISTANCE | awk -F '\t' 'BEGIN { OFS = "\t"}; {print $0, NR}' >$OUTPUT_DIR/$GENES.clusters_dist_${DISTANCE}_neg_strand.bed
bedmap --echo --echo-map-id --delim '\t' $OUTPUT_DIR/$GENES.genes_neg_strand.bed $OUTPUT_DIR/${GENES}.clusters_dist_${DISTANCE}_neg_strand.bed >$OUTPUT_DIR/${GENES}.genes_neg_strand.clusters_dist_${DISTANCE}.bed