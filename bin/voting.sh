#!/bin/bash

export LC_NUMERIC="C"
set -e

if [ $# -lt 2 ]; then
    echo "No options found!"
    exit 1
fi

while getopts "s:q:g:d:e:c:o:" opt; do
    case $opt in
    s)
        SCRIPTSPATH=$OPTARG
        echo "SCRIPTS PATH: $SCRIPTSPATH"
        ;;
    q)
        GENESPATH=$OPTARG
        echo "GENES PATH: $GENESPATH"
        ;;
    g)
        GENES=$OPTARG
        echo "GENES FILE: $GENES"
        ;;
    d)
        DISTANCE=$OPTARG
        echo "DISTANCE: $DISTANCE, clustering distance"
        ;;
    e)
        CONTACTSPATH=$OPTARG
        echo "CONTACTS PATH: $CONTACTSPATH"
        ;;
    c)
        CONTACTS_FILE=$OPTARG
        echo "CONTACTS FILE $CONTACTS_FILE"
        ;;
    o)
        OUTPUT_DIR=$OPTARG
        echo "OUTPUT DIR: $OUTPUT_DIR"
        ;;
    esac
done

if [[ ! -f $CONTACTSPATH/$CONTACTS_FILE ]] ; then
    echo 'File $CONTACTSPATH/$CONTACTS_FILE is not there, aborting.'
    exit
fi

echo Creating clusters..
bash $SCRIPTSPATH/create_clusters.sh -p $GENESPATH -g $GENES -d $DISTANCE -o $OUTPUT_DIR/genes

CLUSTERS_POS_STRAND=$OUTPUT_DIR/genes/$GENES.clusters_dist_${DISTANCE}_pos_strand.bed
CLUSTERS_NEG_STRAND=$OUTPUT_DIR/genes/$GENES.clusters_dist_${DISTANCE}_neg_strand.bed

GENES_CLUSTERS_POS_STRAND=$OUTPUT_DIR/genes/${GENES}.genes_pos_strand.clusters_dist_${DISTANCE}.bed
GENES_CLUSTERS_NEG_STRAND=$OUTPUT_DIR/genes/${GENES}.genes_neg_strand.clusters_dist_${DISTANCE}.bed

echo Sorting contacts and mapping to clusters..

#add new columns to contacts file
awk 'BEGIN{OFS="\t"} NR==1 {
    # Первые 6 колонок
    for(i=1;i<=6;i++) printf "%s\t", $i;
    # Добавляем новые колонки
    printf "gene_name\tgene_type\tfrom_source\t";
    # Оставшиеся колонки до 20-й
    for(i=7;i<=20;i++) {
        if(i==20) print $i;
        else printf "%s\t", $i;
    }
}' $CONTACTSPATH/$CONTACTS_FILE > $OUTPUT_DIR/contacts.voting.UU.bed

awk 'BEGIN{OFS="\t"} NR==1 {
    for(i=1;i<=20;i++) {
        if(i==20) print $i;
        else printf "%s\t", $i;
    }
}' $CONTACTSPATH/$CONTACTS_FILE > $OUTPUT_DIR/singletons.UU.bed

cp $OUTPUT_DIR/contacts.voting.UU.bed $OUTPUT_DIR/contacts.voting.UM.bed
cp $OUTPUT_DIR/singletons.UU.bed $OUTPUT_DIR/singletons.UM.bed

mkdir -p $OUTPUT_DIR/tmp

awk 'BEGIN { OFS = "\t"}; NR == 1 { next } {print $3, $4+int(($5-$4)/2), $4+int(($5-$4)/2)+1, $1, 0, $6, $4 ,$5, $2, $7, $8, $9 , $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' $CONTACTSPATH/$CONTACTS_FILE | sort -k1,1 -k2,2n --temporary-directory=$OUTPUT_DIR/tmp > $OUTPUT_DIR/contacts.sorted.bed

awk '($6=="+")' $OUTPUT_DIR/contacts.sorted.bed | bedmap --echo --echo-map-id --delim '\t' --unmapped-val -1 - $CLUSTERS_POS_STRAND > $OUTPUT_DIR/contacts.pos_strand.clusters.bed
awk '($6=="-")' $OUTPUT_DIR/contacts.sorted.bed | bedmap --echo --echo-map-id --delim '\t' --unmapped-val -1 - $CLUSTERS_NEG_STRAND > $OUTPUT_DIR/contacts.neg_strand.clusters.bed

echo Voting...
python $SCRIPTSPATH/voting.py -o $OUTPUT_DIR -p $GENES_CLUSTERS_POS_STRAND -n $GENES_CLUSTERS_NEG_STRAND

# Split voting results into UU and UM files
echo "Splitting voting results by UU/UM..."

# Combine positive and negative strands into single UU and UM files
awk -F'\t' '($2=="UU" || $2=="U")' $OUTPUT_DIR/contacts.pos_strand.voting.bed $OUTPUT_DIR/contacts.neg_strand.voting.bed >> $OUTPUT_DIR/contacts.voting.UU.bed
awk -F'\t' '($2=="UM")' $OUTPUT_DIR/contacts.pos_strand.voting.bed $OUTPUT_DIR/contacts.neg_strand.voting.bed >> $OUTPUT_DIR/contacts.voting.UM.bed

# Split singletons into UU and UM files
awk -F'\t' '($2=="UU" || $2=="U")' $OUTPUT_DIR/singletons.pos_strand.bed $OUTPUT_DIR/singletons.neg_strand.bed >> $OUTPUT_DIR/singletons.UU.bed
awk -F'\t' '($2=="UM")' $OUTPUT_DIR/singletons.pos_strand.bed $OUTPUT_DIR/singletons.neg_strand.bed >> $OUTPUT_DIR/singletons.UM.bed


# Combine batki files and calculate new value
python $SCRIPTSPATH/combine_batki.py --output-dir $OUTPUT_DIR

# Add missing genes
echo "Adding missing genes..."
python $SCRIPTSPATH/add_missing_genes.py --genes-file $GENESPATH/$GENES --output-dir $OUTPUT_DIR

# Plot gene type distribution
echo "Plotting gene type distribution..."
python $SCRIPTSPATH/make_plots.py -i $OUTPUT_DIR/counts.tsv -o $OUTPUT_DIR

# Clean up intermediate files
rm -f $OUTPUT_DIR/contacts.sorted.bed
rm -f $OUTPUT_DIR/contacts.pos_strand.clusters.bed
rm -f $OUTPUT_DIR/contacts.neg_strand.clusters.bed
rm -f $OUTPUT_DIR/contacts.pos_strand.voting.bed
rm -f $OUTPUT_DIR/contacts.neg_strand.voting.bed
rm -f $OUTPUT_DIR/singletons.pos_strand.bed
rm -f $OUTPUT_DIR/singletons.neg_strand.bed
rm -rf $OUTPUT_DIR/genes
rm -rf $OUTPUT_DIR/tmp

