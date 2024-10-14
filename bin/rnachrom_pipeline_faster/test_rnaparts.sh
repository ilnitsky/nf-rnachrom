#!/bin/bash -l

outdir="./test_files/results"

#разметка для человека из базы данных RNAChrom DB, которая прошла через скрипт валидации
annot="test_annot/All_RNAs_hS_DB_pipe_new.0-corrected_annot.tab"
rnaparts_no_srr="rna_parts/GSM5315226.chr22.no_srr.tab"
dnaparts="dna_parts/GSM5315226.chr22.dna_parts_no_srr.tab"


source ~/.bashrc


# РНК-ЧАСТИ
#по умолчанию файл с header (имена столбцов могут быть любыми), можно указать опцию --no_header, как ниже
#в сплайсинге добавить опцию --rna_parts, тогда от файла ожидаются столбцы ['rna_chr', 'rna_bgn', 'rna_end', 'id', 'rna_strnd', 'rna_cigar']
echo 'splicing'
conda run -n rnachrom_pipe python3 splicing.py test_files/$rnaparts_no_srr --no_header --rna_parts --outdir $outdir

#также указание опции --rna_parts, ожидаются столбцы ['rna_chr', 'rna_bgn', 'rna_end', 'id', 'rna_strnd', 'rna_cigar', 'srr_id']
echo 'voting by rna parts only'
rnaparts_srr="results/GSM5315226.chr22.1-CIGAR_srr.tab"
conda run -n rnachrom_pipe python3 annotation_voting.py test_files/$rnaparts_srr test_files/$annot --rna_parts --outdir $outdir

voted="results/GSM5315226.chr22.4-voted.tab"
#для мержда применить команду: мердж по столбцу 4 ('id'). К ДНК-частям присвоятся результаты голосования и РНК-координаты после поправки на сплайсинг.
join -t $'\t' -j4 -o2.1,2.2,2.3,2.4,2.5,1.1,1.2,1.3,1.5,2.6,2.7,2.8,2.9,2.10 <(sort -k4 test_files/$dnaparts) <(sort -k4 test_files/$voted) > $outdir/GSM5315226.4-voted_rdc.tab