#!/bin/bash -l

conda activate faster_bardic
module load python/python-3.8.2

outdir="./test_files/results"

#разметка для человека из базы данных RNAChrom DB, которая прошла через скрипт валидации
annot="test_annot/All_RNAs_hS_DB_pipe_new.0-corrected_annot.tab"
rnaparts_no_srr="rna_parts/GSM5315226.chr22.no_srr.tab"
dnaparts="dna_parts/GSM5315226.chr22.dna_parts_no_srr.tab"


blk="./test_files/blk/hg38.blacklist.bed"
annot_GTF="./test_files/annot/gencode.v43.annotation.gene.gtf"

Smoother="/home/suvorova/tools/src/Smoother"
cfg="/home/anikolskaya/RDC/rnachrom_pipeline_faster/configs/stereogene/cfg.cfg"
chrsizes="/home/anikolskaya/general/chrsizes/hg38_canonical_chromsizes.tsv"

source ~/.bashrc


# РНК-ЧАСТИ
#по умолчанию файл с header (имена столбцов могут быть любыми), можно указать опцию --no_header, как ниже
#в сплайсинге добавить опцию --rna_parts
echo 'splicing'
conda run -n rnachrom_pipe python3 splicing.py test_files/$rnaparts_no_srr --no_header --rna_parts --outdir $outdir

echo 'voting by rna parts only'

rnaparts_srr="results/GSM5315226.chr22.1-CIGAR_srr.tab"
conda run -n rnachrom_pipe python3 annotation_voting.py test_files/$rnaparts_srr test_files/$annot --rna_parts --outdir $outdir

voted="results/GSM5315226.chr22.4-voted.tab"
#для мержда применить команду 
join -t $'\t' -j4 -o2.1,2.2,2.3,2.4,2.5,1.1,1.2,1.3,1.5,2.6,2.7,2.8,2.9,2.10 <(sort -k4 test_files/$dnaparts) <(sort -k4 test_files/$voted) > $outdir/GSM5315226.4-voted_rdc.tab

#ДНК-ЧАСТИ
echo 'blacklist'
#python3 blacklisting.py $outdir/GSM4041593_no_srr.1-CIGAR.tab $blk --outdir $outdir

echo 'voting'
#аннотировать удобнее похромосомно
#for chr in chrY chr22; do #chr{1.22} chrX chrY
#python3 annotation_voting.py $indir/RedC_K562_${chr}.tab $annot --annot_format BED --outdir $outdir
#не выводить статистику (считается гораздо быстрее)
#python3 annotation_voting.py $indir/RedC_K562_${chr}.tab $annot --annot_format BED --outdir $outdir --no_stat  # --cpus 4
#python3 annotation_voting.py $indir/RedC_K562_${chr}.tab $annot_GTF --annot_format GTF --outdir $outdir
#done



#infile_name="путь/до/файла/после/merge/RedC_K562.tab"

#здесь нужен полный файл
#python3 background.py $full_name $Smoother $cfg $chrsizes --filter_top_n 50 --filter_tail_n 1000 --outdir $outdir

echo 'raw_N2'
#for chr in chr22 chrY; do
#python3 normalize_raw.py $outdir/RedC_K562_${chr}.4-voted.tab $outdir/K562_pipe.5-background_sm.bgr --outdir $outdir #--cpus 4
#done 

#надо объединить скрипты со статистикой
#awk 'FNR>1 || NR==1' $outdir/*N2_raw.stat.tab > $outdir/$(basename "${infile_name%.*}").5-N2_raw.stat.tab

echo 'final N2'
#прохромосомно пересчитать, подать на вход файл со статистикой (параметр --by_chr), полученный выше
#for chr in chr22 chrY; do
#python3 normalize_N2.py $outdir/RedC_K562_${chr}.5-N2_raw.tab $chrsizes --by_chr $outdir/$(basename "${infile_name%.*}").5-N2_raw.stat.tab --outdir $outdir
#done

echo 'scaling'
#for chr in chr22 chrY; do
#python3 scaling.py $outdir/RedC_K562_${chr}.5-N2.tab $chrsizes --outdir $outdir
#done
