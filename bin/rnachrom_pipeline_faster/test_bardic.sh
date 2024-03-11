#!/bin/bash -l
#SBATCH --job-name=bardic_test 
#SBATCH --error=bardic_test.err  
#SBATCH --output=bardic_test.log 
#SBATCH --time=12:00:00  
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=2 
#SBATCH --mem=20G 

scripts="/home/anikolskaya/gpfs/rnachrom_pipeline_faster"

# полный файл с контактами эксперимента RedC K562
rdc="/home/anikolskaya/gpfs/RDC/RedC/3_annot/set8/hg38/K562/voted_merged/K562.4-voted_merged.tab"
# файл с некоторыми РНК для теста
#rdc_selected="/home/anikolskaya/gpfs/rnachrom_pipeline_faster/test_files/for_peaks/K562.selected-for_peaks.bed"
# для теста возьмем аннотацию set8 с энхансерными РНК. Это bed6, который я создаю в скрипте, который валидирует файл с разметкой РНК.
annotation="/home/anikolskaya/gpfs/general/annot_hs/DB_annot/corrected_annot/overlap_0.95_gencode/All_RNAs_hS_DB_pipe_new.0-corrected_annot_no_header.bed"
chromsizes="/home/anikolskaya/gpfs/general/chrsizes/hg38_canonical_chromsizes.tsv"
outdir="/home/anikolskaya/gpfs/rnachrom_pipeline_faster/test_files/for_peaks"

mkdir -p $outdir

source ~/.bashrc

echo 'peak-calling'

# готовим .txt файл с именами белок-кодирующих РНК. На его основе Бардик потом построит фон.
grep -w 'protein_coding' $rdc | awk '{print $11}' | sort | uniq > $outdir/$(basename $rdc .4-voted_merged.tab).4-pc.txt

# из файла с контактами готовим bed6 headerless файл для Бардика
sed 1d $rdc | awk -F"\t" '{OFS=FS} {print $6,$7,$8,$11,".",$9};' > $outdir/$(basename $rdc .4-voted_merged.tab).4-for_peaks.bed

# синтаксис: bardic run $dnaparts $annotation $chromsizes $bgdata $outdir  --qval_threshold 1 --cores 10
# запуск только для ATA!

conda run -n bardic bardic run $outdir/K562.4-for_peaks.bed \
                               $annotation \
                               $chromsizes \
                               $outdir/K562.4-pc.txt \
                               $outdir/peaks  --cores 2 @configs/run.config


# conda run -n bardic bardic run $outdir/$(basename $rdc_selected .tab).4-for_peaks.bed \
#                                $annotation \
#                                $chromsizes \
#                                $outdir/$(basename $rdc .4-voted_merged.tab).4-pc.txt \
#                                $outdir/peaks  --qval_threshold 1 --cores 1

# с конфигом
# conda run -n bardic bardic run $outdir/$(basename $rdc_selected .tab).4-for_peaks.bed \
#                                $annotation \
#                                $chromsizes \
#                                $outdir/$(basename $rdc .4-voted_merged.tab).4-pc.txt \
#                                $outdir/peaks  @configs/run.config --cores 1

# для 