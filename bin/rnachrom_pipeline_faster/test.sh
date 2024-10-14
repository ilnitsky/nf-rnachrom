
basename="SRR9201799_srr"
#basename="SRR9201803"
#basename="SRR9201801"

outdir="/home/anikolskaya/rnachrom_pipeline/final_check"
blcklst="/home/anikolskaya/general/blacklist/mm10.blacklist.bed"
#annot_BED="/home/anikolskaya/rnachrom_pipeline/test_files/genes_set5.pipe.bed"
#annot_BED="/home/anikolskaya/general/annot_DB/All_RNAs_hS_DB_pipe_new.bed"
annot_BED="/home/anikolskaya/rnachrom_pipeline/vanya_code/mus.bed"
annot_GTF="/home/anikolskaya/rnachrom_pipeline/test_files/gencode.v43.annotation.gene.gtf"

Smoother="/home/anikolskaya/tools/stereogene/src/Smoother"
cfg="/home/anikolskaya/rnachrom_pipeline/test_files/stereogene/cfg.cfg"
chrsizes="/home/anikolskaya/general/chrsizes/mm10_canonical_chromsizes.tsv"

#for file in /home/anikolskaya/RedC/check_N2; do
#python3 normalize_N2.py $file /home/anikolskaya/tools/stereogene/src/Smoother /home/anikolskaya/rnachrom_pipeline/test_files/stereogene/cfg.cfg /home/anikolskaya/general/chrsizes/hg38_canonical_chromsizes.tsv --filter_top_n 50 --filter_tail_n 1000 --outdir results_RedC

#python3 annotation_voting.py /home/anikolskaya/rnachrom_pipeline/extra/merged_0_60_results/merged_0h1.2-blacklist_sample.tab  $annot_BED --annot_format BED --cpus 6 --outdir results_annot
#python3 annotation_voting.py /home/anikolskaya/rnachrom_pipeline/extra/merged_0_60_results/merged_0h1.2-blacklist_sample.tab  $annot_GTF --annot_format GTF --cpus 6 --outdir results_ann

#python3 splicing.py -h
#python3 blacklisting.py -h
#python3 annotation_voting.py -h
#python3 normalize_N2.py -h
#python3 scaling.py -h 
#cd /home/anikolskaya/rnachrom_pipeline-new/rnachrom_pipeline
#lida_file="/home/garkul/RNA_parts_project/RNA-seq/data/mapped_SRR4422606_filtered_with_chr_header_andSRR.1-CIGAR.tab"
#vanya_file="/home/ilnitsky/nf-rnachrom/work/c1/43ce1bf596b3de292c0cb84f9b6da9/60h_merged.tab"
vanya_file="/home/anikolskaya/rnachrom_pipeline-new/rnachrom_pipeline_faster/test_files/60h_voted_chr1.4-voted.tab"
#python3 splicing.py /home/anikolskaya/rnachrom_pipeline/test_files/$basename.tab --outdir $outdir
#python3 splicing.py /home/ilnitsky/nf-rnachrom/data/contacts/$basename.tab --outdir $outdir
#python3 blacklisting.py $outdir/$basename.1-CIGAR.tab $blcklst --outdir $outdir
#python3 annotation_voting.py $outdir/$basename.2-blacklist.tab $annot_GTF --annot_format GTF --outdir $outdir
#python3 annotation_voting.py $outdir/$basename.2-blacklist.tab $annot_BED --annot_format BED --outdir $outdir
#python3 annotation_voting_reduced.py $vanya_file $annot_BED --annot_format BED --outdir results
python3 background.py $vanya_file $Smoother $cfg $chrsizes --filter_top_n 50 --filter_tail_n 1000 --outdir results
python3 normalize_raw.py $vanya_file results/60h_voted_chr1.5-background_sm.bgr --outdir results
python3 normalize_N2_reduced.py results/60h_voted_chr1.5-N2.tab results/60h_voted_chr1.5-N2_raw.stat.tab --outdir results
python3 normalize_N2.py $outdir/$basename.4-voted_noribo.tab $Smoother $cfg $chrsizes --filter_top_n 50 --filter_tail_n 1000 --outdir $outdir
python3 scaling_reduced.py results/60h_voted_chr1.5-N2.tab $chrsizes --outdir results

#python3 scaling.py $outdir/$basename.5-N2.tab $chrsizes --threads 6 --outdir $outdir 
#python3 annotation_voting.py radicl_1fa_merged_corr.tab $annot_BED --annot_format BED --cpus 6 --outdir $outdir
#python3 blacklisting.py results_new/SRR9201799.1-CIGAR.tab /home/anikolskaya/rnachrom_pipeline/test_files/blacklist/mm10.blacklist.bed --outdir results_new
#python3 annotation_voting.py results_new/SRR9201799.2-blacklist.tab test_files/genes_set5.pipe.bed --annot_format BED --outdir results_new
#python3 annotation_voting.py results_new/SRR9201799.2-blacklist.tab test_files/gencode.v43.annotation.gene.gtf --annot_format GTF --outdir results_gtf
#python3 normalize_N2.py /home/anikolskaya/rnachrom_pipeline/extra/results/SRR9201799.4-voted_noribo.tab /home/anikolskaya/tools/stereogene/src/Smoother /home/anikolskaya/rnachrom_pipeline/test_files/stereogene/cfg.cfg /home/anikolskaya/general/chrsizes/mm10_canonical_chromsizes.tsv --filter_top_n 50 --filter_tail_n 1000 --outdir results_new
#python3 normalize_N2.py /home/rnachrom2/data/grid_2023_processing/60h_merged.4-voted.tab /home/anikolskaya/tools/stereogene/src/Smoother /home/anikolskaya/rnachrom_pipeline/test_files/stereogene/cfg.cfg /home/anikolskaya/general/chrsizes/mm10_canonical_chromsizes.tsv --filter_top_n 50 --filter_tail_n 1000 --outdir results_new
#python3 scaling.py /home/anikolskaya/rnachrom_pipeline/results_new/SRR9201799.5-N2.tab /home/anikolskaya/chrsizes/mm10_canonical_chromsizes.tsv --threads 6 --outdir results_new

#python3 annotation_voting.py /home/anikolskaya/rnachrom_pipeline/radicl_1fa_results/SRR9201799.2-blacklist.tab $annot_BED --annot_format BED --cpus 6 --outdir /home/anikolskaya/rnachrom_pipeline/radicl_1fa_results 
#outdir="results_cnts"
#basename="SRR9201799_sample"
#endogeneous background track
#python3 peak-calling.py input /home/anikolskaya/rnachrom_pipeline/extra/radicl_1fa_results/${basename}.4-voted_noribo.tab test_files/genes_set5.pipe.bed $chrsizes --annot_file_format BED  --outdir $outdir 
#python3 peak-calling.py input /home/anikolskaya/rnachrom_pipeline/extra/radicl_1fa_results/${basename}.4-voted_noribo.tab test_files/All_RNAs_mM_DB_pipe.bed $chrsizes --annot_file_format BED --outdir $outdir
#python3 peak-calling.py input /home/anikolskaya/rnachrom_pipeline/extra/radicl_1fa_results/${basename}.4-voted_noribo.tab test_files/genes_set5.pipe.bed $chrsizes --annot_file_format BED --outdir $outdir --bgtype rnas_list --bgdata ./test_files/background/rna_list.txt
#python3 peak-calling#.py input /home/anikolskaya/rnachrom_pipeline/extra/radicl_1fa_results/${basename}.4-voted_noribo.tab test_files/genes_set5.pipe.bed $chrsizes --annot_file_format BED --outdir $outdir --bgtype input_contacts --bgdata ./test_files/background/pc_bed3.bed
# python3 peak-calling.py input /home/anikolskaya/rnachrom_pipeline/extra/radicl_1fa_results/${basename}.4-voted_noribo.tab $annot_GTF $chrsizes --annot_file_format GTF --outdir $outdir --bgtype input_contacts --bgdata ./test_files/background/pc_bed3.bed


#python3 peak-calling.py bed2h5 $outdir/${basename}.7-contacts_for_peaks.bed ./results_cnts/All_RNAs_mM_DB_pipe.7-annot_for_peaks.bed $chrsizes ${outdir}/${basename}.8-DnaDataset.dnah5
#python3 peak-calling.py binsizes $outdir/${basename}.8-DnaDataset.dnah5 $outdir/${basename}_sample.8-selection.tsv @configs/binsizes_ATA.config --min_contacts 1500 --trans_min 5000
#python3 peak-calling.py binsizes @configs/binsizes_OTA.config
#python3 peak-calling.py makerdc $outdir/${basename}.8-DnaDataset.dnah5 $outdir/${basename}.7-background.BedGraph $outdir/${basename}.8-contacts.rdc @configs/makerdc.config 
#python3 peak-calling.py scaling $outdir/${basename}.8-contacts.rdc @configs/scaling.config
#python3 peak-calling.py peaks $outdir/${basename}.8-contacts.rdc $outdir/${basename}.8-peaks.narrowPeak @configs/peaks.config
