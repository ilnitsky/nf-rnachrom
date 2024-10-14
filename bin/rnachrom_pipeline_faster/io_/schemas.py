"""
Field names for various genomic tabular files

"""


rdc_BED = [
    'rna_chr',
    'rna_bgn',
    'rna_end',
    'id',
    'rna_strand',
    'rna_cigar',
    'dna_chr',
    'dna_bgn',
    'dna_end',
    'dna_strand',
    'dna_cigar'
]

rdc_BED_sample = [
    'rna_chr',
    'rna_bgn',
    'rna_end',
    'id',
    'rna_strand',
    'rna_cigar',
    'dna_chr',
    'dna_bgn',
    'dna_end',
    'dna_strand',
    'dna_cigar',
    'srr_id'
]

rdc_dtypes = {
    'rna_chr': 'category',
    'rna_bgn': 'int',
    'rna_end': 'int',
    'id': 'str',
    'rna_strand': 'category',
    'dna_chr': 'category',
    'dna_bgn': 'int',
    'dna_end': 'int',
    'dna_strand': 'category',
    'rna_cigar': 'str',
    'dna_cigar': 'str',
    'srr_id': 'category',
    'gene_name': 'str',
    'gene_type': 'category',
    'source': 'category',
    'bg_sm': 'float',
    'N2_raw': 'float',
    'N2': 'float',
    'SCnorm': 'float'
    
}

synonyms_annot = [
    'name'
    'source', 
    'synonym_name',
    'synonym_source',
    'overalap_coef'
]

annotation =[
    'chrom',
    'start',
    'end',
    'name',
    'strand',
    'gene_type',
    'source',
]

annot_BED =[
    'chrom',
    'start',
    'end',
    'name',
    'strand',
    'gene_type',
    'source',
]

anotation_dtypes = {
    'chrom': 'category',
    'start': 'int',
    'end': 'int',
    'name': 'str',
    'strand': 'category',
    'gene_type': 'str',
    'source': 'category'    
}

annot_BED_dtypes = {
    'chrom': 'category',
    'start': 'int',
    'end': 'int',
    'name': 'str',
    'strand': 'category',
    'gene_type': 'str',
    'source': 'category'    
}


voted_BED = [
    'rna_chr',
    'rna_bgn',
    'rna_end',
    'id',
    'rna_strand',
    'dna_chr',
    'dna_bgn',
    'dna_end',
    'dna_strand',
    'srr_id',
    'gene_name',
    'gene_type',
    'source' 
]

normalized_BED_raw = [
    'rna_chr',
    'rna_bgn',
    'rna_end',
    'id',
    'rna_strand',
    'dna_chr',
    'dna_bgn',
    'dna_end',
    'dna_strand',
    'srr_id',
    'gene_name',
    'gene_type' ,
    'source',
    'bg_sm', 
    'N2_raw'
]

normalized_BED = [
    'rna_chr',
    'rna_bgn',
    'rna_end',
    'id',
    'rna_strand',
    'dna_chr',
    'dna_bgn',
    'dna_end',
    'dna_strand',
    'srr_id',
    'gene_name',
    'gene_type' ,
    'source',
    'bg_sm', 
    'N2'
]


normalized_BED_scaling = [
    'rna_chr',
    'rna_bgn',
    'rna_end',
    'id',
    'rna_strand',
    'dna_chr',
    'dna_bgn',
    'dna_end',
    'dna_strand',
    'srr_id',
    'gene_name',
    'gene_type' ,
    'source',
    'bg_sm',
    'N2',
    'SCnorm'
]


BED_Stereogene = [
    'chrom',
    'start',
    'end',
    'bg_sm'
]

BED3 = [
   'chrom',
   'start',
   'end',
]

BED6 = [
    'chrom',
    'start',
    'end',
    'name',
    'score',
    'strand'
]

chrsizes_format = [
    'chr',
    'len'
]

"""
File extensions for processing stages

"""

annot_suffix = '.0-corrected_annot.bed'
annot_bed6_suffix = ".0-corrected_annot.bed6"
annot_old_suffix = ".0-old_names.bed"
bad_annot_suffix = '.0-bad_annot.bed'
synonym_suffix = '.0-synonyms.tab'

CIGAR_suffix = '.1-CIGAR.tab'
CIGAR_suffix_stat = '.1-CIGAR.stat.tab'

blacklist_suffix = '.2-blacklist.tab'
blacklist_suffix_stat = '.2-blacklist.stat.tab'

singleton_suffix = '.3-singletons.tab'
no_annot_suffix = '.3-no_annot.tab'
complement_ann_suffix = '.3-complement_annot.tab'
selected_ann_suffix = '.3-selected_annot.tab'

voted_suffix = '.4-voted.tab'
voted_nr_suffix = '.4-voted_noribo.tab'
voted_stat_cnts_gene_type = '.4-counts_gene_type.stat.tab'
voted_stat_dens = '.4-genes_density.stat.tab'
voted_stat_cnts_by_genes = '.4-counts_by_genes.stat.tab'
voted_stat_cnts_by_chr = '.4-counts_by_chr.stat.tab'
voted_stat_cnts_by_source = '.4-counts_by_source.stat.tab'

N2_suffix = '.5-N2.tab'
N2_raw_suffix = '.5-N2_raw.tab'
bg_suffix_stat = '.5-background.stat.tab'
N2_raw_stats = '.5-N2_raw.stat.tab'
N2_final_stats='.5-N2.stat.tab'
Stereogene_sm_bg = '.5-background_sm.bgr'
Stereogene_bed = '.5-background.bed'

scaling_inner_suffix = '.6-SCinnerCHR.tab'
scaling_inter_suffix = '.6-SCinterCHR.tab'
PC_inner_suffix = '.6-PCinnerCHR.tab'
PC_inter_suffix = '.6-PCinterCHR.tab'

bardic_bg = '.7-background.BedGraph'
bardic_input_rdc = ".7-contacts_for_peaks.bed"
bardic_input_annot=".7-annot_for_peaks.bed"
bardic_bad_rdc = ".7-bad_contacts.bed"
bardic_bad_annot = ".7-bad_annot.bed"
