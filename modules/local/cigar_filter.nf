
process CIGAR_FILTER {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$meta.prefix"

    publishDir (
        path: { "$params.outdir/Splicing_RNA" },
        mode: "copy"
    ) 

    input:
    tuple val(meta), path(beds)

    output:
    tuple val(meta), path("*-CIGAR.{tab,bed}"), emit: tab
    tuple val(meta), path("stat/*.stat"), emit: stat

    script:

    def rna_prefix = meta.RNA
    def dna_prefix = meta.DNA    
    def procedure = ''
    def separate_rna_dna = ''
    def cigar_filtered = rna_prefix + ".1-CIGAR.tab"
    procedure = params.procedure == 'new' ? '--rna_parts' : ''
    columns = params.procedure == 'new' ? '-f1-4,6-7' : '-f1-'
    separate_rna_dna = params.procedure == 'new' ? 'true' : ''
    // Add srr name of rna part to the last column of both files
    """
    cut ${columns} ${beds[0]} > ${beds[0]}.cols
    python3 ${projectDir}/bin/rnachrom_pipeline_faster/splicing.py ${beds[0]}.cols ${procedure} --outdir .
    mkdir stat
    wc -l ${cigar_filtered} | cut -f1 -d' ' > stat/${meta.prefix}.stat
    """


}

    // mkdir res
    // awk -v srr_id="${rna_prefix}" 'BEGIN{FS=OFS='\\t'} {print \$0, "\\t"srr_id}' ${cigar_filtered} > res/${rna_prefix}.tab
    // if [[ "${separate_rna_dna}" == "true" ]]; then
    //     awk -v srr_id="${rna_prefix}" 'BEGIN{FS=OFS='\\t'} {print \$0, "\\t"srr_id}' ${dna_prefix}.dna.bed > res/${dna_prefix}.tab
    // fi