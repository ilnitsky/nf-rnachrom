
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
    tuple val(meta), path("res/*.{tab,bed}"), emit: tab
    tuple val(meta), path("*.stat"), emit: stat

    script:
    def procedure = ''
    def separate_rna_dna = ''
    def cigar_filtered = meta.prefix + "_1.1-CIGAR.tab"
    procedure = params.procedure == 'new' ? '--rna_parts' : ''
    separate_rna_dna = params.procedure == 'new' ? 'true' : ''
    
    """
    cut -f1-4,6-7 ${beds[0]} > ${beds[0]}.cols
    python3 ${projectDir}/bin/rnachrom_pipeline_faster/splicing.py ${beds[0]}.cols ${procedure} --outdir .
    wc -l ${cigar_filtered} | cut -f1 -d' ' > ${meta.prefix}.stat

    mkdir res
    awk -v srr_id="${meta.prefix}" 'BEGIN{FS=OFS='\\t'} {print \$0, "\\t"srr_id}' ${meta.prefix}_1.1-CIGAR.tab > res/${meta.prefix}_1.CIGAR.tab
    if [[ "${separate_rna_dna}" == "true" ]]; then
        awk -v srr_id="${meta.prefix}" 'BEGIN{FS=OFS='\\t'} {print \$0, "\\t"srr_id}' ${meta.prefix}_2.bed > res/${meta.prefix}_2.tab
    fi
    """


}