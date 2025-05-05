
process ADD_SRR {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$meta.id $meta.prefix"

    publishDir (
        path: { "$params.outdir/Splicing_RNA" },
        mode: "copy"
    ) 

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*.{tab,bed}"), emit: tab

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
    awk -v srr_id="${rna_prefix}" 'BEGIN{FS=OFS="\\t"} NR==1{print \$0, OFS"SRR_ID"; next} {print \$0, OFS srr_id}' ${files[0]} > ${rna_prefix}.12_col.tab
    if [[ "${separate_rna_dna}" == "true" ]]; then
        awk -v srr_id="${dna_prefix}" 'BEGIN{FS=OFS="\\t"} NR==1{print \$0, OFS"SRR_ID"; next} {print \$0, OFS srr_id}' ${files[1]} > ${dna_prefix}.12_col.tab
    fi
    """


}


