
process PAIR_AND_SORT_BAMS {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$meta.id $meta.prefix"

    publishDir (
        path: { "$params.outdir/Splicing_RNA" },
        mode: "copy"
    ) 

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.{tab,bed}"), emit: tab

    script:

    def rna_prefix = meta.RNA
    def dna_prefix = meta.DNA    

    """
    samtools merge -nr -f  ${rna_prefix}.merge.bam  ${rna_prefix}.rna.bam ${dna_prefix}.dna.bam 
    samtools sort -n ${rna_prefix}.merge.bam > ${rna_prefix}.star_sorted.bam
    ${projectDir}/bin/bam_merge_v0.3 ${rna_prefix}.star_sorted.bam ${rna_prefix}.star_bwa.bam
    rm  ${rna_prefix}.merge.bam ${rna_prefix}.star_sorted.bam
    samtools sort -n  ${rna_prefix}.star_bwa.bam | samtools fixmate -  ${rna_prefix}.fixmate.star.bam
    """


}


