
process STAR_CHIMERIC_READS {
    conda "bioconda::samtools=1.19.2 bioconda::pysam"
    tag "$meta.id $meta.prefix"

    publishDir (
        path: { "$params.outdir/Splicing_RNA" },
        mode: "copy"
    ) 

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.star_chimeric.bam"), emit: bam

    script:
    def chimSegmentMin = 13
    """
    samtools index ${bams[0]}
    samtools index ${bams[1]}

    python ${projectDir}/bin/star_bam_remap.py ${bams[0]} ${meta.RNA}.star_chimeric.bam ${chimSegmentMin} > ${meta.RNA}.star_bam_remap.log
    python ${projectDir}/bin/star_bam_remap.py ${bams[1]} ${meta.DNA}.star_chimeric.bam ${chimSegmentMin} > ${meta.DNA}.star_bam_remap.log
    """


}


