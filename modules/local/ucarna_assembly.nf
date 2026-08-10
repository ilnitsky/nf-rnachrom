
// Assembles ucaRNAs with StringTie from already-mapped, already strand-corrected
// RNA-part BAMs and scores them with a Poisson p-value. One process call handles
// one biological replicate: bams/id_lists are the technical replicates belonging
// to it (see final_rc.nf, grouped by the same {id, rnaseq} key MERGE_REPLICAS
// uses), merged into one bam by ucarna_assembly.sh before assembly. Strand flips
// flagged by DETECT_STRAND are applied by hand to this process's output
// afterwards, not inside it (see ucarna_assembly.sh header).
//
// Wired into final_rc.nf's ATA workflow only (final_ota.nf is untouched),
// gated behind params.run_ucarna_assembly.

process UCARNA_ASSEMBLY {
    conda "${projectDir}/envs/full_env.yml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

    label 'process_medium'
    publishDir (
        path: { "$params.outdir/ucaRNA_assembly/${meta.id}" },
        mode: "copy"
    )

    input:
    // meta.id groups technical replicates into one biological replicate assembly;
    // bams/id_lists must be given in the same order.
    tuple val(meta), path(bams), path(id_lists)
    path(annot_gtf)

    output:
    tuple val(meta), path('*.ucaRNAs.gtf'),   emit: gtf
    tuple val(meta), path('*.ucaRNAs.bedrc'), emit: bedrc
    tuple val(meta), path('*.ucaRNAs.tab'),   emit: table
    tuple val(meta), path('*.ucaRNAs.pdf'),   emit: pdf

    script:
    def prefix = params.ucarna_prefix ?: 'uca'
    def suffix_arg = params.ucarna_suffix ? "-x ${params.ucarna_suffix}" : ''
    def bam_args = bams.collect { "-b ${it}" }.join(' ')
    def id_args  = id_lists.collect { "-i ${it}" }.join(' ')

    """
    ucarna_assembly.sh \\
        -g ${annot_gtf} \\
        ${bam_args} \\
        ${id_args} \\
        -p ${prefix} \\
        ${suffix_arg} \\
        -o . \\
        -t ${task.cpus}
    """
}
