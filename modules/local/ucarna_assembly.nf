
// Assembles ucaRNAs with StringTie from already-mapped RNA-part BAMs and scores
// them with a Poisson p-value. One process call handles one biological
// replicate: bams/id_lists/vote_files are the technical replicates belonging
// to it (see final_rc.nf, grouped by the same {id, rnaseq} key MERGE_REPLICAS
// uses), merged into one bam by ucarna_assembly.sh before assembly. Each
// replicate's DETECT_STRAND *_wins.tsv vote is passed through so
// ucarna_assembly.sh can correct FLAG-derived strand (ANTI -> flip FLAG 0x10)
// before any strand-aware step - see ucarna_assembly.sh header for why that's
// safe here and DETECT_STRAND.out.strand_vote_result for the vote format.
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
    // bams/id_lists/vote_files must be given in the same order.
    tuple val(meta), path(bams), path(id_lists), path(vote_files)
    path(annot_gtf)

    output:
    tuple val(meta), path('*.ucaRNAs.gtf'),   emit: gtf
    tuple val(meta), path('*.ucaRNAs.bedrc'), emit: bedrc
    tuple val(meta), path('*.ucaRNAs.tab'),   emit: table
    tuple val(meta), path('*.ucaRNAs.pdf'),   emit: pdf

    script:
    def prefix = params.ucarna_prefix ?: 'uca'
    def suffix_arg = params.ucarna_suffix ? "-x ${params.ucarna_suffix}" : ''
    def bam_args  = bams.collect { "-b ${it}" }.join(' ')
    def id_args   = id_lists.collect { "-i ${it}" }.join(' ')
    def vote_args = vote_files.collect { "-v ${it}" }.join(' ')

    """
    ucarna_assembly.sh \\
        -g ${annot_gtf} \\
        ${bam_args} \\
        ${id_args} \\
        ${vote_args} \\
        -p ${prefix} \\
        ${suffix_arg} \\
        -o . \\
        -t ${task.cpus}
    """
}
