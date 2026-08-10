process CHROMATIN_POTENTIAL {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    conda "${projectDir}/envs/full_env.yml"


    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

    publishDir (
        path: { "$params.outdir/chromatin_potential" },
        mode: "copy"
    )

    input:
    tuple val(meta), path(uu_contacts), path(um_contacts)
    path(rnaseq_data)
    path chrom_sizes

    output:
    // Original declaration named files (*.normalized.tab / *.stats.txt) that
    // count_contacts_all.sh/RD_chP.py never actually produce -- real output is
    // 4 result tables (one per UU/UU_UM x all/filter_dist combination) plus
    // matching plots, all under ./chP/, named chP_<type>.tab / chP_<type>.png.
    // tuple val(meta), path("*.normalized.tab"), emit: normalized_contacts
    // tuple val(meta), path("*.stats.txt"), emit: stats
    tuple val(meta), path("chP/chP_*.tab"), emit: chp_tables
    tuple val(meta), path("chP/chP_*.png"), emit: chp_plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    
    count_contacts_all.sh \
            -d 500000 \
            -i . \
            -o ./chP \
            -u ${uu_contacts} \
            -m ${um_contacts} \
            -n 100 \
            -f 0.05 \
            -r ${rnaseq_data} \
            -s ${projectDir}/bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}

    // count_contacts_all.sh \
    //         -d 500000 \
    //         -i . \
    //         -o ./chP \
    //         -u contacts.voting.UU.bed \
    //         -m contacts.voting.UM.bed \
    //         -n 100 \
    //         -f 0.05 \
    //         -r ${rnaseq_data}
    
