process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "${projectDir}/envs/full_env.yml"
    // conda "${projectDir}/envs/secondary_processing.yml"
  
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'ilnitsky/nf-rnachrom:latest' : '' }"
        
   
    // conda "conda-forge::python=3.8.3"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/python:3.8.3' :
    //     'biocontainers/python:3.8.3' }"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnachrom/bin/

    if( params.bridge_processing || params.exp_type in ['chart', 'rap', 'chirp']) {
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed \'s/Python //g\')
    END_VERSIONS
    """
    } else {
    """
    check_samplesheet_rna_dna_parts.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed \'s/Python //g\')
    END_VERSIONS
    """
    }
}
