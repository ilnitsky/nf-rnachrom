
// conda "conda-forge::gawk=5.1.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
    //     'quay.io/biocontainers/gawk:5.1.0' }"

process GET_RAW_CONTACTS_A2A {
    tag "$config"
    label 'process_single'
    
    input:
    path samplesheet
    path config

    output:
    path '*', emit: json

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    jq -nR '[inputs|split(",")|.[2]]' \\
    <(tail -n+2 $samplesheet) | jq --argjson ids "\$(< $config)" '\$ids + {rna_ids: .}' \\
    > ${workflow.projectDir}/assets/imargi_test.json 

    """

}





