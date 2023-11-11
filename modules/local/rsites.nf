
// conda "conda-forge::gawk=5.1.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
    //     'quay.io/biocontainers/gawk:5.1.0' }"

process RSITES {
    // publishDir ${params.outdir}/'processing', mode: 'copy'
    tag "$config"
    label 'process_single'

    
    input:
    path samplesheet
    path config
    path reads

    output:
    path '*.final.json', emit: json
    path 'rsites/*', emit: rsites


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    
    jq -nR '[inputs|split(",")|.[1]|split("/")|.[-1]|split(".")|.[0]]' \\
    <(tail -n+2 $samplesheet) | jq --argjson ids "\$(< $config)" '\$ids + {rna_ids: .}' \\
    > config.json

    jq -nR '[inputs|split(",")|.[2]|split("/")|.[-1]|split(".")|.[0]]' \\
    <(tail -n+2 $samplesheet) | jq --argjson ids "\$(< $config)" '\$ids + {dna_ids: .}' \\
    > config.json

    jq '.hisat.genome="${params.genome_path}"' config.json \\
    | jq '.rsites.type="${params.exp_type}"' \\
    | jq '.trim.tool="${params.trim_tool}"' \\
    | jq '.trim.tool_path="${params.trim_path}"' \\
    > config.final.json

    rnachromprocessing -c config.final.json -s rsites -v
    """

}





