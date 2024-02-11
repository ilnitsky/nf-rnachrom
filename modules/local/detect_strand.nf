process MAPPING {
    tag "$config"
    label 'process_high'

    conda '/home/ryabykhgrigory/gpfs/nf-core-rnachrom/envs/rnachromprocessing.yaml'

    input:
    path config
    path contacts

    output:
    path 'contacts/*', emit: contacts
    stdout emit: stdout_ch


    when:
    task.ext.when == null || task.ext.when

    script:
    """

    jq '.hisat.known_splice="${params.splice_sites}"' $config \\
    | jq '.rsites.type="skip"' \\
    | jq '.trim.tool="skip"' \\
    > config.map.json

    rnachromprocessing -c config.final.json -v
    """

}