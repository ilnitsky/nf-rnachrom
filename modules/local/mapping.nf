process MAPPING {
    // publishDir ${params.outdir}/'processing', mode: 'copy'
    tag "$config"
    label 'process_high'

    conda '/home/ryabykhgrigory/gpfs/nf-core-rnachrom/envs/rnachromprocessing.yaml'

    input:
    path config
    path reads

    output:
    path 'contacts/*', emit: contacts
    stdout emit: stdout_ch


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pip install ${projectDir}/bin/rnachromprocessing_allvsall/
    which hisat2
    hisat2 --version
    jq '.hisat.known_splice="${params.splice_sites}"' $config \\
    | jq '.rsites.type="skip"' \\
    | jq '.trim.tool="skip"' \\
    > config.map.json

    rnachromprocessing -c config.final.json -v
    """

}





