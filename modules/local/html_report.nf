process HTML_REPORT {
    conda "${projectDir}/envs/secondary_processing.yml"

    publishDir (
        path: "${params.outdir}/html_report",
        mode: "copy"
    )

    input:
    tuple val(meta), path(adaptersfastp), path(fastqc), path(debridged), path(restrsites) 

    output:
    path "${meta[0]}/", emit: folders

    script:
    """
    mkdir -p ${meta[0]}

    cp ${adaptersfastp} ${meta[0]}/
    cp -r ${fastqc} ${meta[0]}/
    cp ${debridged} ${meta[0]}/
    cp ${restrsites} ${meta[0]}/

    python3 ${projectDir}/bin/generate_html_report.py --sample-id ${meta[0]} --out-dir ${meta[0]} --output report_${meta[0]}.html
    """
}
