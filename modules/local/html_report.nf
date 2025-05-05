process HTML_REPORT {
    conda "${projectDir}/envs/full_env.yml"
    // conda "${projectDir}/envs/secondary_processing.yml"
  
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'ilnitsky/nf-rnachrom:latest' : '' }"
        
   
    publishDir (
        path: "${params.outdir}/html_report",
        mode: "copy"
    )

    input:
    val(sample_paths)
    
    output:
    path "*.zip", emit: zip

    script:
    """
    python3 ${projectDir}/bin/generate_html_report.py --sample-paths ${sample_paths.join(' ')} --output-dir ./
    """
}
// python3 ${projectDir}/bin/generate_html_report.py --sample-id ${meta[0]} --out-dir ${meta[0]} --output ${meta[0]}/report_${meta[0]}.html
        // zip -r ${meta[0]}.zip ${meta[0]}/