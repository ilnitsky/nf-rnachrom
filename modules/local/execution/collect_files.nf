process COLLECT_FILES {
    
    conda "${projectDir}/envs/full_env.yml"
    
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'ilnitsky/nf-rnachrom:latest' : '' }"
        
    publishDir (
        path: "${params.outdir}/html_report",
        mode: "copy"
    )

    input:
    tuple val(meta), path(adaptersfastp), path(fastqc_initial), path(trim_log), path(fastqc_final), 
            path(debridged), path(restrsites),  path(rna_align_log), path(dna_align_log), path(filter_png), path(detect_strand_png)


    output:
    path "${meta[0]}-${meta[1]}/", emit: folders

    script:
    """
    mkdir -p ${meta[0]}-${meta[1]}


    if [ -f "${adaptersfastp}" ] && [ "${adaptersfastp}" != "NO_FILE" ]; then
        cp ${adaptersfastp} ${meta[0]}-${meta[1]}/
    fi
    
    for fqc_file in ${fastqc_initial}; do
        if [ -f "\$fqc_file" ] && [ "\$fqc_file" != "NO_FILE" ]; then
            cp \$fqc_file ${meta[0]}-${meta[1]}/
        fi
    done

    for fqc_file in ${fastqc_final}; do
        if [ -f "\$fqc_file" ] && [ "\$fqc_file" != "NO_FILE" ]; then
            cp \$fqc_file ${meta[0]}-${meta[1]}/
        fi
    done
    
    for deb_file in ${debridged}; do
        if [ -f "\$deb_file" ] && [ "\$deb_file" != "NO_FILE" ]; then
            cp \$deb_file ${meta[0]}-${meta[1]}/
        fi
    done
    
    for restr_file in ${restrsites}; do
        if [ -f "\$restr_file" ] && [ "\$restr_file" != "NO_FILE" ]; then
            cp \$restr_file ${meta[0]}-${meta[1]}/
        fi
    done
    
    for trim_log_file in ${trim_log}; do
        if [ -f "\$trim_log_file" ] && [ "\$trim_log_file" != "NO_FILE" ]; then
            cp \$trim_log_file ${meta[0]}-${meta[1]}/
        fi
    done

    for rna_align_log_file in ${rna_align_log}; do
        if [ -f "\$rna_align_log_file" ] && [ "\$rna_align_log_file" != "NO_FILE" ]; then
            cp \$rna_align_log_file ${meta[0]}-${meta[1]}/
        fi
    done

    for dna_align_log_file in ${dna_align_log}; do
        if [ -f "\$dna_align_log_file" ] && [ "\$dna_align_log_file" != "NO_FILE" ]; then
            cp \$dna_align_log_file ${meta[0]}-${meta[1]}/
        fi
    done

    for filter_png_file in ${filter_png}; do
        if [ -f "\$filter_png_file" ] && [ "\$filter_png_file" != "NO_FILE" ]; then
            cp \$filter_png_file ${meta[0]}-${meta[1]}/
        fi
    done

    for detect_strand_png_file in ${detect_strand_png}; do
        if [ -f "\$detect_strand_png_file" ] && [ "\$detect_strand_png_file" != "NO_FILE" ]; then
            cp \$detect_strand_png_file ${meta[0]}-${meta[1]}/
        fi
    done

    """
}
// python3 ${projectDir}/bin/generate_html_report.py --sample-id ${meta[0]} --out-dir ${meta[0]} --output ${meta[0]}/report_${meta[0]}.html
        // zip -r ${meta[0]}.zip ${meta[0]}/