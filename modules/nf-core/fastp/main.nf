process FASTP {
    tag "$meta.id,$meta.prefix"
    label 'process_medium'

    conda (params.use_nfcore_env ? "${moduleDir}/environment.yml" : "${projectDir}/envs/full_env.yml")
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
    //     'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)
    path(adapter_fasta)// val adapter_fasta
    val save_trimmed_fail
    val save_merged
    val only_remove_adapters

    output:
    tuple val(meta), path('*.{adapt,fastp}.fastq') , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    path "versions.yml"                       , emit: versions
    tuple val(meta), path('*.fail.fastq')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq'), optional:true, emit: reads_merged

    when:
    task.ext.when == null || task.ext.when

    script:
    
    // def adapter_fasta  = "${projectDir}/bin/adapters/TruSeq3-PE.fa"

    def args = only_remove_adapters ? (task.ext.args_adapters ?: '') : (task.ext.args ?: '')

    // def adapter_list = params.adapters_file ? "--adapter_fasta ${params.adapters_file}" : "--adapter_fasta ${projectDir}/assets/adapters_redc.fa"
    def detect_adapters = params.disable_adapter_autodetect ? "" : "--detect_adapter_for_pe" 
    def postfix = only_remove_adapters ? 'adapt' : 'fastp'
    def adapters = only_remove_adapters ? "${detect_adapters}" : "-A"

    def prefix = task.ext.prefix ?: "${meta.id}"
    def fail_fastq = save_trimmed_fail && meta.single_end ? "--failed_out ${prefix}.fail.fastq.gz" : save_trimmed_fail && !meta.single_end ? "--failed_out ${prefix}.paired.fail.fastq.gz --unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    // Use single ended for interleaved. Add --interleaved_in in config.

    //   | gzip -c > ${prefix}.fastp.fastq.gz

    if ( task.ext.args?.contains('--interleaved_in') ) {
        """
        [ ! -f  ${prefix}.fastq ] && ln -sf $reads ${prefix}.fastq

        fastp \\
            --stdout \\
            --in1 ${prefix}.fastq \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            --adapter_fasta ${adapter_fasta}  \\
            $fail_fastq \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2) \\
        > ${prefix}.fastp.fastq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq ] && ln -sf $reads ${prefix}.fastq

        fastp \\
            --in1 ${prefix}.fastq \\
            --out1  ${prefix}.${postfix}.fastq \\
            --thread $task.cpus \\
            --json ${prefix}.${postfix}.json \\
            --html ${prefix}.${postfix}.html \\
            --adapter_fasta ${adapter_fasta}  \\
            $adapters \\
            $fail_fastq \\
            $args \\
            2> >(tee ${prefix}.${postfix}.log >&2)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq" : ''

        def input1  = (params.bridge_processing || meta.method == "OTA" || meta.method == "RNA-seq") ? "${prefix}_1.fastq" : "${meta.RNA}.fastq"
        def input2  = (params.bridge_processing || meta.method == "OTA" || meta.method == "RNA-seq") ? "${prefix}_2.fastq" : "${meta.DNA}.fastq"
        def output1 = (params.bridge_processing || meta.method == "OTA" || meta.method == "RNA-seq") ? "${prefix}_1.${postfix}.fastq" : "1_${meta.RNA}.${postfix}.fastq"
        def output2 = (params.bridge_processing || meta.method == "OTA" || meta.method == "RNA-seq") ? "${prefix}_2.${postfix}.fastq" : "2_${meta.DNA}.${postfix}.fastq"

        """
        [ ! -f  ${input1} ] && ln -sf ${reads[0]} ${input1}
        [ ! -f  ${input2} ] && ln -sf ${reads[1]} ${input2}
        fastp \\
            --in1 ${input1} \\
            --in2 ${input2} \\
            --out1 ${output1} \\
            --out2 ${output2} \\
            --json ${prefix}.${postfix}.json \\
            --html ${prefix}.${postfix}.html \\
            $adapters \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }

    stub:
    def prefix              = task.ext.prefix ?: "${meta.id}"
    def is_single_output    = task.ext.args?.contains('--interleaved_in') || meta.single_end
    def touch_reads         = is_single_output ? "${prefix}.fastp.fastq.gz" : "${prefix}_1.fastp.fastq.gz ${prefix}_2.fastp.fastq.gz"
    def touch_merged        = (!is_single_output && save_merged) ? "touch ${prefix}.merged.fastq.gz" : ""
    """
    touch $touch_reads
    touch "${prefix}.fastp.json"
    touch "${prefix}.fastp.html"
    touch "${prefix}.fastp.log"
    $touch_merged

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
