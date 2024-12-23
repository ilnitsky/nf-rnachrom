process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.paired.trim*.fastq")   , emit: trimmed_reads
    tuple val(meta), path("*.unpaired.trim_*.fastq"), optional:true, emit: unpaired_reads
    tuple val(meta), path("*.log")                     , emit: log
    tuple val(meta), path("*.summary")                 , emit: summary
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def adapters = params.adapters_file ? "${params.adapters_file}" : ""
    def adapters = params.adapters_file ?  "${params.adapters_file}" : (meta.single_end ? "TruSeq3-SE" : "TruSeq3-PE.fa")
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "SE" : "PE"
    def output = meta.single_end ?
        "${prefix}.SE.paired.trim.fastq" // HACK to avoid unpaired and paired in the trimmed_reads output
        : "${prefix}.paired.trim_1.fastq ${prefix}.unpaired.trim_1.fastq ${prefix}.paired.trim_2.fastq ${prefix}.unpaired.trim_2.fastq"
    // TODO Give better error output
    def qual_trim = task.ext.args2 ?: ''
    """
    trimmomatic \\
        $trimmed \\
        -threads $task.cpus \\
        -trimlog ${prefix}.log \\
        -summary ${prefix}.summary \\
        $reads \\
        $output \\
        $qual_trim \\
        ILLUMINACLIP:$adapters:2:30:10 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}
