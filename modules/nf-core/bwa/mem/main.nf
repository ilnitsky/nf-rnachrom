process BWA_MEM {
    tag "$meta.id"
    // label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:a34558545ae1413d94bde4578787ebef08027945-0' :
        'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:a34558545ae1413d94bde4578787ebef08027945-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rna_prefix = meta.RNA
    def dna_prefix = meta.DNA
    def samtools_command = sort_bam ? 'sort' : 'view'
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | tee >(samtools view $args2 -f 64 --threads $task.cpus -o ${rna_prefix}.rna.bam) \\
              >(samtools view $args2 -f 128 --threads $task.cpus -o ${dna_prefix}.dna.bam) > /dev/null 



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
        // | samtools $samtools_command $args2 -h --threads $task.cpus - \\

    //         bwa mem \\
    //     $args \\
    //     -t $task.cpus \\
    //     \$INDEX \\
    //     ${reads[0]} \\
    //     | samtools view -f 256 --threads $task.cpus -o ${prefix}.rna.bam

    // bwa mem \\
    //     $args \\
    //     -t $task.cpus \\
    //     \$INDEX \\
    //     ${reads[1]} \\
    //     | samtools view -f 256 --threads $task.cpus -o ${prefix}.rna.bam