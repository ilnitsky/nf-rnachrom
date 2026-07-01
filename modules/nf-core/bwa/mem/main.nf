process BWA_MEM {
    tag "$meta.id"
    // label 'process_high'
    conda (params.use_nfcore_env ? "${moduleDir}/environment.yml" : "${projectDir}/envs/full_env.yml")
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:a34558545ae1413d94bde4578787ebef08027945-0' :
    //     'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:a34558545ae1413d94bde4578787ebef08027945-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = meta.method == 'RNA-seq' ? task.ext.args_rna : task.ext.args
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rna_prefix = meta.RNA
    def dna_prefix = meta.DNA
    def samtools_command = sort_bam ? 'sort' : 'view'
    // TO DO Fix iMARGI mapping
    if ( meta.method == "OTA" || meta.method == "RNA-seq" ){
        if (meta.single_end) {
            """
            INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

            bwa mem \\
            -a \\
            $args \\
            -t $task.cpus \\
            \$INDEX \\
            $reads \\
            2> sorted_${meta.id}_${prefix}.bwa.log \\
            | samtools  sort -n --threads $task.cpus -O BAM - > sorted_${meta.id}_${prefix}.bam

            ln -s sorted_${meta.id}_${prefix}.bam sorted_${meta.id}_${prefix}.COPY.bam

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
            """
        } else {
        """
            INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

            bwa mem \\
            -a \\
            $args \\
            -t $task.cpus \\
            \$INDEX \\
            $reads \\
            2> ${meta.id}.bwa.log \\
            | tee >(samtools view -@ ${task.cpus} -f 64  -b - | samtools sort -n -@ ${task.cpus} -o sorted_${meta.id}.r1.bam -) \\
            | samtools view -@ ${task.cpus} -f 128 -b - | samtools sort -n -@ ${task.cpus} -o sorted_${meta.id}.r2.bam -
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
            """
        }
    } else if (meta.method == 'ATA') {
        if (params.exp_type == 'imargi' || params.exp_type == 'margi') {
            // iMARGI: R1=RNA, R2=DNA are pre-split (circularized construct, no bridge).
            // Align each file independently as single-end; RNA_ALIGN and DNA_ALIGN call this separately.
            def imargi_args = meta.RNA ? (task.ext.args_rna ?: task.ext.args ?: '') : (task.ext.args_dna ?: task.ext.args ?: '')
            def imargi_prefix = meta.RNA ?: meta.DNA
            def imargi_postfix = meta.RNA ? 'rna' : 'dna'
            """
            INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

            bwa mem \\
                -a \\
                -5 \\
                -T 1 \\
                $imargi_args \\
                -t $task.cpus \\
                \$INDEX \\
                ${reads[0]} \\
                2> ${imargi_prefix}.bwa.log \\
                | samtools sort -n --threads $task.cpus -O BAM - > sorted_${meta.id}_${imargi_prefix}.${imargi_postfix}.bam

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
            """
        } else {
            args = meta.RNA ? (task.ext.args_rna ?: '') : (task.ext.args_dna ?: '')

            if (!task.ext.args_rna && meta.RNA) {
                log.warn "RNA aligner args not found, using empty string"
            }
            if (!task.ext.args_dna && !meta.RNA) {
                log.warn "DNA  aligner args not found, using empty string"
            }

            prefix = meta.RNA ? meta.RNA : meta.DNA
            def postfix = meta.RNA ? 'rna' : 'dna' 

        
            """
            INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

            bwa mem \\
                -a \\
                -T 10 \\
                $args \\
                -t $task.cpus \\
                \$INDEX \\
                $reads \\
                2> ${prefix}.bwa.log \\
                | samtools sort -n --threads $task.cpus -O BAM - > sorted_${meta.id}_${prefix}.${postfix}.bam

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
            """
        }
    
    }
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