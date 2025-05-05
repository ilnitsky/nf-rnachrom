process HISAT2_ALIGN {
    tag "$meta.id"
    // label 'process_high'
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.use_nfcore_env ? "${moduleDir}/environment.yml" : "${projectDir}/envs/full_env.yml")
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'ilnitsky/nf-rnachrom:latest' : '' }"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
    //     'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(splicesites)

    output:
    tuple val(meta), path("*.bam")                   , emit: bam
    tuple val(meta), path("*.log")                   , emit: summary
    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def samtools_view_args = '-F 4 -F 8 -F 256'
    def samtools_view_args = '-f 2'  // read is properly paired
    def VERSION = '2.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = meta.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (meta.strandedness == 'reverse') {
        strandedness = meta.single_end ? '--rna-strandness R' : '--rna-strandness RF'
    }
    ss = "$splicesites" ? "--known-splicesite-infile $splicesites" : ''
    // def seq_center = params.seq_center ? "--rg-id ${prefix} --rg SM:$prefix --rg CN:${params.seq_center.replaceAll('\\s','_')}" : "--rg-id ${prefix} --rg SM:$prefix"
    // if (meta.method == "OTA" || meta.method == "RNA-seq" ){    
    if ( meta.method == "OTA" || meta.method == "RNA-seq" ){
        if (meta.single_end) {
            // def unaligned = params.save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
            """
            INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`
            hisat2 \\
                -x \$INDEX \\
                -U $reads \\
                $strandedness \\
                $ss \\
                --summary-file ${prefix}.hisat2.summary.log \\
                --threads $task.cpus \\
                $args \\
                | samtools view -bS - > ${meta.id}_${prefix}.bam

            ln -s ${meta.id}_${prefix}.bam ${meta.id}_${prefix}.COPY.bam 
            
            if [ -f *.tmp.* ]; then
                rm *.tmp.*
            fi
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                hisat2: $VERSION
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
            """
        } else {
            // def unaligned = params.save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
            """
            INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`
            hisat2 \\
                -x \$INDEX \\
                -1 ${reads[0]} \\
                -2 ${reads[1]} \\
                $strandedness \\
                $ss \\
                --summary-file ${prefix}.hisat2.summary.log \\
                --threads $task.cpus \\
                --no-mixed \\
                --no-discordant \\
                $args \\
                | tee >(samtools view -bS -f 64 ${samtools_view_args} - > ${meta.id}_${prefix}.R1.bam) \\
                | samtools view -bS -f 128 ${samtools_view_args} - > ${meta.id}_${prefix}.R2.bam

            if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
                mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
            fi
            if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
                mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
            fi

            if [ -f *.tmp.* ]; then
                rm *.tmp.*
            fi

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                hisat2: $VERSION
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
            """
        }
    } else if (meta.method == "ATA") {
        //TODO fix dna/rna determination for alignment
            args = meta.rna ? (task.ext.args_rna ?: '') : (task.ext.args_dna ?: '')

            //  logging
            if (!task.ext.args_rna && meta.rna) {
                log.warn "RNA aligner args not found, using empty string"
            }
            if (!task.ext.args_dna && !meta.rna) {
                log.warn "DNA  aligner args not found, using empty string"
            }

            // def args_rna = task.ext.args_rna ?: ''
            // def args = meta.rna ? task.ext.args_rna : task.ext.args_dna 
            // def rna_prefix = meta.RNA
            // def dna_prefix = meta.DNA
            prefix = meta.RNA ? meta.RNA : meta.DNA
            def postfix = meta.RNA ? 'rna' : 'dna' 


            """
            INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`

            hisat2 \\
                -x \$INDEX \\
                -U ${reads[0]} \\
                $ss \\
                --summary-file ${prefix}.hisat2.summary.log \\
                --threads $task.cpus \\
                $args \\
                $strandedness \\
                | samtools sort -n --threads $task.cpus -O BAM - > sorted_${meta.id}_${prefix}.${postfix}.bam

            if [ -f *.tmp.* ]; then
                rm *.tmp.*
            fi

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                hisat2: $VERSION
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
            """
    }
}




// } else if (meta.method == "ATA") {
//         //TODO fix dna/rna determination for alignment
//             def rna_prefix = meta.RNA
//             def dna_prefix = meta.DNA

//             """
//             INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`

//             hisat2 \\
//                 -x \$INDEX \\
//                 -U ${reads[0]} \\
//                 $ss \\
//                 --summary-file ${rna_prefix}.hisat2.summary.log \\
//                 --threads $task.cpus \\
//                 $args_rna \\
//                 $strandedness \\
//                 | samtools view -bSh - > ${meta.id}_${rna_prefix}.rna.bam

//             hisat2 \\
//                 -x \$INDEX \\
//                 -U ${reads[1]} \\
//                 --summary-file ${dna_prefix}.hisat2.summary.log \\
//                 --threads $task.cpus \\
//                 $args_dna \\
//                 $strandedness \\
//                 | samtools view -bSh - > ${meta.id}_${dna_prefix}.dna.bam

//             cat <<-END_VERSIONS > versions.yml
//             "${task.process}":
//                 hisat2: $VERSION
//                 samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
//             END_VERSIONS
//             """
//     }
