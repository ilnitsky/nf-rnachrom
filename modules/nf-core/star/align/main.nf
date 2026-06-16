process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.use_nfcore_env ? "${moduleDir}/environment.yml" : "${projectDir}/envs/full_env.yml")
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' :
    //     'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' }"

    input:
    // tuple val(meta), path(reads, stageAs: "input*/*")
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    val star_ignore_sjdbgtf
    val seq_platform
    val seq_center

    output:
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*d.out.bam')              , optional:true, emit: bam
    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab
    tuple val(meta), path('*.SJ.out.tab')            , optional:true, emit: spl_junc_tab
    tuple val(meta), path('*.ReadsPerGene.out.tab')  , optional:true, emit: read_per_gene_tab
    tuple val(meta), path('*.out.junction')          , optional:true, emit: junction
    tuple val(meta), path('*.out.sam')               , optional:true, emit: sam
    tuple val(meta), path('*.wig')                   , optional:true, emit: wig
    tuple val(meta), path('*.bg')                    , optional:true, emit: bedgraph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ignore_gtf   = star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    def seq_platform = seq_platform ? "'PL:$seq_platform'" : ""
    def seq_center   = seq_center   ? "'CN:$seq_center'"   : ""
    def attrRG       = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:$prefix' $seq_center 'SM:$prefix' $seq_platform"
    def read_files_command = reads[0].toString().endsWith('.gz') ? '--readFilesCommand zcat' : ''

    if (meta.method == "OTA" || meta.method == "RNA-seq") {
        if (meta.single_end) {
            """
            STAR \\
                --genomeDir $index \\
                --readFilesIn ${reads[0]} \\
                --runThreadN $task.cpus \\
                --outFileNamePrefix ${prefix}. \\
                $read_files_command \\
                $ignore_gtf \\
                $attrRG \\
                $args

            ln -s ${prefix}.Aligned.sortedByCoord.out.bam ${prefix}.COPY.bam || true

            if [ -f ${prefix}.Unmapped.out.mate1 ]; then
                mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
                gzip ${prefix}.unmapped_1.fastq
            fi

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                star: \$(STAR --version | sed -e "s/STAR_//g")
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
                gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
            END_VERSIONS
            """
        } else {
            """
            STAR \\
                --genomeDir $index \\
                --readFilesIn ${reads[0]} ${reads[1]} \\
                --runThreadN $task.cpus \\
                --outFileNamePrefix ${prefix}. \\
                $read_files_command \\
                $ignore_gtf \\
                $attrRG \\
                $args \\
                | tee >(samtools view -@ ${task.cpus} -f 64  -b - | samtools sort -n -@ ${task.cpus} -o sorted_${meta.id}.r1.bam -) \\
                | samtools view -@ ${task.cpus} -f 128 -b - | samtools sort -n -@ ${task.cpus} -o sorted_${meta.id}.r2.bam -

            if [ -f ${prefix}.Unmapped.out.mate1 ]; then
                mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
                gzip ${prefix}.unmapped_1.fastq
            fi
            if [ -f ${prefix}.Unmapped.out.mate2 ]; then
                mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
                gzip ${prefix}.unmapped_2.fastq
            fi

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                star: \$(STAR --version | sed -e "s/STAR_//g")
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
                gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
            END_VERSIONS
            """
        }
    } else if (meta.method == "ATA") {
        args = meta.rna ? (task.ext.args_rna ?: task.ext.args ?: '') : (task.ext.args_dna ?: task.ext.args ?: '')

        if (!task.ext.args_rna && meta.rna) {
            log.warn "RNA aligner args not found, using empty string"
        }
        if (!task.ext.args_dna && !meta.rna) {
            log.warn "DNA aligner args not found, using empty string"
        }

        prefix = meta.RNA ? meta.RNA : meta.DNA
        def postfix = meta.RNA ? 'rna' : 'dna'

        """
        STAR \\
            --genomeDir $index \\
            --readFilesIn ${reads[0]} \\
            --runThreadN $task.cpus \\
            --outFileNamePrefix ${prefix}. \\
            --outSAMtype BAM SortedByCoordinate \\
            # --outSAMtype BAM Unsorted \\
            # --outStd BAM_Unsorted \\
            $read_files_command \\
            $ignore_gtf \\
            $attrRG \\
            $args 
            # | samtools sort -n --threads $task.cpus -O BAM - > sorted_${meta.id}_${prefix}.${postfix}.bam

        if [ -f ${prefix}.Unmapped.out.mate1 ]; then
            mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
            gzip ${prefix}.unmapped_1.fastq
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}Xd.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.sortedByCoord.out.bam
    touch ${prefix}.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.unmapped_1.fastq.gz
    touch ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}.tab
    touch ${prefix}.SJ.out.tab
    touch ${prefix}.ReadsPerGene.out.tab
    touch ${prefix}.Chimeric.out.junction
    touch ${prefix}.out.sam
    touch ${prefix}.Signal.UniqueMultiple.str1.out.wig
    touch ${prefix}.Signal.UniqueMultiple.str1.out.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}
