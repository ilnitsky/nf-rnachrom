include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/gunzip/main'
include { HISAT2_ALIGN                  } from '../../modules/nf-core/hisat2/align/main'  
include { HISAT2_EXTRACTSPLICESITES     } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                  } from '../../modules/nf-core/hisat2/build'
include { STAR_ALIGN                    } from '../../modules/nf-core/star/align/main'
include { STAR_GENOMEGENERATE           } from '../../modules/nf-core/star/genomegenerate' 
include { STAR_CHIMERIC_READS           } from '../../modules/local/star_chimeric'    
include { BOWTIE2_ALIGN                 } from '../../modules/nf-core/bowtie2/align/main' 
include { BOWTIE2_BUILD                 } from '../../modules/nf-core/bowtie2/build/main'
include { SAM_TO_FASTQ                  } from '../../modules/local/sam_to_fastq'
include { BWA_MEM                       } from '../../modules/nf-core/bwa/mem/main' 
include { BWA_INDEX                     } from '../../modules/nf-core/bwa/index/main'  


workflow ALIGN {
    take:
    ch_input_align // file: /path/to/samplesheet.csv
    ch_splicesites

    main:
    ch_versions     = Channel.empty()
    ch_hisat2_index   = params.hisat2_index ? Channel.fromPath(params.hisat2_index) : Channel.empty()
    ch_bwa_sort     = Channel.empty()
    ch_bwa_index    = params.bwa_index ? Channel.fromPath(params.bwa_index) : Channel.empty()
    ch_align_log    = Channel.empty()
    ch_gtf = Channel.value(params.annot_GTF)

    //
    // generate HISAT2 index from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()

    if (params.genome_fasta.endsWith('.gz')) {
        ch_genome_fasta    = GUNZIP_FASTA ( [ [:], params.genome_fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_genome_fasta = Channel.value(params.genome_fasta)
    }

    if (params.align_tool == 'hisat2') {
        if (params.splice_sites == null) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            ch_splicesites = Channel.fromPath(params.splice_sites , checkIfExists: true)
        }
        if (params.hisat2_index != null) {
                ch_hisat2_index = Channel.fromPath(params.hisat2_index , checkIfExists: true )
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_genome_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }

        HISAT2_ALIGN( 
            ch_input_align,
            ch_hisat2_index.map { [ [:], it ] }.collect(),
            ch_splicesites.map { [ [:], it ] }.collect()
        )
        ch_align_bam       = HISAT2_ALIGN.out.bam
        ch_align_log       = HISAT2_ALIGN.out.summary
        ch_versions        = ch_versions.mix(HISAT2_ALIGN.out.versions)
    }


    if (params.align_tool == 'star') {
        if (params.star_index != null) {
                ch_star_index = Channel.fromPath(params.star_index, checkIfExists: true)
        } else {
            ch_star_index = STAR_GENOMEGENERATE ( ch_genome_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }
        STAR_ALIGN ( 
            ch_input_align,
            ch_star_index.map { [ [:], it ] }.collect(), 
            ch_gtf.map { [ [:], it ] }.collect(), 
            false, 
            false, 
            false
        )
        ch_align_bam      = STAR_ALIGN.out.bam
        ch_log_final      = STAR_ALIGN.out.log_final
        ch_log_out        = STAR_ALIGN.out.log_out
        ch_log_progress   = STAR_ALIGN.out.log_progress
        ch_bam_sorted     = STAR_ALIGN.out.bam_sorted
        ch_bam_transcript = STAR_ALIGN.out.bam_transcript
        ch_fastq          = STAR_ALIGN.out.fastq
        ch_tab            = STAR_ALIGN.out.tab
        ch_versions       = ch_versions.mix(STAR_ALIGN.out.versions.first())


        ch_align_bam.view()

        // if (params.exp_type == 'imargi') {
        //     STAR_CHIMERIC_READS ( ch_align_bam )
        //     ch_align_bam = STAR_CHIMERIC_READS.out.bam
        // }
    }

    if (params.align_tool == 'bowtie2') {

        if (params.bowtie2_index != null) {
                ch_bowtie2_index = Channel.fromPath(params.bowtie2_index , checkIfExists: true )
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_genome_fasta.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions     = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
        ch_bowtie2_index.view()
        BOWTIE2_ALIGN( 
            ch_input_align,
            ch_bowtie2_index.map { [ [:], it ] }.collect(),
            ch_genome_fasta.map { [ [:], it ] }.collect(),
            false,
            false
        )
        ch_align_bam       = BOWTIE2_ALIGN.out.bam
        ch_align_log       = BOWTIE2_ALIGN.out.log
        ch_versions        = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
    }

    if (params.align_tool == 'bwa_mem') {

        if (params.bwa_index != null) {
                ch_bwa_index = Channel.fromPath(params.bwa_index , checkIfExists: true )
        } else {
            ch_bwa_index = BWA_INDEX ( ch_genome_fasta.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions     = ch_versions.mix(BWA_INDEX.out.versions)
        }
        BWA_MEM (
            ch_input_align,
            ch_bwa_index.map { [ [:], it ] }.collect(),
            false
        )
        ch_align_bam         = BWA_MEM.out.bam
        ch_versions        = ch_versions.mix(BWA_MEM.out.versions)
    }

    //TODO: IMARGI
    if (params.align_tool == 'bwa_mem_hisat' && params.exp_type == 'imargi') {

        if (params.bwa_index != null) {
            ch_bwa_index = Channel.fromPath(params.bwa_index , checkIfExists: true )
        } else {
            ch_bwa_index = BWA_INDEX ( ch_genome_fasta.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions     = ch_versions.mix(BWA_INDEX.out.versions)
        }

        BWA_MEM (
            ch_input_align,
            ch_bwa_index.map { [ [:], it ] }.collect(),
            false
        )
        ch_bwa_bam         = BWA_MEM.out.bam
        ch_versions        = ch_versions.mix(BWA_MEM.out.versions)
        ch_bwa_bam.view()
        SAM_TO_FASTQ ( ch_bwa_bam )
        ch_input_align = SAM_TO_FASTQ.out.fastq
        ch_input_align.view()
        if (params.splice_sites == null) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            ch_splicesites = Channel.fromPath(params.splice_sites , checkIfExists: true)
        }
        if (params.hisat2_index != null) {
                ch_hisat2_index = Channel.fromPath(params.hisat2_index , checkIfExists: true )
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_genome_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }

        HISAT2_ALIGN( 
            ch_input_align,
            ch_hisat2_index.map { [ [:], it ] }.collect(),
            ch_splicesites.map { [ [:], it ] }.collect()
        )
        
        ch_align_bam      = HISAT2_ALIGN.out.bam
        ch_align_log       = HISAT2_ALIGN.out.summary
        ch_versions        = ch_versions.mix(HISAT2_ALIGN.out.versions)

        ch_align_bam.view()
    }



    // ch_hisat2_bam.view()
    //TODO: fix
    if (!(params.exp_type in ['rap', 'chirp', 'chart'])) {
        ch_align_bam
        | flatMap { meta, bam -> bam.collect { data -> tuple(meta, data) } } 
        | set { ch_bam }
    }  else {
        ch_bam = ch_align_bam
    }
    // ch_input_bam_filter.view()
    emit:
    bam                = ch_bam
    logs                = ch_align_log                               
    versions           = ch_versions // channel: [ versions.yml ]

}