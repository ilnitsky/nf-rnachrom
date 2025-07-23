include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/gunzip/main'
include { HISAT2_ALIGN                  } from '../../modules/nf-core/hisat2/align/main'  
include { STAR_ALIGN                    } from '../../modules/nf-core/star/align/main'
include { BOWTIE2_ALIGN                 } from '../../modules/nf-core/bowtie2/align/main' 
include { SAM_TO_FASTQ                  } from '../../modules/local/sam_to_fastq'
include { BWA_MEM                       } from '../../modules/nf-core/bwa/mem/main' 

// include { STAR_CHIMERIC_READS           } from '../../modules/local/star_chimeric'    

workflow ALIGN {
    take:
    ch_input_align // file: /path/to/samplesheet.csv
    ch_hisat2_index // HISAT2 index (required if using HISAT2)
    ch_star_index   // STAR index (required if using STAR)
    ch_bowtie2_index // Bowtie2 index (required if using Bowtie2)
    ch_bwa_index    // BWA index (required if using BWA)
    ch_splicesites  // HISAT2 splice sites (required if using HISAT2)
    ch_genome_fasta // Genome FASTA file
    ch_gtf          // GTF annotation
    ch_aligner      // Aligner tool

    main:
    ch_versions     = Channel.empty()
    ch_align_log    = Channel.empty()

    if (ch_aligner == 'hisat2') {
        HISAT2_ALIGN ( 
            ch_input_align,
            ch_hisat2_index.map { [ [:], it ] }.collect(),
            ch_splicesites.map { [ [:], it ] }.collect()
        )
        ch_align_bam       = HISAT2_ALIGN.out.bam
        ch_align_log       = HISAT2_ALIGN.out.summary
        ch_versions        = ch_versions.mix(HISAT2_ALIGN.out.versions)
    }

    if (ch_aligner == 'star') {
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

    }

    if (ch_aligner == 'bowtie2') {
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

    if (ch_aligner == 'bwa_mem') {
        BWA_MEM (
            ch_input_align,
            ch_bwa_index.map { [ [:], it ] }.collect(),
            false
        )
        ch_align_bam         = BWA_MEM.out.bam
        ch_versions        = ch_versions.mix(BWA_MEM.out.versions)
    }

    //TODO: IMARGI
    if (ch_aligner == 'bwa_mem_hisat' && params.exp_type == 'imargi') {
        BWA_MEM (
            ch_input_align,
            ch_bwa_index.map { [ [:], it ] }.collect(),
            false
        )
        ch_bwa_bam         = BWA_MEM.out.bam
        ch_versions        = ch_versions.mix(BWA_MEM.out.versions)
        
        SAM_TO_FASTQ ( ch_bwa_bam )
        ch_input_align = SAM_TO_FASTQ.out.fastq

        HISAT2_ALIGN( 
            ch_input_align,
            ch_hisat2_index.map { [ [:], it ] }.collect(),
            ch_splicesites.map { [ [:], it ] }.collect()
        )
        
        ch_align_bam      = HISAT2_ALIGN.out.bam
        ch_align_log       = HISAT2_ALIGN.out.summary
        ch_versions        = ch_versions.mix(HISAT2_ALIGN.out.versions)
    }

    emit:
    bam                = ch_align_bam
    logs                = ch_align_log                               
    versions           = ch_versions // channel: [ versions.yml ]
}