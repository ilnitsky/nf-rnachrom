include { FASTUNIQ                    } from '../../modules/local/dedup/fastuniq'
include { FASTQ_DUPAWAY               } from '../../modules/local/dedup/fastq_dupaway'
include { SEQKIT_RMDUP                } from '../../modules/local/dedup/seqkit_rmdup'
include { BBMAP_CLUMPIFY              } from '../../modules/nf-core/bbmap/clumpify/main'

workflow DEDUP {
    take:
    reads // file: /path/to/samplesheet.csv

    main:
    ch_versions     = Channel.empty()
    ch_deduplicated = Channel.empty()
    if (params.dedup_tool == "fastq-dupaway") {
        FASTQ_DUPAWAY ( reads )
        ch_deduplicated = FASTQ_DUPAWAY.out.reads
        ch_versions     = ch_versions.mix(FASTQ_DUPAWAY.out.versions)
    } else if (params.dedup_tool == "fastuniq") {  
        FASTUNIQ ( reads )
        ch_deduplicated = FASTUNIQ.out.reads
        ch_versions     = ch_versions.mix(FASTUNIQ.out.versions)
//    } else {params.dedup_tool == "fastuniq" && meta.single_end} {
//        error "Error: Fastuniq only processes paired-end files, please use 'fastq-dupaway' or other tool."
    } else if (params.dedup_tool == "clumpify") {  
        BBMAP_CLUMPIFY ( reads )
        ch_deduplicated     = BBMAP_CLUMPIFY.out.reads
        ch_deduplicated_log = BBMAP_CLUMPIFY.out.log
        ch_versions         = ch_versions.mix(BBMAP_CLUMPIFY.out.versions)
    } else if (params.dedup_tool == "seqkit_rmdup") {
        SEQKIT_RMDUP ( reads )
        ch_deduplicated = SEQKIT_RMDUP.out.reads
        ch_versions     = ch_versions.mix(SEQKIT_RMDUP.out.versions)
    }

    // ch_deduplicated.view{"ch_deduplicated: $it"}

    emit:
    reads    =    ch_deduplicated                                  
    versions =    ch_versions // channel: [ versions.yml ]
    // logs     =    ch_deduplicated_log

}