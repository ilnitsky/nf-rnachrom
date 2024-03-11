include { BBMAP_BBDUK                 } from '../../modules/nf-core/bbmap/bbduk/main'
include { TRIMMOMATIC                 } from '../../modules/nf-core/trimmomatic/main'

workflow TRIM {
    take:
    reads // file: /path/to/samplesheet.csv

    main:
    ch_versions      = Channel.empty()
    ch_trimmed_reads = Channel.empty()
    ch_adapters      = Channel.fromPath( "$projectDir/bin/adapters/TruSeq3-PE.fa", checkIfExists: true)

    if (params.trim_tool == "trimmomatic") {
        TRIMMOMATIC ( reads )
        ch_trimmed_reads    = TRIMMOMATIC.out.trimmed_reads
        ch_unpaired_reads   = TRIMMOMATIC.out.unpaired_reads
        ch_trimmed_summary  = TRIMMOMATIC.out.summary
        // ch_stats            = ch_stats.mix(TRIMMOMATIC.out.summary)
        ch_versions         = ch_versions.mix(TRIMMOMATIC.out.versions)
    } else if (params.trim_tool == "bbduk") {
        BBMAP_BBDUK ( reads, ch_adapters )
        ch_trimmed_reads    = BBMAP_BBDUK.out.reads
        ch_trim_log         = BBMAP_BBDUK.out.log
        // ch_stats            = ch_stats.mix(BBMAP_BBDUK.out.log)
        ch_versions         = ch_versions.mix(BBMAP_BBDUK.out.versions)
    }

    emit:
    reads    =    ch_trimmed_reads                                  
    versions =    ch_versions // channel: [ versions.yml ]

}


