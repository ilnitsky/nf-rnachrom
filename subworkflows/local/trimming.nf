include { BBMAP_BBDUK                 } from '../../modules/nf-core/bbmap/bbduk/main'
include { TRIMMOMATIC                 } from '../../modules/nf-core/trimmomatic/main'
include { TRIMGALORE                  } from '../../modules/nf-core/trimgalore/main' 
include { FASTP                       } from '../../modules/nf-core/fastp/main' 

workflow TRIM {
    take:
    ch_reads // file: /path/to/samplesheet.csv
    ch_adapters_file // file: /path/to/adapters_file.fa

    main:
    ch_stats = Channel.empty()
    ch_versions      = Channel.empty()
    ch_trimmed_reads = Channel.empty()
    ch_adapters      = Channel.fromPath( "$projectDir/bin/adapters/TruSeq3-PE.fa", checkIfExists: true)
    ch_adapters_redc = Channel.fromPath( "$projectDir/assets/adapters_redc.fa", checkIfExists: true)
    ch_adapters      = ch_adapters_redc

    ch_reads_with_adapters = ch_reads.combine(ch_adapters_file)
    reads =  ch_reads_with_adapters.map { meta, reads, adapters -> [meta, reads] }
    ch_adapters_file = ch_reads_with_adapters.map { meta, reads, adapters -> adapters }

    if (params.trim_tool == "trimmomatic") {
        TRIMMOMATIC ( reads )
        ch_trimmed_reads    = TRIMMOMATIC.out.trimmed_reads
        ch_unpaired_reads   = TRIMMOMATIC.out.unpaired_reads
        ch_trim_log  = TRIMMOMATIC.out.summary
        // ch_stats            = ch_stats.mix(TRIMMOMATIC.out.summary)
        ch_versions         = ch_versions.mix(TRIMMOMATIC.out.versions)

    } else if (params.trim_tool == "bbduk") {
        BBMAP_BBDUK ( 
            reads, 
            ch_adapters_file 
        )
        ch_trimmed_reads    = BBMAP_BBDUK.out.reads
        ch_trim_log         = BBMAP_BBDUK.out.log
        // ch_stats            = ch_stats.mix(BBMAP_BBDUK.out.log)
        ch_versions         = ch_versions.mix(BBMAP_BBDUK.out.versions)

    } else if (params.trim_tool == "fastp") {
        FASTP ( 
            reads, 
            ch_adapters_file, 
            true, 
            false, 
            false 
        )
        // FASTP ( reads, true, false, false )
        ch_trimmed_reads    = FASTP.out.reads
        ch_trim_log         = FASTP.out.log
        ch_stats            = FASTP.out.html
        ch_versions         = ch_versions.mix(FASTP.out.versions)

    }  else if (params.trim_tool == "trimgalore") {
        TRIM_GALORE ( reads )
        ch_trimmed_reads    = TRIM_GALORE.out.reads
        ch_unpaired_reads   = TRIM_GALORE.out.unpaired_reads
        ch_trim_log         = TRIM_GALORE.out.log
        ch_stats            = TRIM_GALORE.out.html
        ch_versions         = ch_versions.mix(TRIM_GALORE.out.versions)
    } 

    emit:
    reads    =    ch_trimmed_reads  
    stats    =    ch_stats                                
    versions =    ch_versions // channel: [ versions.yml ]
    logs     =    ch_trim_log

}


