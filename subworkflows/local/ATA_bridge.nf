
include { DEDUP          } from '../../modules/local/dedup'
include { TRIMMOMATIC    } from '../../modules/nf-core/trimmomatic/main'
include { PEAR           } from '../../modules/nf-core/pear/main'
include { BITAP_DEBRIDGE } from '../../modules/local/debridge'
include { HISAT2_ALIGN   } from '../../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_VIEW as BAM_FILTER } from '../../modules/nf-core/samtools/view/main'
include { BEDTOOLS_BAMTOBED           } from '../../modules/nf-core/bedtools/bamtobed/main'
include { JULIA_DEBRIDGE_CHARTOOLS; RNA_AND_DNA_PARTS } from '../../modules/local/bridge_redc'
include { CONFIG; RSITES; NUCLEOTIDE_DISTRIBUTION_RSITES; TRIM; ALIGN; BED_FILES;  CONTACTS } from '../../modules/local/rnachromprocessing'
include { XRNA_CONFIG; INFER_XRNA } from '../../modules/local/xrna_assembly'

ch_hisat_index   = params.hisat2_index ? Channel.fromPath(params.hisat2_index) : Channel.empty()
ch_splicesites   = params.splice_sites ? Channel.fromPath(params.splice_sites) : Channel.empty()

ch_config_detect_strand =  Channel.fromPath( "$projectDir/assets/detect_strand.json", checkIfExists: true)
ch_config_xrna          =  Channel.fromPath( "$projectDir/assets/xrna.json", checkIfExists: true)
ch_config               =  Channel.fromPath( "$projectDir/assets/new_config.json", checkIfExists: true)


workflow ATA_BRIDGE {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    reads

    main:
    ch_versions = Channel.empty()

    DEDUP ( reads ) 
    ch_deduplicated = DEDUP.out.deduplicated

    // bioawk -c fastx '{print $seq}' ~/nf-rnachrom/data/imargi/SRR9900120_2.fastq | head -n200000  | cut -c1-2 | sort | uniq -c

    TRIMMOMATIC ( ch_deduplicated )
    ch_trimmed_reads    = TRIMMOMATIC.out.trimmed_reads
    ch_unpaired_reads   = TRIMMOMATIC.out.unpaired_reads
    ch_trimmed_summary  = TRIMMOMATIC.out.summary
    ch_versions         = ch_versions.mix(TRIMMOMATIC.out.versions)
    
    // NUCLEOTIDE_DISTRIBUTION_RSITES( DEDUP4REDC.out.map{ meta, fastq -> fastq }.flatten() )


    //TO DO: If not single end
    PEAR ( ch_trimmed_reads )
    ch_pear_assembled      = PEAR.out.assembled
    ch_pear_unassembled_f  = PEAR.out.unassembled_forward
    ch_pear_unassembled_r  = PEAR.out.unassembled_reverse
    ch_pear_discarded      = PEAR.out.discarded
    ch_pear_stats          = PEAR.out.stats
    ch_versions            = ch_versions.mix(PEAR.out.versions)

    ch_pear = Channel.empty()



    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FIND BRIDGE AND SPLIT RNA AND DNA PARTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    

    if (params.debridge_tool.contains('chartools')){
        JULIA_DEBRIDGE_CHARTOOLS( 
            ch_pear_assembled,
            ch_pear_unassembled_f,
            ch_pear_unassembled_r
        )
        ch_bridge_coords_as      = JULIA_DEBRIDGE_CHARTOOLS.out.coords_as
        ch_bridge_coords_unF     = JULIA_DEBRIDGE_CHARTOOLS.out.coords_unF
        ch_bridge_coords_unR     = JULIA_DEBRIDGE_CHARTOOLS.out.coords_unR
    } else if (params.debridge_tool.contains('bitap')){
        BITAP_DEBRIDGE(
            ch_pear_assembled,
            ch_pear_unassembled_f,
            ch_pear_unassembled_r
        )
        ch_bridge_coords_as      = BITAP_DEBRIDGE.out.coords_as
        ch_bridge_coords_unF     = BITAP_DEBRIDGE.out.coords_unF
        ch_bridge_coords_unR     = BITAP_DEBRIDGE.out.coords_unR
        ch_separated_fastq   = BITAP_DEBRIDGE.out.fastq
        // ch_separated_dna_fastq   = BITAP_DEBRIDGE.out.dna_fastq
        ch_bridge_not_found      = BITAP_DEBRIDGE.out.bridge_not_found_fastq
    }


    // RNA_AND_DNA_PARTS( 
    //     ch_bridge_coords_as,
    //     ch_bridge_coords_unF,
    //     ch_bridge_coords_unR,
    //     ch_pear_assembled,
    //     ch_pear_unassembled_f,
    //     ch_pear_unassembled_r
    //  )
    
    // RNA_AND_DNA_PARTS.out
    // | set { parts_ch }


    HISAT2_ALIGN( 
        // parts_ch,
        ch_separated_fastq,
        ch_hisat_index.map { [ [:], it ] }.collect(),
        ch_splicesites.map { [ [:], it ] }.collect()
    )
    ch_hisat2_bam      = HISAT2_ALIGN.out.bam
    ch_hisat2_summary  = HISAT2_ALIGN.out.summary
    ch_versions        = ch_versions.mix(HISAT2_ALIGN.out.versions)

    ch_hisat2_bam
    | flatMap { meta, bam -> bam.collect { data -> tuple(meta, data) } } 
    | set { ch_input_bam_filter }

    //  [meta, [rna.bam, dna.bam]]   ->    [meta, rna.bam]
    //                                     [meta, dna.bam]


    BAM_FILTER( ch_input_bam_filter )
    ch_filtered_bam     = BAM_FILTER.out.bam
    ch_bam_filter_stat  = BAM_FILTER.out.stat
    ch_versions         = ch_versions.mix(BAM_FILTER.out.versions)

    BEDTOOLS_BAMTOBED( ch_filtered_bam )
    ch_bed_files        = BEDTOOLS_BAMTOBED.out.bed
    ch_versions         = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)



    ch_bed_files                                    //   [meta, rna.bam]    ->    [meta, [rna.bam, dna.bam]] 
    | groupTuple (sort: true)                       //   [meta, dna.bam]   
    | map { meta,bed -> tuple( meta, bed.sort{it.name})}
    | map { meta,bed -> [meta, bed[0], bed[1]] }
    | set {ch_join_bed_raw}
 
    ch_bed_files                                    //   [meta, rna.bam]    ->    [meta, [rna.bam, dna.bam]] 
    | groupTuple (sort: true)                       //   [meta, dna.bam]   
    | map { meta,bed -> tuple( meta, bed.sort{it.name})}
    | map { meta,bed -> [meta, bed[0]] }
    | set { ch_rna_beds }

    ch_bed_files                                    //   [meta, rna.bam]    ->    [meta, [rna.bam, dna.bam]] 
    | groupTuple (sort: true)                       //   [meta, dna.bam]   
    | map { meta,bed -> tuple( meta, bed.sort{it.name})}
    | map { meta,bed -> [meta, bed[1]] }
    | set { ch_dna_beds }

    

    ch_bed_files
    | groupTuple (sort: true)                         
    | map { meta,bed -> tuple( meta, bed.sort{it.name})}
    | set { ch_rna_dna_bed }

    ch_bam_filter_stat
    | groupTuple (sort: true)                         
    | map { meta,file -> tuple( meta, file.sort{it.name})}
    | set {ch_bam_filter_stat_grouped}

    emit:
    trimmed_summary   = TRIMMOMATIC.out.summary
    pear_stats        = PEAR.out.stats
    bridge_coords_as  = ch_bridge_coords_as
    bridge_coords_unF = ch_bridge_coords_unF
    bridge_coords_unR = ch_bridge_coords_unR
    hisat2_summary    = HISAT2_ALIGN.out.summary
    bam_filter_stat   = ch_bam_filter_stat_grouped
    hisat2_bam        = HISAT2_ALIGN.out.bam
    join_bed_raw      = ch_join_bed_raw
    rna_beds          = ch_rna_beds
    dna_beds          = ch_dna_beds
    rna_dna_bed       = ch_rna_dna_bed
    versions          = ch_versions

    // contacts = CONTACTS.out.contacts           // channel: [ config.json ]
    // up = MAPPING.out.stdout_ch
       
}


    // MERGEPEAR4REDC.out.fastq
    // | map { meta, assembled, unassembled_F, unassembled_R -> [meta, assembled] }
    // | set { pear_ch }

    

    // JULIA_DEBRIDGE_CHARTOOLS.out
    // | combine( pear_ch, by: 0 )
    // | set { bridge_coords_ch }


    // ch_pear
    // | mix(ch_pear_assembled)
    // | mix(ch_pear_unassembled_f)
    // | groupTuple
    // | map { meta, left, right -> [meta, left[0], left[1], right[0]] }
    // | set { pear_ch }