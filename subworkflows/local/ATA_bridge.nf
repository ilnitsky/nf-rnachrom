
include { DEDUP          } from '../../modules/local/dedup'
include { TRIMMOMATIC    } from '../../modules/nf-core/trimmomatic/main'
include { PEAR           } from '../../modules/nf-core/pear/main'
include { HISAT2_ALIGN   } from '../../modules/nf-core/hisat2/align/main'
include { JULIA_DEBRIDGE_CHARTOOLS; BITAP_DEBRIDGE; RNA_AND_DNA_PARTS } from '../../modules/local/bridge_redc'
include { CONFIG; RSITES; NUCLEOTIDE_DISTRIBUTION_RSITES; TRIM; ALIGN; BAM_FILTER; BED_FILES;  CONTACTS } from '../../modules/local/rnachromprocessing'
include { SPLICING_RNA; MERGE_REPLICAS; SPLIT_BY_CHRS;  ANNOTATE_RNA } from '../../modules/local/rnachromprocessing'
include { BACKGROUND; NORMALIZE_RAW; NORMALIZE_N2; SCALING } from '../../modules/local/rnachromprocessing'
include { VALIDATE_ANNOT; BARDIC } from '../../modules/local/rnachromprocessing'
include { XRNA_CONFIG; DETECT_STRAND; INFER_XRNA } from '../../modules/local/xrna_assembly'

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

    // ch_pear
    // | mix(ch_pear_assembled)
    // | mix(ch_pear_unassembled_f)
    // | groupTuple
    // | map { meta, left, right -> [meta, left[0], left[1], right[0]] }
    // | set { pear_ch }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FIND BRIDGE AND SPLIT RNA AND DNA PARTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    JULIA_DEBRIDGE_CHARTOOLS( 
        ch_pear_assembled,
        ch_pear_unassembled_f,
        ch_pear_unassembled_r
     )
    ch_bridge_coords_as      = JULIA_DEBRIDGE_CHARTOOLS.out.coords_as
    ch_bridge_coords_unF     = JULIA_DEBRIDGE_CHARTOOLS.out.coords_unF
    ch_bridge_coords_unR     = JULIA_DEBRIDGE_CHARTOOLS.out.coords_unR

    // JULIA_DEBRIDGE_CHARTOOLS.out
    // | combine( pear_ch, by: 0 )
    // | set { bridge_coords_ch }

    BITAP_DEBRIDGE(
        ch_pear_assembled,
        ch_pear_unassembled_f,
        ch_pear_unassembled_r
    )
    // BITAP_DEBRIDGE.out.coords_as.view()

    RNA_AND_DNA_PARTS( 
        ch_bridge_coords_as,
        ch_bridge_coords_unF,
        ch_bridge_coords_unR,
        ch_pear_assembled,
        ch_pear_unassembled_f,
        ch_pear_unassembled_r
     )
    
    CONFIG( samplesheet, ch_config )
    XRNA_CONFIG( samplesheet, ch_config_detect_strand, ch_config_xrna )


    RNA_AND_DNA_PARTS.out
    // | flatten
    // | collect
    | set { parts_ch }


    parts_ch.view()


    HISAT2_ALIGN( 
        parts_ch,
        ch_hisat_index.map { [ [:], it ] }.collect(),
        ch_splicesites.map { [ [:], it ] }.collect()
    )
    ch_hisat2_bam      = HISAT2_ALIGN.out.bam
    ch_hisat2_summary  = HISAT2_ALIGN.out.summary
    ch_versions        = ch_versions.mix(HISAT2_ALIGN.out.versions)

    HISAT2_ALIGN.out.summary.view()

    emit:
    trimmed_summary   = TRIMMOMATIC.out.summary
    pear_stats        = PEAR.out.stats
    bridge_coords_as  = JULIA_DEBRIDGE_CHARTOOLS.out.coords_as
    bridge_coords_unF = JULIA_DEBRIDGE_CHARTOOLS.out.coords_unF
    bridge_coords_unR = JULIA_DEBRIDGE_CHARTOOLS.out.coords_unR
    hisat2_summary    = HISAT2_ALIGN.out.summary
    hisat2_bam        = HISAT2_ALIGN.out.bam
    versions          = ch_versions

    // contacts = CONTACTS.out.contacts           // channel: [ config.json ]
    // up = MAPPING.out.stdout_ch
       
}


    // MERGEPEAR4REDC.out.fastq
    // | map { meta, assembled, unassembled_F, unassembled_R -> [meta, assembled] }
    // | set { pear_ch }