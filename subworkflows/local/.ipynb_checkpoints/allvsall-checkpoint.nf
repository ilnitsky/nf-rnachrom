include { RSITES } from '../../modules/local/rsites'
include { MAPPING } from '../../modules/local/mapping'
include { CIGAR_FILTER; BLACKLISTING; MERGE_REPLICAS } from '../../modules/local/filter'
include { ANNOTATION; NORMALIZE_N2 } from '../../modules/local/annotation'

workflow ALLVSALL {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    config      // file: /path/to/config.json
    reads

    main:
    RSITES( samplesheet, config, reads )
    MAPPING( RSITES.out.json, RSITES.out.rsites )
    CIGAR_FILTER( MAPPING.out.contacts.flatten() )
    BLACKLISTING( CIGAR_FILTER.out )
    MERGE_REPLICAS( samplesheet, BLACKLISTING.out )

    
    // ANNOTATION.out
    // | flatten()
    // | filter{ file -> file.name =~ /.*4-voted.tab/ }
    // | set { ch_voted }

    // NORMALIZE_N2( ch_voted )
    

    emit:
    contacts = MAPPING.out.contacts           // channel: [ config.json ]
    // up = MAPPING.out.stdout_ch
       
}