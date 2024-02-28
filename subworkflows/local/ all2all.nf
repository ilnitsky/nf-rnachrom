include { GET_RAW_CONTACTS_A2A               } from '../../modules/local/get_raw_contacts_a2a'


workflow{
    take:
    samplesheet
    config

    main:
    ch_versions = Channel.empty()
    // Channel.fromPath( params.reads, checkIfExists:true )
    GET_RAW_CONTACTS_A2A(
        
    )



}