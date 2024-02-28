include { CONFIG; DEDUP; RSITES; NUCLEOTIDE_DISTRIBUTION_RSITES; TRIM; ALIGN; BAM_FILTER; BED_FILES;  CONTACTS } from '../../modules/local/rnachromprocessing'
include { SPLICING_RNA; MERGE_REPLICAS; SPLIT_BY_CHRS;  ANNOTATE_RNA } from '../../modules/local/rnachromprocessing'
include { BACKGROUND; NORMALIZE_RAW; NORMALIZE_N2; SCALING } from '../../modules/local/rnachromprocessing'
include { VALIDATE_ANNOT; BARDIC } from '../../modules/local/rnachromprocessing'
include { XRNA_CONFIG; DETECT_STRAND } from '../../modules/local/xrna_assembly'


ch_config_detect_strand =  Channel.fromPath( "$projectDir/assets/detect_strand.json", checkIfExists: true)
ch_config_xrna          =  Channel.fromPath( "$projectDir/assets/xrna.json", checkIfExists: true)
ch_config               =  Channel.fromPath( "$projectDir/assets/new_config.json", checkIfExists: true)


workflow ATA {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    reads

    main:

    // rna_samples = Channel.fromPath( params.input )
    // .splitCsv(header: true, sep:',')
    // .map{row -> tuple(file(row.rna).baseName, file(row.rna))}

    CONFIG( samplesheet, ch_config )
    XRNA_CONFIG( samplesheet, ch_config_detect_strand, ch_config_xrna )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PREPROCESSING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //TO DO: margi ; DEDUP -> TRIM ?

    if (params.exp_type == "imargi"){
        RSITES( reads, CONFIG.out.json )
        DEDUP( RSITES.out.rsites, CONFIG.out.json )
        TRIM( DEDUP.out.dedup, CONFIG.out.json )
    } else {
        DEDUP( reads, CONFIG.out.json )
        RSITES( DEDUP.out.dedup, CONFIG.out.json )
        TRIM( RSITES.out.rsites, CONFIG.out.json )
    }
    
    NUCLEOTIDE_DISTRIBUTION_RSITES( RSITES.out.rsites.flatten() )

    ALIGN( TRIM.out.trim, CONFIG.out.json )
    BAM_FILTER( ALIGN.out.align, CONFIG.out.json )
    BED_FILES( BAM_FILTER.out.bam, CONFIG.out.json )


    if (params.procedure == 'new'){
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // NEW RNACHROMPROCESSING PROCEDURE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
        BED_FILES.out.bed
        | flatten
        | map { file -> tuple(file.baseName, file) }
        | set { bed_ch }

        // rna_beds_ch = bed_ch.join(rna_samples).map{name, bed, fastq -> tuple( name, bed)}.view()

        SPLICING_RNA(rna_beds_ch)

        //ungroup sample and chromosome splits
        MERGE_REPLICAS( samplesheet, SPLICING_RNA.out.map{name, cigar -> cigar}.collect() )

        SPLIT_BY_CHRS( MERGE_REPLICAS.out.map{merged -> tuple(file(merged).baseName, file(merged))} )

        ANNOTATE_RNA( SPLIT_BY_CHRS.out.map{name, split -> split}.flatten() )
        // ANNOTATE_RNA.out.view()


    } else {
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Old rnachromprocessing procedure
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
        
        CONTACTS( BED_FILES.out.bed, CONFIG.out.json )
        CONTACTS.out.contacts
        | set { final_table_ch }

        SPLICING_RNA( CONTACTS.out.contacts.flatten().map{ file -> tuple(file.baseName, file) } )

        DETECT_STRAND( final_table_ch, XRNA_CONFIG.out.strand_json )
        MERGE_REPLICAS( samplesheet, DETECT_STRAND.out.strand_vote_result, SPLICING_RNA.out.map{name, cigar -> cigar}.collect() )
        SPLIT_BY_CHRS(  MERGE_REPLICAS.out.flatten().map{merged -> tuple(file(merged).baseName, file(merged))}  )

        ANNOTATE_RNA(  SPLIT_BY_CHRS.out.transpose()  )
        ANNOTATE_RNA.out
        | transpose
        | map { group, annotated -> [annotated.name.split('_')[0], group, annotated.name.split('merged*...')[1], annotated]}
        | branch { 
            voted:      it[2] == 'voted.tab'
            singletons: it[2] == 'singletons.tab'
            selected:   it[2] == 'selected_annot.tab'
            complement: it[2] == 'complement_annot.tab'
            }
        | set { cfiles_ch }

        cfiles_ch.voted
        | map { chr, group, extension, file -> [group, file] }
        | collectFile(storeDir: "$params.outdir/voted", keepHeader: true, sort: true) { group, file -> [ "${group}.voted.tab", file.text] }
        | map { voted -> tuple(file(voted).name.split('.voted.tab')[0], file(voted))}
        | set { voted_ch }


        BACKGROUND( voted_ch )
        | map { bgr -> tuple(file(bgr).name.split('.5-background_sm.bgr')[0], file(bgr))}
        | set { bgr_ch }

        cfiles_ch.voted
        | map{ chr,  group, extension, file -> [group, file] }
        | combine( bgr_ch, by: 0 )
        | set { norm_raw_ch }

        NORMALIZE_RAW( norm_raw_ch )

        NORMALIZE_RAW.out.raw_stat
        | collectFile(storeDir: "$params.outdir/Normalize_raw", keepHeader: true) { group, file -> [ "${group}.5-N2_raw_merged.stat.tab", file.text] }
        | map{stat -> tuple(file(stat).name.split('.5-N2_raw_merged.stat')[0], file(stat))}
        | combine( NORMALIZE_RAW.out.raw_norm, by: 0 )
        | NORMALIZE_N2 
        | SCALING
        // | view

        Channel.fromPath(params.annot_BED).ifEmpty { exit 1, "Input file not found: ${params.annot_BED}" }
        | set { annot_file_ch }

        VALIDATE_ANNOT( annot_file_ch )

        bed6_annot_files_ch = VALIDATE_ANNOT.out.flatten().filter { it.toString().endsWith('.0-corrected_annot.bed6') }
        bed_annot_files_ch  = VALIDATE_ANNOT.out.flatten().filter { it.toString().endsWith('.0-corrected_annot.bed') }

        BARDIC( voted_ch, bed6_annot_files_ch )
        
 
    }

    





    // CIGAR_FILTER( CONTACTS.out.contacts.flatten() )
    // BLACKLISTING( CIGAR_FILTER.out )
    // MERGE_REPLICAS( samplesheet, BLACKLISTING.out )

    
    // ANNOTATION.out
    // | flatten()
    // | filter{ file -> file.name =~ /.*4-voted.tab/ }
    // | set { ch_voted }
    

    emit:
    contacts = CONTACTS.out.contacts           // channel: [ config.json ]
    // up = MAPPING.out.stdout_ch
       
}


    //    ANNOTATE_RNA.out
    //     | branch { 
    //         complement: it.name.contains('complement')
    //         singletons: it.name.contains('singletons')
    //         selected: it.name.contains('selected')
    //         voted: it.name.contains('voted')
    //      }
    //     | set { collected_annot_ch }
    //     collected_annot_ch.complement.collectFile(name: 'TC.txt', storeDir: '/Users/mribeirodantas/sandbox')