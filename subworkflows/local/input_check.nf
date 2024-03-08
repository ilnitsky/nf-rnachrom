//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { params.bridge_processing || ['chart', 'rap', 'chirp'].contains(params.exp_type) ? create_fastq_channel(it) : create_fastq_channel_rna_dna(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
    csv = SAMPLESHEET_CHECK.out.csv           // channel: [ samplesheet.csv ]
}

// Match common formats for paired end reads filenames to extract prefix 
def extractPrefix(String filename) {
    def matcher = filename =~ /^(.+?)(_[R12]|\.fq[12]|\.fastq|_sequence|_read[12])/
    return matcher ? matcher[0][1] : null
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    filename = new File(row.fastq_1).getName()

    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()
    meta.prefix     = extractPrefix(filename)

    //Parse input/treatment metadata in OTA libraries
    if (row.containsKey("control")) {
        meta.control = row.control
    }

    if (params.bridge_processing) {
        meta.method = "ATA"
    } else if (['chart', 'rap', 'chirp'].contains(params.exp_type)) {
        meta.method = "OTA"
    }



    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

def create_fastq_channel_rna_dna(LinkedHashMap row) {
    filename = new File(row.rna).getName()
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.prefix     = extractPrefix(filename)
    meta.single_end = false
    meta.method     = "ATA"                      // All-to-all type of methods

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.rna).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> RNA FastQ file does not exist!\n${row.rna}"
    }
   
    if (!file(row.dna).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> DNA FastQ file does not exist!\n${row.dna}"
    }
    fastq_meta = [ meta, [ file(row.rna), file(row.dna) ] ]
    return fastq_meta
}

