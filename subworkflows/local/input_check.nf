//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
// include { SAMPLESHEET_CHECK_RNA } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    binaries

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
    // def matcher = filename =~ /^(.+?)(_[R12]|\.fq[12]|\.fastq|_sequence|_read[12])/
    def matcher = filename =~  /^(.+?)(?:_1|_2|_[rR]1|_[rR]2)?(\.fastq|\.fq|\.fq1|\.fq2)(\.gz)?$/
    return matcher ? matcher[0][1] : null
}

def extractPrefix2(String filename) {
    // Adjusted regex to capture the entire prefix including optional _1 or _2 or _r1 or _R1 before .fastq
    def matcher = filename =~ /^(.+?)((\.fastq|\.fq|\.fq1|\.fq2)(\.gz)?)$/
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

    if (meta.id == "rnaseq") {
        meta.method = "RNA-seq"
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
    rna_filename = new File(row.rna).getName()
    dna_filename = new File(row.dna).getName()

    def meta = [:]
    meta.id         = row.sample
    meta.prefix     = extractPrefix(rna_filename)
    meta.RNA        = extractPrefix2(rna_filename)
    meta.DNA        = extractPrefix2(dna_filename)
    meta.single_end = false
    meta.method     = "ATA"                      // All-to-all type of methods


    if (meta.id == "rnaseq") {
        meta.method = "RNA-seq"
    }
    
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

def create_rnaseq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()
    meta.prefix     = extractPrefix(new File(row.fastq_1).getName())
    meta.method     = "RNA-seq"

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