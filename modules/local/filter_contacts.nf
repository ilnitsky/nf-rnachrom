process FILTER_CONTACTS {
    conda "${projectDir}/envs/secondary_processing.yml"
    label 'process_single'
    publishDir (
        path: { "$params.outdir/Filter_contacts" },
        mode: "copy"
    ) 
        
    input:
    tuple val(meta), path(unique_tab)   
    // tuple val(meta) 

    output:
    tuple val(meta), path('filtered_*.tab.rc'), emit: filtered_contacts
    tuple val(meta), path('out*.tab.rc'), emit: filtered_out
    tuple val(meta), path('id_reads_*.tab.rc'), emit: ucarna_id
    tuple val(meta), path('*png'), emit: png
    // tuple val(meta), path('*_wins.tsv'),  emit: strand_vote_result


    script:
    // sample = params.procedure == 'new' ? meta.prefix + '_1' : meta.prefix

    def mode = meta.method == "OTA" ? 
        (meta.single_end ? "OTA_SE" : "OTA_PE") : 
        (meta.method == "ATA" ? "ATA, not iMARGI" : null)

    def ucarna_assembly = params.ucarna_assembly ? "no" : "yes" 
    // def max_insert_size = params.pe_insert_size

    """
    python3 ${projectDir}/bin/EditDistance_CIGAR_filter.py \\
        "NM + N_softClipp_bp" 2 2 0 0 300 "${ucarna_assembly}" "not explorer" "${mode}" "${unique_tab}" "./" "./"

    mv filtered_out_"${unique_tab}" out_"${unique_tab}"
    """

}



    // if (meta.method == "OTA" ){
    //     if (meta.single_end) {
    //         """
    //         python3 ${projectDir}/bin/EditDistance_CIGAR_filter.py \\
    //             "NM + N_softClipp_bp" 2 2 0 0 200 "yes" "not explorer" "OTA_SE" "${unique_tab}" "./" "./"
    //         """
    //     } else {
    //         """
    //         python3 ${projectDir}/bin/EditDistance_CIGAR_filter.py \\
    //             "NM + N_softClipp_bp" 2 2 0 0 200 "yes" "not explorer" "OTA_PE" "${unique_tab}" "./" "./"
    //         """
    //     } 
    // } else if (meta.method == "ATA" ){
    //         """
    //         python3 ${projectDir}/bin/EditDistance_CIGAR_filter.py \\
    //             "NM + N_softClipp_bp" 2 2 0 0 200 "yes" "not explorer" "ATA, not iMARGI" "${unique_tab}" "./" "./"
    //         """
    // }

       

    // samtools view -h -F 256 ${rna_prefix}.rna.bam | python3 ${projectDir}/bin/extract_sam_file_matched_seq_to_fastq.py > ${rna_prefix}.rna.fastq
    // samtools view -h -F 256 ${dna_prefix}.dna.bam | python3 ${projectDir}/bin/extract_sam_file_matched_seq_to_fastq.py > ${dna_prefix}.dna.fastq
