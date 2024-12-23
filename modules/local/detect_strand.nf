process DETECT_STRAND {
    // tag "$params.trim_tool"
    conda "${projectDir}/envs/rnachromprocessing.yaml"
    label 'process_single'
    errorStrategy 'ignore'
    publishDir (
        path: { "$params.outdir/detect_strand" },
        mode: "copy"
    ) 
        
    input:
    tuple val(meta), path(contacts)

    output:
    tuple val(meta), path('res/*.{bed,tab}'), emit: files_fixed_strand
    tuple val(meta), path('*_wins.tsv'),  emit: strand_vote_result
    tuple val(meta), path('*.png'),  emit: strand_vote_png


    script:
    // def sample = ''
    def separate_rna_dna = ''
    def prefix = "${meta.prefix}"

    // String filePath = contacts.get(0) // Adjust this line according to your needs   
    String filename = contacts[0]
    String sample = extractPrefix2(filename)

    // filename = new File(contacts[0]).getName()
    // sample = extractPrefix(filename)
    // sample = params.procedure == 'new' ? meta.prefix + '_1' : meta.prefix
    separate_rna_dna = params.procedure == 'new' ? 'true' : ''

    """
    echo -e "rna_chr\\trna_bgn\\trna_end\\tid\\trna_strand\\trna_cigar\\tdna_chr\\tdna_bgn\\tdna_end\\tdna_strand\\tdna_cigar\\tSRR_ID" > blacklist_${prefix}.tab
    awk  'BEGIN{FS="\\t"; OFS=FS}{print \$3, \$4, \$5, \$1, \$6, \$7, \$10, \$11, \$12, \$13, \$14, \$1}' ${contacts[0]} >> blacklist_${prefix}.tab



    cat <<-END_JSON > config.json
    {
      "input_dir":".",
      "output_dir":".",
      "gtf_annotation":"${params.annot_GTF}",
      "genes_list":"${params.detect_strand_genes_list}",
      "prefix":"blacklist_${prefix}",
      "exp_groups":{"${params.exp_type}":["blacklist_${prefix}"]}
    }
    END_JSON

    detect-strand -c config.json -v
    strand=\$(grep "${prefix}" blacklist_${prefix}_wins.tsv | awk -F"\\t" '{print \$5}')

    case \$strand in
        "ANTI")
            sed -i '0,/-/{s/-/±/}; 0,/+/{s/+/-/}; s/±/+/' ${contacts[0]}
            ;;
        "SAME")
            #Do nothing
            ;;
        *)
            echo "Strand orientation has not been deduced!" 
            exit 1
            ;;
    esac
    mkdir res
    mv ${contacts[0]} res/${contacts[0]}

    """

}


def extractPrefix2(String filename) {
    def matcher = filename =~ /^(.+?)(\.bed|\.tab)(\.gz)?$/
    return matcher ? matcher[0][1] : null
}


    // if [[ "${separate_rna_dna}" == "true" ]]; then
    //     mv ${contacts[1]} res/${contacts[1]}
    // fi