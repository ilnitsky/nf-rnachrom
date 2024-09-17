process DETECT_STRAND {
    // tag "$params.trim_tool"
    conda "${projectDir}/envs/rnachromprocessing.yaml"
    label 'process_single'
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

    // String filePath = contacts.get(0) // Adjust this line according to your needs   
    String filename = contacts[0]
    String sample = extractPrefix2(filename)

    // filename = new File(contacts[0]).getName()
    // sample = extractPrefix(filename)
    // sample = params.procedure == 'new' ? meta.prefix + '_1' : meta.prefix
    separate_rna_dna = params.procedure == 'new' ? 'true' : ''

    """
    cat <<-END_JSON > config.json
    {
      "input_dir":".",
      "output_dir":".",
      "gtf_annotation":"${params.annot_GTF}",
      "genes_list":"${params.detect_strand_genes_list}",
      "prefix":"${meta.prefix}",
      "exp_groups":{"${params.exp_type}":["${sample}"]}
    }
    END_JSON

    detect-strand -c config.json -v
    strand=\$(grep "${sample}" ${meta.prefix}_wins.tsv | awk -F"\\t" '{print \$5}')

    case \$strand in
        "ANTI")
            sed -i 's/+/temp/g; s/-/+/g; s/temp/-/g' ${contacts[0]}
            ;;
        "SAME")
            #Do nothing
            ;;
        *)
            echo "Strand orientation has not been deduced!" 
            
            ;;
    esac
    mkdir res
    mv ${contacts[0]} res/${contacts[0]}
    if [[ "${separate_rna_dna}" == "true" ]]; then
        mv ${contacts[1]} res/${contacts[1]}
    fi
    """

}


def extractPrefix2(String filename) {
    def matcher = filename =~ /^(.+?)(\.bed|\.tab)(\.gz)?$/
    return matcher ? matcher[0][1] : null
}