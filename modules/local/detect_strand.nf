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


    script:
    def sample = ''
    def separate_rna_dna = ''
    sample = params.procedure == 'new' ? meta.prefix + '_1' : meta.prefix
    separate_rna_dna = params.procedure == 'new' ? 'true' : ''

    """
    cat <<-END_JSON > ${meta.prefix}.config.json
    {
      "input_dir":".",
      "output_dir":".",
      "gtf_annotation":"${params.annot_GTF}",
      "genes_list":"${params.detect_strand_genes_list}",
      "prefix":"${meta.prefix}",
      "exp_groups":{"${params.exp_type}":["${sample}"]}
    }
    END_JSON
    
    detect-strand -c ${meta.prefix}.config.json -v
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
            exit 1
            ;;
    esac
    mkdir res
    mv ${contacts[0]} res/${contacts[0]}
    if [[ "${separate_rna_dna}" == "true" ]]; then
        mv ${contacts[1]} res/${contacts[1]}
    fi
    """

}