process DETECT_STRAND {
    // tag "$params.trim_tool"
    conda "${projectDir}/envs/full_env.yml"
    // conda "${projectDir}/envs/rnachromprocessing.yaml"
     
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
        
    label 'process_single'
    // errorStrategy 'ignore'
    publishDir (
        path: { "$params.outdir/detect_strand" },
        mode: "copy"
    ) 
        
    input:
    // tuple val(meta), path(contacts)
    tuple val(meta), path(contacts)
    path(gtf_file)
    path(genes_list_file)

    output:
    tuple val(meta), path('res/*.{bed,tab,tab.rc}'), emit: files_fixed_strand
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
    mkdir res
    mv ${contacts[0]} res/${contacts[0]}

    ln -s res/${contacts[0]} ${sample}.tab
    

    cat <<-END_JSON > config.json
    {
      "input_dir":".",
      "output_dir":".",
      "gene_annotation":"${gtf_file}",
      "genes_list":"${genes_list_file}",
      "prefix":"${sample}",
      "exp_groups":{"${params.exp_type}":["${sample}"]}
    }
    END_JSON

    detect-strand -c config.json -v
    strand=\$(grep "${sample}" ${sample}_wins.tsv | awk -F"\\t" '{print \$5}')

    case \$strand in
        "ANTI")
            sed -i 'y/+-/-+/' res/${contacts[0]}
            ;;
        "SAME")
            #Do nothing
            ;;
        *)
            echo "Strand orientation has not been deduced!" 
            exit 0
            ;;
    esac
    
    

    """

}
    // mkdir res
    // mv ${contacts[0]} res/${contacts[0]}
    // echo -e "rna_chr\\trna_start\\trna_end\\trna_strand\\trna_cigar\\tdna_chr\\tdna_start\\tdna_end\\tdna_strand\\tdna_cigar\\tSRR_ID\\tid" > ${prefix}.tab
    // awk  'NR>1{FS="\\t"; OFS=FS}{print \$3, \$4, \$5, \$6, \$7, \$10, \$11, \$12, \$13, \$14, \$1, \$1}' ${contacts[0]} >> ${prefix}.tab


def extractPrefix2(String filename) {
    def matcher = filename =~ /^(.+?)(\.bed|\.tab|\.tab\.rc|\.rc)(\.gz)?$/
    return matcher ? matcher[0][1] : null
}


    // if [[ "${separate_rna_dna}" == "true" ]]; then
    //     mv ${contacts[1]} res/${contacts[1]}
    // fi

    
    // [ ! -f ${prefix}.tabrc ] && ln -s ${contacts[0]} ${prefix}.tabrc
