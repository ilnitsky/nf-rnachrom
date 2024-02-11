
process XRNA_CONFIG {
    // publishDir ${params.outdir}/'processing', mode: 'copy'


    label 'process_single'

    
    input:
    path samplesheet
    path config_strand
    path config_xrna

    output:
    path '*.strand.json', emit: strand_json
    path '*.xrna.json', emit: xrna_json


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    jq \\
    --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[1]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))" '.exp_groups.${params.exp_type} = \$rna_ids'  $config_strand \\
    | jq '.gtf_annotation="${params.annot_GTF}"' \\
    | jq '.genes_list="${params.detect_strand_genes_list}"' \\
    | jq '.prefix="detect"' \\
    > config.strand.json

    jq \\
    --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[1]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
     '. + {rna_ids: \$rna_ids}' $config_xrna > config.xrna.json

    """

}

process DETECT_STRAND {
    // tag "$params.trim_tool"
    label 'process_single'
    publishDir (
        path: { "$params.outdir/detect_strand" },
        mode: "copy"
    ) 
        
    input:
    path(contacts)
    path(config)

    output:
    path '*'
    path 'detect_wins.tsv', emit: strand_vote_result


    script:
    """
    detect-strand -c $config -v
    """

}

// process INFER_XRNA {

// }