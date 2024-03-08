
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
    if (!params.bridge_processing){  
    """
    jq \\
    --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[1]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))" '.exp_groups.${params.exp_type} = \$rna_ids'  $config_strand \\
    | jq '.gtf_annotation="${params.annot_GTF}"' \\
    | jq '.genes_list="${params.detect_strand_genes_list}"' \\
    | jq '.prefix="detect"' \\
    > config.strand.json

    jq \\
    --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[1]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
     '. + {rna_ids: \$rna_ids}' $config_xrna \\
    | jq '.annotation.gtf_annotation="${params.annot_GTF}"' \\
    | jq '.annotation.strand_info="detect_wins.tsv"' \\
    | jq '.hisat.genome_path="${params.genome_path}"' \\
    | jq '.hisat.known_splice="${params.splice_sites}"' \\
    | jq '.prefix="detect"' \\
      > config.xrna.json
    """
    } else{
    """
    jq \\
    --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[2]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))" '.exp_groups.${params.exp_type} = \$rna_ids'  $config_strand \\
    | jq '.gtf_annotation="${params.annot_GTF}"' \\
    | jq '.genes_list="${params.detect_strand_genes_list}"' \\
    | jq '.prefix="detect"' \\
    > config.strand.json

    jq \\
    --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[2]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
     '. + {rna_ids: \$rna_ids}' $config_xrna \\
    | jq '.annotation.gtf_annotation="${params.annot_GTF}"' \\
    | jq '.annotation.strand_info="detect_wins.tsv"' \\
    | jq '.hisat.genome_path="${params.genome_path}"' \\
    | jq '.hisat.known_splice="${params.splice_sites}"' \\
    | jq '.prefix="detect"' \\
     > config.xrna.json
    """   
    }


}

process INFER_XRNA {
    // tag "$params.trim_tool"

    publishDir (
        path: { "$params.outdir/XRNA_assembly" },
        mode: "copy"
    ) 
        
    input:
    path(bed)
    path(fastq)
    path(config)
    path(strand_vote_result)

    output:
    path '*'

    script:
    """
    infer-xrna -c $config -vv
    """

}