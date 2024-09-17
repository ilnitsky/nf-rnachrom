
// TODO: params.splicesites -> splicesites variable

process INFER_XRNA {
    // tag "$params.trim_tool"
    conda "${projectDir}/envs/rnachromprocessing.yaml"
    publishDir (
        path: { "$params.outdir/XRNA_assembly" },
        mode: "copy"
    ) 
        
    input:
    path(bed)
    path(fastq)
    path(strand_vote_result)

    output:
    path '*'

    script:
    def rna_ids = ''

    """
    cat <<-END_JSON > xrna_config.json
    {
        "bed_input_dir": ".",
        "fq_input_dir": ".",
        "output_dir": "xrnas",
        "ouputs_prefix": "xrna",
        "keep_extras": [
            "ids",
            "merged_bam"
        ],
        "cpus": ${task.cpus},
        "rna_ids": [
            ${rna_ids}
        ],
        "annotation": {
            "gtf_annotation": "${params.annot_GTF}",
            "strand_info": "${strand_vote_result}"
        },
        "hisat": {
            "genome_path": "${params.genome_path}",
            "known_splice": "${params.splice_sites}",
            "hisat_threads": ${ 2 * task.cpus },
            "cpus": ${task.cpus}
        },
        "stringtie": {
            "stringtie_threads": ${ 2 * task.cpus },
            "cpus": ${task.cpus}
        }
    }
    END_JSON


    infer-xrna -c xrna_config.json -vv
    """

}
        // "rna_ids": [
        //     "SRR5035944",
        //     "SRR5035946",
        //     "SRR9201799",
        //     "SRR9201801"
        // ],

// process XRNA_CONFIG {
//     // publishDir ${params.outdir}/'processing', mode: 'copy'

//     conda "${projectDir}/envs/rnachromprocessing.yaml"
//     label 'process_single'

    
//     input:
//     path samplesheet
//     path config_strand
//     path config_xrna

//     output:
//     path '*.strand.json', emit: strand_json
//     path '*.xrna.json', emit: xrna_json


//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     if (!params.bridge_processing){  
//     """

//     jq \\
//     --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[1]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
//      '. + {rna_ids: \$rna_ids}' $config_xrna \\
//     | jq '.annotation.gtf_annotation="${params.annot_GTF}"' \\
//     | jq '.annotation.strand_info="detect_wins.tsv"' \\
//     | jq '.hisat.genome_path="${params.genome_path}"' \\
//     | jq '.hisat.known_splice="${params.splice_sites}"' \\
//     | jq '.prefix="detect"' \\
//       > config.xrna.json
//     """
//     } else{
//     """

//     jq \\
//     --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[2]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
//      '. + {rna_ids: \$rna_ids}' $config_xrna \\
//     | jq '.annotation.gtf_annotation="${params.annot_GTF}"' \\
//     | jq '.annotation.strand_info="detect_wins.tsv"' \\
//     | jq '.hisat.genome_path="${params.genome_path}"' \\
//     | jq '.hisat.known_splice="${params.splice_sites}"' \\
//     | jq '.prefix="detect"' \\
//      > config.xrna.json
//     """   

//     }


// }