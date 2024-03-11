process RSITES {
    // tag "$params.trim_tool"
    conda "${projectDir}/envs/rnachromprocessing.yaml"
    label 'process_single'
    publishDir (
        path: { "$params.outdir/rsites" },
        mode: "copy",
        pattern: "rsites/*.fastq",
        saveAs: { fn -> file(fn).name }
    ) 
        
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('rsites/*.fastq'), emit: fastq
    // tuple val(meta), path('*_wins.tsv'),  emit: strand_vote_result


    script:
    def sample = ''
    def sample_rna = meta.RNA
    def sample_dna = meta.DNA
    // sample = params.procedure == 'new' ? meta.prefix + '_1' : meta.prefix

    """
    cat <<-END_JSON > config.json
    {
        "rna_ids": ["${sample_rna}"],
        "dna_ids": ["${sample_dna}"],
        "base_dir" : ".",
        "input_dir": ".",
        "output_dir": ".",
        "cpus": ${task.cpus},
        "keep" : ["rsites"],
        "rsites": {
            "type": "${params.exp_type}"
        }
    }
    END_JSON

    rnachromprocessing -c config.json -s rsites -v
    

    """

}