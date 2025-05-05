ANSI_RESET = "\u001B[0m";
ANSI_PURPLE = "\u001B[35m";

def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }


process FASTQ_DUPAWAY {
    //TODO: CHECK FOR CORRECTNESS OF DEDUP TOOL NAME
    // TODO: Add default compare seq mode for fastq-dupaway
    //TODO: conda create -n gcc-boost -c conda-forge gcc make boost    export BOOST_ROOT=\$CONDA_PREFIX/include
    tag "$meta.id,$meta.prefix"
    publishDir (
        path: { "$params.outdir/dedup" },
        mode: "copy",
        pattern: "dedup/*",
        saveAs: { fn -> file(fn).name } 
    )
    conda "${projectDir}/envs/full_env.yml"
    
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'ilnitsky/nf-rnachrom:latest' : '' }"
    // conda 'conda-forge::boost=1.84.0'
    

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("dedup/*.fastq"), emit: reads
    path "versions.yml"                  , emit: versions


    script:
    println print_purple("Started deduplication " + meta.prefix  + " with fastq-dupaway" )
    def args                 = task.ext.args ?: ''
    def prefix               = task.ext.prefix ?: "${meta.prefix}"
    def fastq_dupaway_input  = ''
    def fastq_dupaway_output = ''
    
    fastq_dupaway_input = meta.single_end ? "-i ${reads[0]}" : "-i ${reads[0]} -u ${reads[1]}"

    fastq_dupaway_output = ( params.bridge_processing == false && meta.method == "ATA" ) ? "-o dedup/1_${meta.RNA}.rna.dedup.fastq -p dedup/2_${meta.DNA}.dna.dedup.fastq" :
        (meta.single_end ? "-o dedup/${prefix}.dedup.fastq" : "-o dedup/${prefix}_1.dedup.fastq -p dedup/${prefix}_2.dedup.fastq")

    
    """
    export BOOST_ROOT=\$CONDA_PREFIX/include
    [[ ! -d "dedup" ]] && mkdir dedup
    
    ${projectDir}/bin/fastq-dupaway \
            ${fastq_dupaway_input} \
            ${fastq_dupaway_output} \
            -m ${task.memory.mega} \
            $args

            

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       fastq-dupaway: 1.0
    END_VERSIONS
    """
    
}


        // """
        // ${projectDir}/bin/fastq-dupaway \
        //     -i ${reads[0]} -u ${reads[1]} \
        //     -o ${meta.prefix}_1.uniq.fastq -p ${meta.prefix}_2.uniq.fastq \
        //     -m ${task.memory.mega} \
        //     --format fastq \
        //     --compare-seq ${params.fastq_dupaway_seq_compare}
        // """

// process DEDUP4REDC {
//   // log.info nfcoreHeader()
//   // tag "$meta"
    
//   publishDir (
//       path: { "$params.outdir/dedup" },
//       mode: "copy"
//   )

//   input:
//   tuple val(meta), path(reads)

//   output:
//   tuple val(meta), path("*.uniq.fastq"), emit: fastq
//   path("dedup.stats"), emit: stats

//   script:
//   println print_purple("Started deduplication " + meta.prefix  )
//   """
//   printf "${reads[0]}\\n${reads[1]}" > input.list
//   fastuniq -i input.list -t q -c 0 -o ${meta.prefix}_1.uniq.fastq -p ${meta.prefix}_2.uniq.fastq

//   """

// }


// else if (params.dedup_tool == "fastuniq" && !meta.single_end) {
//         """
//         printf "${reads[0]}\\n${reads[1]}" > input.list
//         fastuniq -i input.list -t q -c 0 -o ${meta.prefix}_1.uniq.fastq -p ${meta.prefix}_2.uniq.fastq
//         """
//     } 