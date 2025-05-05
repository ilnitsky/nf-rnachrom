
process SMOOTH_INPUT {
    tag "$meta.id"

    conda "${projectDir}/envs/full_env.yml"
    // conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"
    // tag "$name"
    publishDir "${params.outdir}/smooth_input", mode: 'copy'

    input:
    tuple val(meta), path(input)
    path(bins)

    output:
    tuple val(meta), path("*_sm.bins.bgr"), emit: smoothed
    tuple val(meta), path("*.log"), emit: log

    script:
    """
    cat <<-END > smoother.cfg  
    chrom=${params.chromsizes}             
    profPath =./profiles              
    trackPath=.
    resPath=./res
    log=${meta.id}.log
    #============ Prepare parameters
    bin=${params.binsize}                                                
    #============ Statistics parameters
    wSize      =1000000                             
    flankSize  =10000                                       
    kernelSigma=3000.                               
    kernelType =NORMAL                              
    BufSize=40000000
    verbose=1   
    END


    LC_NUMERIC="C" awk '{print \$1, int((\$2+\$3)/2), int((\$2+\$3)/2+1)}' ${input} | \\
	    sort-bed --tmpdir . --max-mem ${task.memory.mega}M - > ${meta.id}.tmp

    bedmap  --echo --count --delim '\\t' ${bins} ${meta.id}.tmp  > ${meta.id}.bgr

    ${projectDir}/bin/Smoother cfg=smoother.cfg  ${meta.id}.bgr
    bedmap --echo --echo-map-id --delim '\\t' ${bins} <(sort-bed ${meta.id}_sm.bgr) > ${meta.id}_sm.bins.bgr
    """
}
