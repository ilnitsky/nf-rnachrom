
process SMOOTH_INPUT {
    tag "$meta.id"

    conda "${projectDir}/envs/full_env.yml"
    // conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
    // tag "$name"
    publishDir "${params.outdir}/smooth_input", mode: 'copy'

    input:
    tuple val(meta), path(input)
    path(bins)
    path(chromsizes_file)

    output:
    tuple val(meta), path("*_sm.bins.bgr"), emit: smoothed
    tuple val(meta), path("*.log"), emit: log

    script:
    """
    cat <<-END > smoother.cfg
    chrom=${chromsizes_file}
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
