process BACKGROUND {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$voted.baseName"

    publishDir (
        path: { "$params.outdir/background" },
        mode: "copy"
    ) 
    // memory '15 GB'

    input:
    tuple val(name), path(voted)

    output:
    path("*.5-background_sm.bgr")

    script:

    """
    cat <<-END > smoother.cfg
    chrom    =CHRSIZE               # Chromosomes names/sizes
    profPath =REPLACE/profiles              # path to binary profiles (output for prepare)
    trackPath=REPLACE
    resPath=REPLACE/res
    log=log
    #============ Prepare parameters
    bin=500                                                 # step for bynary profile
    #============ Statistics parameters
    wSize      =1000000                             # size of widow (nucleotides)
    flankSize  =10000                                       # size of flanks(nucleotides)
    kernelSigma=3000.                               # kernel width (nucleotides)
    kernelType =NORMAL                              # type of the kernel: NORM | LEFT_EXP | RIGHT_EXP
    BufSize=40000000
    verbose=1   
    END

    python3 ${projectDir}/bin/rnachrom_pipeline_faster/background.py \\
    ${voted} ${projectDir}/bin/Smoother smoother.cfg ${params.chromsizes} \\
    --filter_top_n 50 --filter_tail_n 1000 --outdir .
    """

}

// Rscript =1 