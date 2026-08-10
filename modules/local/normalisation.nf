process NORMALISATION {
    conda "${projectDir}/envs/full_env.yml"
    // conda "${projectDir}/envs/secondary_processing.yml"
    
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
        
     tag "$voted.baseName"

    publishDir (
        path: { "$params.outdir/background" },
        mode: "copy"
    ) 
    // memory '15 GB'

    input:
    tuple val(meta), path(voted)
    path(chromsizes)

    output:
    tuple val(meta), path("*.intersected.N2.bed"),     emit: normalized

    script:

    def prefix     = task.ext.prefix ?: "${meta.id}"
    def binsize    = params.binsize ?: '500'
    def windowsize = params.windowsize ?: '1000000'
    

    """
    cat <<-END > smoother.cfg
    chrom    =${chromsizes}                  # Chromosomes names/sizes
    profPath =./profiles                     # path to binary profiles (output for prepare)
    trackPath=.
    resPath=./res
    log=log
    #============ Prepare parameters
    bin=${binsize}                           # step for bynary profile
    #============ Statistics parameters
    wSize      =${windowsize}                # size of window (nucleotides)
    flankSize  =10000                        # size of flanks(nucleotides)
    kernelSigma=3000                         # kernel width (nucleotides)
    kernelType =NORMAL                       # type of the kernel: NORM | LEFT_EXP | RIGHT_EXP
    BufSize=40000000
    verbose=1   

    END


    # Background model needs trans-contacts (RNA and DNA on different chromosomes)
    # as a proxy for random/non-specific ligation. Single-chromosome genomes
    # (bacteria) can never satisfy \$3 != \$13, so the background is always
    # empty there and Smoother fails with "profile contains only zeros".
    # Fix: fall back to far-apart cis-contacts (same chromosome, but farther
    # apart than normalisation_cis_bg_mindist) as an equivalent background
    # proxy -- at large enough genomic distance, contact frequency decays to
    # the same background level trans-contacts represent. Opt-in only
    # (params.normalisation_cis_bg_mindist defaults to 0 => disabled): with
    # mindist=0 the added clause is always false, so this is byte-for-byte
    # identical to the original line for every organism that doesn't set it.
    # Original line, kept for reference:
    # awk '\$8 == "protein_coding" && \$3 != \$13 { print \$13, \$14, \$15, \$1 }' OFS='\\t' ${voted} > ${voted}.bg.bed
    
    awk -v mindist="${params.normalisation_cis_bg_mindist ?: 0}" '\$8 == "protein_coding" && (\$3 != \$13 || (mindist > 0 && \$3 == \$13 && (\$14 > \$4 + mindist || \$4 > \$14 + mindist))) { print \$13, \$14, \$15, \$1 }' OFS='\\t' ${voted} > ${voted}.bg.bed

    Smoother cfg="smoother.cfg" ${voted}.bg.bed

    {
        # awk 'NR==1 {print \$0, "bg_sm", "N2_raw"}' OFS='\\t' ${voted}
        awk 'NR>1 { print \$13, \$14, \$14 + 1, \$0 }' OFS='\\t' ${voted}  | \\
        bedtools intersect -a - -b ${voted}.bg_sm.bgr -loj | \\
        cut -f4-26,30 | \\
        awk 'BEGIN { OFS="\\t" } NR>1 { if (\$24 != -1) { n2_raw = 1 / (\$24 + 0.5) } else { n2_raw = 0 } print \$0, n2_raw }'
    } > ${prefix}.intersected.N2_raw.bed

    library_size=\$((\$(wc -l < ${voted} ) - 1 ))
    sum_of_weights=\$(awk '{sum += \$NF} END {print sum}' ${prefix}.intersected.N2_raw.bed)

    awk -v ls="\$library_size" -v sr="\$sum_of_weights" 'BEGIN {OFS="\\t"}
        
        NR>1 { N2 = \$NF * (ls / sr);
               print \$0, N2
        }' ${prefix}.intersected.N2_raw.bed > ${prefix}.intersected.N2.bed
    """

}