process NUCL_DISTR_RSITES {
//TODO:  folders
    conda "${projectDir}/envs/secondary_processing.yml"

    publishDir (
        path: { "$params.outdir/nucleotide_distributions" },
        mode: "copy"
        // pattern: "dedup/*.fastq",
        // saveAs: { fn -> file(fn).name } 
    )

    input:
    tuple val(meta), path(reads)

    output:
    path '*.png'

    script:
    """
    python3 ${projectDir}/bin/LastKNuclPlot.py ${reads} 2 f
    python3 ${projectDir}/bin/LastKNuclPlot.py ${reads} 3 f
    python3 ${projectDir}/bin/LastKNuclPlot.py ${reads} 2 r
    python3 ${projectDir}/bin/LastKNuclPlot.py ${reads} 3 r
    """

}

// bioawk -c fastx '{print $seq}' ~/nf-rnachrom/data/imargi/SRR9900120_2.fastq | head -n200000  | cut -c1-2 | sort | uniq -c