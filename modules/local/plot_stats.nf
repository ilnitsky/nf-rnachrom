
process PLOT_STATS {
    conda "${projectDir}/envs/secondary_processing.yml"

    publishDir (
        path: { "$params.outdir/stats" },
        mode: "copy"
    ) 

    input:
    val(table)

    output:
    path("*.png"), emit: png

    script:
    """
    echo 'Starting A' `date +%H-%M-%S`
    echo "${table}" >  stats.tsv
    python ${projectDir}/bin/plot_stats.py
    
    """


}


