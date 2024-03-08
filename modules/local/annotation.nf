
process ANNOTATION_VOTING {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$cigar.baseName"
    // errorStrategy 'ignore'

    publishDir (
        path: { "$params.outdir/annotation" },
        mode: "copy"
    ) 
    input:
    tuple val(id), path(cigar)

    output:
    tuple val(id), path("voted/*.tab"),            emit: voted
    tuple val(id), path("singletons/*.tab"),       emit: singletons
    tuple val(id), path("complement_annot/*.tab"), emit: complement_annot
    tuple val(id), path("selected_annot/*.tab"),   emit: selected_annot

    script:
    def procedure = ''
    procedure = params.procedure == 'new' ? '--rna_parts' : ''

    """
    python3 ${projectDir}/bin/rnachrom_pipeline_faster/annotation_voting.py \\
    ${cigar} ${params.annot_BED} ${procedure} --no_stat --cpus $task.cpus --outdir .

    mkdir voted singletons selected_annot complement_annot
    mv *${id}.4-voted.tab  voted/
    mv *${id}.3-complement_annot.tab  complement_annot/
    mv *${id}.3-singletons.tab  singletons/
    mv *${id}.3-selected_annot.tab  selected_annot/

    """
    

}


    // wc -l ${id}.4-voted.tab | cut -f1 -d' ' > ${id}.voted.stat
    // wc -l ${id}.3-singletons.tab | cut -f1 -d' ' > ${id}.singletons.stat