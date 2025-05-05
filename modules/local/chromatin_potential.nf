process CHROMATIN_POTENTIAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/envs/full_env.yml"

  
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
        
   
    input:
    tuple val(meta), path(contacts_dir)
    path(rnaseq_data)
    path chrom_sizes

    output:
    tuple val(meta), path("*.normalized.tab"), emit: normalized_contacts
    tuple val(meta), path("*.stats.txt"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    
    ${projectDir}/bin/chromatin_potential/count_contacts_all.sh \
                        -d 500000 \
                        -i ${contacts_dir} \
                        -o ./chP \
                        -u contacts.voting.UU.bed \
                        -m contacts.voting.UM.bed \
                        -n 100 \
                        -f 0.05 \
                        -r ${rnaseq_data}
    

    """
}