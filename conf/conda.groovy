process {

    withLabel: 'PrepareSoftware' {
        conda = "${projectDir}/envs/rnachromprocessing.yaml"
    }

    withName: 'FASTQ_DUPAWAY' {
        conda = "${projectDir}/envs/rnachromprocessing.yaml"
    }
    
    withName: 'DETECT_STRAND' {
        conda = "${projectDir}/envs/rnachromprocessing.yaml"
    }
    
    withName: 'RSITES' {
        conda = "${projectDir}/envs/rnachromprocessing.yaml"
    }


    withName: 'HISAT2_BUILD' {
        conda = "${projectDir}/envs/primary_processing.yml"
    }



    withName: 'FASTUNIQ' {
        conda = "${projectDir}/envs/primary_processing.yml"
    }

    withName: 'BBMAP_CLUMPIFY' {
        conda = "${projectDir}/envs/primary_processing.yml"
    }


    
    withName: 'TRIMMOMATIC' {
        conda = "${projectDir}/envs/primary_processing.yml"
    }
    
    withName: 'BBMAP_BBDUK' {
        conda = "${projectDir}/envs/primary_processing.yml"
    }


    withName: 'SAMTOOLS_VIEW' {
        conda = "${projectDir}/envs/primary_processing.yml"
    }

    withName: 'BEDTOOLS_BAMTOBED' {
        conda = "${projectDir}/envs/primary_processing.yml"
    }

    withName: 'CIGAR_FILTER' {
        conda = "${projectDir}/envs/secondary_processing.yml"
    }

    withName: 'MERGE_REPLICAS' {
        conda = "${projectDir}/envs/secondary_processing.yml"
    }
}