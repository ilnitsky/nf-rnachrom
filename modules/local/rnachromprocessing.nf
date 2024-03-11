
// conda "conda-forge::gawk=5.1.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
    //     'quay.io/biocontainers/gawk:5.1.0' }"

ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";


def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }



process CONTACTS {
    tag "$params.exp_type"
    label 'process_single'
    
    publishDir (
        path: { "$params.outdir/contacts" },
        mode: "copy",
        pattern: "contacts/*",
        saveAs: { fn -> file(fn).name } 
    )

    input:
    path(beds)
    path(config)

    output:
    path 'contacts/*', emit: contacts


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    rnachromprocessing -c config.final.json -s contacts --input_dir . --output_dir . -v
    """

}











process NORMALIZE_RAW {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$voted.baseName"
    publishDir (
        path: { "$params.outdir/Normalize_raw" },
        mode: "copy"
    ) 
   // memory '15 GB'

    input:
    tuple val(name), path(voted), path(background)

    output:
    tuple val(name), path("*.5-N2_raw.tab"), emit: raw_norm
    tuple val(name), path("*.5-N2_raw.stat.tab"), emit: raw_stat

    script:
    """
    python3 ${projectDir}/bin/rnachrom_pipeline_faster/normalize_raw.py \\
    ${voted} ${background} --outdir .
    """

}

process NORMALIZE_N2 {
    //TODO: Add combining stats into single file 
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$norm_raw.baseName"
    publishDir (
        path: { "$params.outdir/Normalize_N2" },
        mode: "copy"
    ) 
   // memory '15 GB'

    input:
    tuple val(name), path(norm_raw_stat), path(norm_raw)

    output:
    tuple val(name), path("*.5-N2.tab")

    script:

    """
    python3 ${projectDir}/bin/rnachrom_pipeline_faster/normalize_N2.py \\
    ${norm_raw} ${params.chromsizes} --by_chr ${norm_raw_stat} --outdir .
    """

}

process SCALING {
    conda "${projectDir}/envs/secondary_processing.yml"
    //TODO: Collect all scaling files by chr
    tag "$norm_n2.baseName"
    publishDir (
        path: { "$params.outdir/Scaling" },
        mode: "copy"
    ) 
   // memory '15 GB'

    input:
    tuple val(name), path(norm_n2)

    output:
    tuple val(name), path("*")

    script:

    """
    python3 ${projectDir}/bin/rnachrom_pipeline_faster/scaling.py \\
    ${norm_n2} ${params.chromsizes} --outdir .
    """

}

process VALIDATE_ANNOT {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$annot.baseName"
    // Remove header for Bardic
    // TODO: preparation of BED or GTF files

    publishDir (
        path: { "$params.outdir/Validated_annot" },
        mode: "copy"
    ) 
   // memory '15 GB'

    input:
    path(annot)

    output:
    path("*")

    script:

    """
    echo \$CONDA_PREFIX
    python3 ${projectDir}/bin/rnachrom_pipeline_faster/validate_annotation.py \\
    ${annot} ${params.chromsizes} --annot_format BED --outdir .

    sed -i '1d' ${annot.baseName}.0-corrected_annot.bed6
    """

}




// process BLACKLISTING {
//   publishDir "$params.outdir/blkres", mode: 'copy'
//   input:
//   path(splice)

//   output:
//   path("*-blacklist.tab")

//   script:
//   """
//   python3 ${projectDir}/bin/rnachrom_pipeline_faster/blacklisting.py ${splice} ${params.blacklist} --outdir .
//   """
// }




//     script:
//     """
//     jq \\
//     --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[1]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
//     --argjson dna_ids "\$(jq -nR '[inputs|split(",")|.[2]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
//      '. + {rna_ids: \$rna_ids, dna_ids: \$dna_ids}' $config > config.json
    

//     jq '.align.dna_genome_path="${params.genome_path}"' config.json \\
//     | jq '.align.known_splice="${params.splice_sites}"' \\
//     | jq '.align.tool="${params.align_tool}"' \\
//     | jq '.rsites.type="${params.exp_type}"' \\
//     | jq '.trim.tool="${params.trim_tool}"' \\
//     | jq '.trim.tool_path="${params.trim_path}"' \\
//     > config.final.json






// process CONFIG {
//     // publishDir ${params.outdir}/'processing', mode: 'copy'
//     // log.info print_green("Complete!")
//     log.info print_yellow("Experiment:           ") + print_purple(params.exp_type)

//     label 'process_single'

    
//     input:
//     path samplesheet
//     path config

//     output:
//     path '*.final.json', emit: json


//     when:
//     task.ext.when == null || task.ext.when


    
//     script:

//     if (params.exp_type in ['radicl', 'char', 'redc', 'margi']) {
//         rsites = 'skip'
//     } else {
//         rsites = params.exp_type
//     }
//     if (!params.bridge_processing){    
//     """
//     jq \\
//     --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[1]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
//     --argjson dna_ids "\$(jq -nR '[inputs|split(",")|.[2]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
//      '. + {rna_ids: \$rna_ids, dna_ids: \$dna_ids}' $config > config.json
    
//     jq '.align.dna_genome_path="${params.genome_path}"' config.json \\
//     | jq '.align.known_splice="${params.splice_sites}"' \\
//     | jq '.align.tool="${params.align_tool}"' \\
//     | jq '.align.tool_path="${params.align_tool_path}"' \\
//     | jq '.rsites.type="${rsites}"' \\
//     | jq '.trim.tool="${params.trim_tool}"' \\
//     | jq '.trim.tool_path="${params.trim_path}"' \\
//     | jq '.dedup.tool="${params.dedup_tool}"' \\
//     | jq '.dedup.tool_path="${params.dedup_tool_path}"' \\
//     > config.final.json
//     """
//     } else {
//     """
//     jq \\
//     --argjson rna_ids "\$(jq -nR '[inputs|split(",")|.[2]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
//     --argjson dna_ids "\$(jq -nR '[inputs|split(",")|.[3]|split("/")|.[-1]|split(".")|.[0]]' <(tail -n+2 $samplesheet))"   \\
//      '. + {rna_ids: \$rna_ids, dna_ids: \$dna_ids}' $config > config.json

//     jq '.align.dna_genome_path="${params.genome_path}"' config.json \\
//     | jq '.align.known_splice="${params.splice_sites}"' \\
//     | jq '.align.tool="${params.align_tool}"' \\
//     | jq '.align.tool_path="${params.align_tool_path}"' \\
//     | jq '.rsites.type="skip"' \\
//     | jq '.trim.tool="skip"' \\
//     | jq '.trim.tool_path="skip"' \\
//     | jq '.dedup.tool="skip"' \\
//     | jq '.dedup.tool_path="skip"' \\
//     > config.final.json
//     """
//     }

// }

// process DEDUP {

//     tag "$params.dedup_tool"
    
//     publishDir (
//         path: { "$params.outdir/dedup" },
//         mode: "copy",
//         pattern: "dedup/*.fastq",
//         saveAs: { fn -> file(fn).name } 
//     )

//     tag "$params.exp_type"
//     // label 'process_single'
        
//     input:
//     path reads
//     path config 

//     output:
//     path 'dedup/*.fastq', emit: dedup


//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     """
//     rnachromprocessing -c ${config} -s dedup --input_dir . --output_dir . -v
//     """


// }

// process RSITES {
//     // withLabel : x
//     // executor "local"
//     // maxForks 48
//     // cpus 1
//     // errorStrategy 'ignore'
//     // executor = 'slurm'
//     // beforeScript 'module load cuda/11.4.2; export SINGULARITY_NV=1'
//     // errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
//     // echo ${task.memory}
//     // echo ${task.time}
//     // echo ${task.cpus}


    
//     publishDir (
//         path: { "$params.outdir/rsites" },
//         mode: "copy",
//         pattern: "rsites/*.fastq",
//         saveAs: { fn -> file(fn).name } 
//     )

//     tag "$params.exp_type"
//     label 'process_single'
        
//     input:
//     path(reads)
//     path(config)

//     output:
//     path 'rsites/*.fastq', emit: rsites


//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     """
//     rnachromprocessing -c config.final.json -s rsites --input_dir . --output_dir . -v
//     """


// }



// process TRIM {
//     tag "$params.trim_tool"
//     label 'process_single'

//     publishDir (
//         path: { "$params.outdir/trim" },
//         mode: "copy",
//         pattern: "trim/*.fastq",
//         saveAs: { fn -> file(fn).name } 
//     )
        
//     input:
//     path(reads)
//     path(config)

//     output:
//     path 'trim/*', emit: trim


//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     if (params.trim_tool == 'trimmomatic') {
//     """
//     rnachromprocessing -c config.final.json -s trim --input_dir . --output_dir . -v
//     """
//     } else {
    

//     }

// }

// process ALIGN {
//     tag "$params.genome_path"
//     // label 'process_single'


//     publishDir (
//         path: { "$params.outdir/align" },
//         mode: "copy",
//         pattern: "align/*",
//         saveAs: { fn -> file(fn).name } 
//     )
        
//     input:
//     path(trimmed)
//     path(config)

//     output:
//     path('align/*'), emit: align


//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     println print_purple("Aligning files with " + params.align_tool + " tool" )
//     """
//     which python
//     rnachromprocessing -c config.final.json -s align --input_dir . --output_dir . -v
//     """

// }

// process BAM_FILTER {
//     tag "$config"
//     label 'process_single'

//     publishDir (
//         path: { "$params.outdir/bam" },
//         mode: "copy",
//         pattern: "bam/*",
//         saveAs: { fn -> file(fn).name } 
//     )
        
//     input:
//     path(aligned)
//     path(config)

//     output:
//     path 'bam/*.bam', emit: bam


//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     """
//     rnachromprocessing -c config.final.json -s bam --input_dir . --output_dir . -v
//     """

// }


// process BED_FILES {
//     tag "$config"
//     label 'process_single'

//     publishDir (
//         path: { "$params.outdir/bed" },
//         mode: "copy",
//         pattern: "bed/*",
//         saveAs: { fn -> file(fn).name } 
//     )  
//     input:
//     path(bams)
//     path(config)

//     output:
//     path 'bed/*', emit: bed


//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     """
//     rnachromprocessing -c config.final.json -s bed --input_dir . --output_dir . -v
//     """

// }
