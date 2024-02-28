// Check for fastq compressed inputs

    shell:
    ext_file = f_ext.getExtension()
    base_name_file = f_ext.getBaseName()
    if (ext_file == "fasta" | ext_file == "fa"){
	    '''
	    cp -n !{f_ext} !{base_name_file}.fasta
	    '''
    }else if(ext_file == "zip"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
	    '''
      gunzip -f -S .zip !{f_ext}
	    cp -n !{old_base_name_file} !{base_name_file}.fasta
	    '''
    }else if(ext_file == "gz"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
	    '''
	    gunzip -f !{f_ext}
	    cp -n !{old_base_name_file} !{base_name_file}.fasta
	    '''
    }else {
      '''
      echo "Your pathogen genome files appear to have the wrong extension. \n Currently, the pipeline only supports .fasta or .fa, or compressed files with .zip or .gz extensions."
      '''
    }
}






log.info print_green("LncPipe Pipeline Complete!")


workflow.onComplete {

        log.info print_green("LncPipe Pipeline Complete!")

        //email information
        if (params.mail) {
            recipient = params.mail
            def subject = 'My LncPipe execution'

            ['mail', '-s', subject, recipient].execute() <<
                    """

            LncPipe execution summary
            ---------------------------
            Your command line: ${workflow.commandLine}
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            Error report: ${workflow.errorReport ?: '-'}
        
            """
        }


}

workflow.onError {
    println print_yellow("Oops... Pipeline execution stopped with the following message: ")+print_red(workflow.errorMessage)
}


def helpMessage() {
    log.info nfcoreHeader()
    log.info "foo"
}

def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}


//help information
params.help = null
if (params.help) {
    log.info ''
    log.info print_purple('------------------------------------------------------------------------')
    log.info "LncPipe: a Nextflow-based Long non-coding RNA analysis Pipeline v$version"
    log.info "LncPipe integrates several NGS processing tools to identify novel long non-coding RNAs from"
    log.info "un-processed RNA sequencing data. To run this pipeline, users either need to install required tools manually"
    log.info "or use the docker image for LncPipe that comes with all tools pre-installed. (note: docker needs to be installed on your system). More information on usage can be found at https://github.com/likelet/LncPipe ."
    log.info "Bugs or new feature requests can be reported by opening issues in our github repository."
    log.info print_purple('------------------------------------------------------------------------')
    log.info ''
    log.info print_yellow('Usage: ')
    log.info print_yellow('    The typical command for running the pipeline is as follows (we do not recommend users passing configuration parameters through command line, please modify the config.file instead):\n') +
            print_purple('       Nextflow run LncRNAanalysisPipe.nf \n') +

            print_yellow('    General arguments:             Input and output setting\n') +
            print_cyan('      --inputdir <path>         ') + print_green('Path to input data(optional), current path default\n') +
            print_cyan('      --reads <*_fq.gz>         ') + print_green('Filename pattern for pairing raw reads, e.g: *_{1,2}.fastq.gz for paired reads\n') +
            print_cyan('      --out_folder <path>           ') + print_green('The output directory where the results will be saved(optional), current path is default\n') +
            print_cyan('      --aligner <hisat>             ') + print_green('Aligner for reads mapping (optional),"hisat"(defalt)/"star"/"tophat"\n') +
            print_cyan('      --qctools <fastp>            ') + print_green('Tools for assess reads quality, fastp(default)/afterqc/fastqc/none(skip QC step)\n') +
            print_cyan('      --detools <edger>             ') + print_green('Tools for differential analysis, edger(default)/deseq/noiseq\n') +
            print_cyan('      --quant <kallisto>            ') + print_green('Tools for estimating abundance of transcript, kallisto(default)/htseq\n') +
            '\n' +
            print_yellow('    Options:                         General options for run this pipeline\n') +
            print_cyan('      --merged_gtf <gtffile>        ') + print_green('Start analysis with assemblies already produced and skip fastqc/alignment step, DEFAOUL NULL\n') +
            print_cyan('      --design <file>               ') + print_green('A flat file stored the experimental design information ( required when perform differential expression analysis)\n') +
            print_cyan('      --singleEnd                   ') + print_green('Reads type, True for single ended \n') +
            print_cyan('      --unstrand                    ') + print_green('RNA library construction strategy, specified for \'unstranded\' library \n') +
            '\n' +
            print_yellow('    References:                      If not specified in the configuration file or you wish to overwrite any of the references.\n') +
            print_cyan('      --fasta                       ') + print_green('Path to Fasta reference(required)\n') +
            print_cyan('      --gencode_annotation_gtf      ') + print_green('An annotation file from GENCODE database in GTF format (required)\n') +
            print_cyan('      --lncipedia_gtf               ') + print_green('An annotation file from LNCipedia database in GTF format (required)\n') +
            '\n' +
            print_yellow('    LncPipeReporter Options:         LncPipeReporter setting  \n') +
            print_cyan('      --lncRep_Output                ') + print_green('Specify report file name, \"report.html\" default.\n') +
            print_cyan('      --lncRep_theme                 ') + print_green('Plot theme setting in interactive plot, \"npg\" default.\n') +
            print_cyan('      --lncRep_min_expressed_sample  ') + print_green('Minimum expressed gene allowed in each sample, 50 default.\n') +
            '\n' +
            print_yellow('    Other options:                   Specify the email and \n') +
            print_cyan('      --sam_processor                ') + print_green('program to process samfile generated by hisat2 if aligner is hisat2. Default \"sambamba\". \n') +
            print_cyan('      --mail                         ') + print_green('email info for reporting status of your LncPipe execution  \n') +



            log.info '------------------------------------------------------------------------'
    log.info print_yellow('Contact information: zhaoqi@sysucc.org.cn')
    log.info print_yellow('Copyright (c) 2013-2017, Sun Yat-sen University Cancer Center.')
    log.info '------------------------------------------------------------------------'
    exit 0
}



// read file
fasta = file(params.fasta)
if (!fasta.exists()) exit 1, "Reference genome not found: ${params.fasta}"
if(params.aligner=='star'){
    star_index = file(params.star_index)
    if (!star_index.exists()) exit 1, "STAR index not found: ${params.star_index}"
}else if(params.aligner =='hisat'){
    hisat2_index = Channel.fromPath("${params.hisat2_index}*")
            .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
}else if(params.aligner =='tophat'){
    bowtie2_index = Channel.fromPath("${params.bowtie2_index}*")
            .ifEmpty { exit 1, "bowtie2 index for tophat not found: ${params.bowtie2_index}" }
}


if (params.species=="human") {
    gencode_annotation_gtf = file(params.gencode_annotation_gtf)
    if (!gencode_annotation_gtf.exists()) exit 1, "GENCODE annotation file not found: ${params.gencode_annotation_gtf}"
    lncipedia_gtf = file(params.lncipedia_gtf)
    if (!lncipedia_gtf.exists()) exit 1, "lncipedia annotation file not found: ${params.lncipedia_gtf}"





    if (params.aligner == 'star' && params.star_index == false && fasta) {
        process Make_STARindex {
            tag fasta

            storeDir { params.outdir + "/STARIndex" }

            input:
            file fasta from fasta
            file gencode_annotation_gtf

            output:
            file "star_index" into star_index

            shell:
            star_threads = ava_cpu- 1
            """
                mkdir star_index
                STAR \
                    --runMode genomeGenerate \
                    --runThreadN ${star_threads} \
                    --sjdbGTFfile $gencode_annotation_gtf \
                    --sjdbOverhang 149 \
                    --genomeDir star_index/ \
                    --genomeFastaFiles $fasta
                """
        }
    } else if (params.aligner == 'star' && params.star_index == false && !fasta) {
        println print_red("No reference fasta sequence loaded! please specify ") + print_red("--fasta") + print_red(" with reference.")

    } else if (params.aligner == 'hisat' && !fasta) {
        process Make_hisat_index {

            tag fasta

            storeDir { params.outdir + "/hisatIndex" }

            input:
            file fasta from fasta
            file gencode_annotation_gtf

            output:
            file "genome_ht2.*" into hisat2_index

            shell:
            hisat2_index_threads = ava_cpu- 1
            """
                #for human genome it will take more than 160GB memory and take really  long time (6 more hours), thus we recommand to down pre-build genome from hisat website
                extract_splice_sites.py !{gencode_annotation_gtf} >genome_ht2.ss
                extract_exons.py !{gencode_annotation_gtf} > genome_ht2.exon
                hisat2-build -p !{hisat2_index_threads} --ss genome_ht2.ss --exo genome_ht2.exon !{fasta} genome_ht2
                """
        }
    } else if (params.aligner == 'tophat' && params.hisat_index == false && !fasta) {
        println print_red("No reference fasta sequence loaded! please specify ") + print_red("--fasta") + print_red(" with reference.")
    }

trimmomatic_params: 'SLIDINGWINDOW:5:26 MINLEN:0'

include { TRIMMOMATIC as FASTQ_TRIM } from '../../modules/local/trimmomatic/main' addParams( options: [args: trimmomatic_params, suffix:'.trim', args2: [gzip: false]])











def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/dualrnaseq --input '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --input [file]                  Path to input data (must be surrounded with quotes)
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, podman, awsbatch, <institute> and more

    References and annotative files can be specified in the configuration file.
    Alternatively, the following params can be edited directly.
    
    Library type and genome files:
      --single_end      [bool]  Specifies that the input is single-end reads (default: false)
      --fasta_host      [file]  Host genome ("folder/file.fa")
      --fasta_pathogen  [file]  Pathogen genome
      
    Annotation files:
      --gff_host        [file]   Host GFF
      --gff_pathogen    [file]   Pathogen GFF
      --gff_host_tRNA   [file]   Host tRNA (optional)

    The pipeline will automatically generate transcriptome files for both the host and pathogen.
    These parameters should only be used when using custom transcriptome files 

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode, [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    Transcriptome files:
      --read_transcriptome_fasta_host_from_file        [bool]   Include custom host transcriptome
                                                                (Default: false)
      --read_transcriptome_fasta_pathogen_from_file    [bool]   Include custom pathogen transcriptome
                                                                (Default: false)
      --transcriptome_host                             [file]    Custom host transcriptome
                                                                (Default: "")
      --transcriptome_pathogen                         [file]    Custom pathogen transcriptome
                                                                (Default: "")  
      
    Trimming is performed by either Cutadapt or BBDuk with the following related options
    
    Cutadapt:
      --run_cutadapt    [bool]  To run cutadapt (Default: false)
      --a               [str]   Adapter sequence for single-end reads or first reads of paired-end data
                                (Default: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
      --A               [str]   Adapter sequence for second reads of paired-end data
                                (Default: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
      --quality-cutoff  [int]   Cutoff to remove low-quality ends of reads. (Default: 10)
                                A single cutoff value is used to trim the 3’ end of reads. 
                                If two comma-separated cutoffs are defined, the first value reprerents 5’ cutoff, 
                                and the second value defines the 3’ cutoff.
      --cutadapt_params [str]   Set of additional cutadapt parameters
    
    BBDuk:
      --run_bbduk    [bool]    To run BBDuk (Default: false)
      --minlen       [int]     Reads shorter than this after trimming will be discarded.
                               Pairs will be discarded if both are shorter.
                               (Default: 18)
      --qtrim        [str]     To trim read ends to remove bases with quality below trimq. 
                               Possible options:rl (trim both ends), f (neither end), r (right end only), 
                               l (left end only), w (sliding window). 
                               (Default: "r")
      --trimq        [int]     Cutoff to trim regions with average quality BELOW given value. 
                               Option is avaiable if qtrim is set to something other than f. 
                               (Default: 10)
      --ktrim        [str]     To trim reads to remove bases matching reference kmers. 
                               Avaiable options: f (don't trim), r (trim to the right - 3' adapters), 
                               l (trim to the left - 5' adapters). 
                               (Default: "r")
      --k            [int]     Kmer length used for finding contaminants (adapters). Contaminants 
                               shorter than k will not be found. k must be at least 1.
                               (Default: 17)
      --mink        [int]      Look for shorter kmers at read tips down to this length, 
                               when k-trimming or masking. 0 means disabled.  Enabling
                               this will disable maskmiddle.
                               (Default: 11)
      --hdist        [int]     Maximum Hamming distance for ref kmers (subs only).
                               (Default: 1)
      --adapters     [file]    Fasta file with adapter sequences (Default: $projectDir/data/adapters.fa)
      --bbduk_params [str]     Set of additional BBDuk parameters 
    
    Basic quality control is reported through FastQC, which is run on raw reads and trimmed reads.
    
    FastQC:
      --skip_fastqc      [bool]  Option to skip running FastQC
      --fastqc_params   [str]   Set of additional fastqc parameters
      
    The following options are related to the three main methods to extract gene expression:
    
    Salmon:
      --libtype                     [str]   To define the type of sequencing library of your data 
                                            (Default:'')
      --kmer_length                 [int]   To define the k-mer length (-k parameter in Salmon)
                                            (Default: 21)
      --writeUnmappedNames          [bool]  By default the pipeline does not save names of unmapped reads
                                            (Default: false)
      --softclipOverhangs           [bool]  By default, the pipeline does not allow soft-clipping of reads 
                                            (Default: false)
      --incompatPrior               [int]   This is set to 0.0, to ensure that only mappings or alignments that 
                                            are compatible with the specified library type are considered by Salmon 
                                            (Default: 0.0)
      --dumpEq                      [bool]  To save the equivalence classes and their counts, change this option to True 
                                            (Default: false)
      --writeMappings               [bool]  If set to True, the pipeline will create a files named mapping.sam 
                                            containing mapping information
                                            (Default: false)
      --keepDuplicates              [bool]  Option to remove/collapse identical transcripts during the indexing stage 
                                            (Default: false)
      --generate_salmon_uniq_ambig  [bool]  Option to extract all of the unique and ambigious reads after quantification
                                            Works for both Selective alignment and alignment-based modes 
                                            (Default: false)					      
      
    Salmon selective alignment:
      --run_salmon_selective_alignment                      [bool]   Run this mode
                                                                     (Default: false)
      --gene_attribute_gff_to_create_transcriptome_host     [str]    Host transcriptome  - gene attributes
                                                                     (Default: transcript_id) 
      --gene_feature_gff_to_create_transcriptome_host       [str]    Host transcriptome  - gene feature
                                                                     (Default: ["exon", "tRNA"])
      --gene_attribute_gff_to_create_transcriptome_pathogen [str]    Pathogen transcriptome - gene attribute
                                                                     (Default: locus_tag)
      --gene_feature_gff_to_create_transcriptome_pathogen   [str]    Pathogen transcriptome - gene features
                                                                     (Default: ["gene","sRNA","tRNA","rRNA"] )
      --salmon_sa_params_index                              [str]    Set of additional parameters for creating an index with Salmon Selective Alignment
      --salmon_sa_params_mapping                            [str]    Set of additional parameters for mapping with Salmon Selective Alignment
      
    STAR - general options available for both modes, genome mapping with HTSeq quantification and salmon - alignment-based mode:
      --run_star                          [bool]   Run STAR
                                                   (Default: false)
      --outSAMunmapped                    [str]    By default, the pipeline saves unmapped reads 
                                                   within the main BAM file. If you want to switch off this option, 
                                                   set the --outSAMunmapped flag to None
                                                   (Default: Within)
      --outSAMattributes                  [str]    To specify the attributes of the output BAm file
                                                   (Default: Standard)
      --outFilterMultimapNmax             [int]    To specify the maximum number of loci a read is allowed to map to
                                                   (Default: 999)
      --outFilterType                     [str]    By default, the pipeline keeps reads containing junctions that 
                                                   passed filtering into the file SJ.out.tab. This option reduces 
                                                   the number of ”spurious” junctions
                                                   (Default: BySJout)
      --alignSJoverhangMin                [int]    The number of minimum overhang for unannotated junctions
                                                   (Default: 8)
      --alignSJDBoverhangMin              [int]    The number of minimum overhang for annotated junctions
                                                   (Default: 1)
      --outFilterMismatchNmax             [int]    To define a threshold for the number of mismatches to be allowed.
                                                   The pipeline uses a large number to switch this filter off 
                                                   (Default: 999)
      --outFilterMismatchNoverReadLmax    [int]    Here, you can define a threshold for a ratio of mismatches to 
                                                   read length. The alignment will be considered if the ratio is 
                                                   less than or equal to this value. For 2x100b, max number of 
                                                   mismatches is 0.04x200=8 for paired-end reads
                                                   (Default: 1)
      --alignIntronMin                    [int]    By default, the nf-core dualrnaseq pipeline uses 20 as a 
                                                   minimum intron length. If the genomic gap is smaller than this
                                                   value, it is considered as a deletion
                                                   (Default: 20)
      --alignIntronMax                    [int]    The maximum intron length
                                                   (Default: 1000000)
      --alignMatesGapMax                  [int]    The maximum genomic distance between mates is 1,000,000
                                                   (Default: 1000000)
      --limitBAMsortRAM                   [int]    Option to limit RAM when sorting BAM file. 
                                                   If 0, will be set to the genome index size, which can be quite 
                                                   large when running on a desktop or laptop
                                                   (Default: 0)
      --winAnchorMultimapNmax             [int]    By default, the nf-core dualrnaseq pipeline uses 999 as a 
                                                   maximum number of loci anchors that are allowed to map to
                                                   (Default: 999)
      --sjdbOverhang                      [int]    To specify the length of the donor/acceptor sequence on each side of the junctions
                                                   used in constructing the splice junctions database.
                                                   ideally = (mate length - 1)
                                                   (Default: 100)
      
    STAR - additional options available only for genome mapping with HTSeq quantification mode
      --outWigType                        [str]    Used to generate signal outputs, such as "wiggle" and "bedGraph"
                                                   (Default: None)
      --outWigStrand                      [str]    Options are Stranded or Unstranded when defining 
                                                   the strandedness of wiggle/bedGraph output
                                                   (Default: Stranded)
      --star_index_params                 [str]    Set of additional parameters for creating an index with STAR
      --star_alignment_params             [str]    Set of additional parameters for alignment with STAR
      
    STAR - additional options available only for Salmon - alignment-based mode:
      --quantTranscriptomeBan             [str]    The nf-core/dualrnaseq pipeline runs STAR to generate a 
                                                   transcriptomic alignments. By default, it allows for insertions, 
                                                   deletions and soft-clips (Singleend option). To prohibit this 
                                                   behaviour, specify IndelSoftclipSingleend
                                                   (Default: Singleend)
      --star_salmon_index_params          [str]    Set of additional parameters for creating an index with STAR in salmon alignment-based mode
      --star_salmon_alignment_params      [str]    Set of additional parameters for alignment with STAR in salmon alignment-based mode
      
    Salmon - alignment-based mode:
      --run_salmon_alignment_based_mode   [bool]   Option to run Salmn in alignment mode
                                                   (Default: false)
      -- salmon_alignment_based_params    [str]    Set of additional parameters for Salmon in alignment-based mode
      
    HTSeq:
      --run_htseq_uniquely_mapped               [bool]   Option to run HTSeq
                                                         (Default: false)
      --stranded                                [char]   Is library type stranded (yes/no)
                                                         (Default: yes)
      --max_reads_in_buffer                     [int]    To define the number of maximum reads allowed to
                                                         stay in memory until the mates are found
                                                         Has an effect for paired-end reads
                                                         (Default: 30000000) 
      --minaqual                                [int]    To specify a threshold for a minimal MAPQ alignment quality.(Default: 10) 
                                                         Reads with the MAPQ alignment quality below the given number will be removed
                                                         (Default: 10)
      --htseq_params                            [str]    Set of additional parameters for HTSeq
      --gene_feature_gff_to_quantify_host       [str]    Host - gene features to quantify from GFF
                                                         (Default: ["exon","tRNA"] )
      --gene_feature_gff_to_quantify_pathogen   [str]    Pathogen - gene features to quantify from GFF
                                                         (Default: ["gene", "sRNA", "tRNA", "rRNA"] )
      --host_gff_attribute                      [str]    Host - attribute from GFF
                                                         (Default: gene_id )
      --pathogen_gff_attribute                  [str]    Pathogen - attribute from GFF
                                                         (Default: locus_tag)
     
    RNA mapping statistics:
      --mapping_statistics               [bool]   Option to generate mapping statistics. This will create the following:
                                                  - Count the total number of reads before and after trimming
                                                  - Scatterplots comparing all replicates (separate for both host and pathogen reads)
                                                  - Plots of the % of mapped/quantified reads
                                                  - Plots of RNA-class statistics (as many types can be identified, 
                                                    the parameter below --rna_classes_to_replace_host can help to summarise these)
                                                  (Default: false)
      --rna_classes_to_replace_host      [file]   Located within the data/ directory, this tab delimited file contains headers which 
                                                  groups similar types of RNA classes together. This helps to keep the RNA-class 
                                                  names simplified for plotting purposes.
                                                  (Default: $projectDir/data/RNA_classes_to_replace.tsv)
   
    Report options:
      --email                   [email]   Set this parameter to your e-mail address to get a summary e-mail with details of the 
                                          run sent to you when the workflow exits
                                          (Default: false)
      --email_on_fail           [email]   Same as --email, except only send mail if the workflow is not successful
                                          (Default: false)  
      --max_multiqc_email_size  [str]     Theshold size for MultiQC report to be attached in notification email. 
                                          If file generated by pipeline exceeds the threshold, it will not be attached
                                          (Default: 25MB)
      --plaintext_email         [bool]    Set to receive plain-text e-mails instead of HTML formatted
                                          (Default: false)
      --monochrome_logs         [bool]    Set to disable colourful command line output and live life in monochrome
                                          (Default: false)
      --multiqc_config          [bool]    Specify path to a custom MultiQC configuration file.
                                          (Default: false)

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}