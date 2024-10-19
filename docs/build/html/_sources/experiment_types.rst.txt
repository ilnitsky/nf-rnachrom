.. _rnachrom_configuration:

RNAchrom Configuration
======================

The RNAchrom workflow uses a configuration file to specify various parameters and options. This document describes the key sections and parameters in the configuration file.
We have provided a list of configurations most suitable for data processing of each experiment type. While core tools and parameters are constant across pipelines to ensure reliability, 
each configuration uniquely tailors settings to align with the specific analysis type. Those options are the ones used while preparing raw files
for `RNA-Chrom <https://rnachrom2.bioinf.fbb.msu.ru>`_  database:

.. toctree::
  :maxdepth: 3
  :caption: Tutorials
  :titlesonly:

  ./examples/iMARGI.ipynb

Below you can explore options in details.

Input Options
-------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Description
   * - ``input``
     - Path to input file(s)
   * - ``exp_type``
     - Type of experiment (e.g. 'grid', 'radicl', 'imargi')
   * - ``procedure``
     - Processing procedure ('old' or 'new')
   * - ``split_by_chromosomes``
     - Whether to split processing by chromosomes (true/false)

Processing Tools
----------------
Those options determine which specific tools are used at each step of the workflow, affecting processing speed and output quality.
All available tools from modules are listed below:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Description
   * - ``dedup_tool``
     - Tool for deduplication (options: "fastq-dupaway", "fastuniq", "climpify")
   * - ``trim_tool``
     - Tool for read trimming (options: "trimmomatic", "bbduk", "cutadapt", "fastp")
   * - ``align_tool``
     - Tool for alignment (options: "hisat2", "bowtie2", "star")
   * - ``merge_pairedend_tool``
     - Tool for merging paired-end reads (options: "bbmerge", "pear")

Reference Genome
----------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Description
   * - ``genome``
     - Reference genome identifier ("GRCh38", "GRCm38" "TAIR10",  ... )
   * - ``genome_fasta``
     - Path to reference genome FASTA file
   * - ``hisat2_index``; ``bwa_index``; ``star_index``
     - Path to aligner index files
   * - ``splice_sites``
     - Path to splice sites file (optional)

Annotation Files
----------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Description
   * - ``annot_BED``
     - Path to BED annotation file
   * - ``annot_GTF``
     - Path to GTF annotation file
   * - ``blacklist``
     - Path to blacklist regions file
   * - ``chromsizes``
     - Path to chromosome sizes file
   * - ``detect_strand_genes_list``
     - Path to gene list for strand detection

Restriction Sites & Bridge Processing
-------------------------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Description
   * - ``dna_part_processing``
     - Pattern for DNA part processing
   * - ``rna_part_processing``
     - Pattern for RNA part processing
   * - ``bridge_processing``
     - Enable bridge processing (true/false)
   * - ``debridge_tool``
     - Tool for debridge processing (options: "bitap", "chartools")
   * - ``forward_bridge_seq``
     - Forward bridge sequence
   * - ``reverse_bridge_seq``
     - Reverse bridge sequence
   * - ``max_mismatches``
     - Maximum allowed mismatches
   * - ``min_rna_dna_parts_length``
     - Minimum length of RNA/DNA parts
   * - ``description_sequence``
     - Description of sequence processing

Other Options
-------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Description
   * - ``smartseq_filter``
     - Enable Smart-seq filtering (true/false)
   * - ``max_memory``
     - Maximum memory allocation
   * - ``max_cpus``
     - Maximum CPU cores to use
   * - ``max_time``
     - Maximum run time
   * - ``outdir``
     - Output directory path
   * - ``publish_dir_mode``
     - Mode for publishing output files
   * - ``email``
     - Email for notifications

Command Flags
-------------

The configuration file includes specific flags for various tools:

- ``fastq_dupaway``: Flags for fastq-dupaway
- ``trimmomatic``: Flags for Trimmomatic
- ``fastp``: Flags for fastp
- ``pear``: Flags for PEAR
- ``bam_filter``: Flags for BAM filtering

Environment Variables
---------------------

The config file also sets various environment variables like ``PYTHONPATH``, ``R_PROFILE_USER``, etc.

Included Configs
----------------

The main config file includes additional config files:

- ``base.config``: Base configuration for all pipelines
- ``modules.config``: DSL2 module specific options

This configuration setup allows flexible customization of the RNAchrom workflow for different experimental designs and processing requirements.