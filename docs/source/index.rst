Overview
========

`RNAchrom` is a comprehensive and flexible Nextflow pipeline designed to process RNA-DNA interactome sequencing data. It efficiently handles large-scale data from various experiments such as GRID-seq, RADICL-seq, and iMARGI, providing a streamlined workflow for analyzing RNA-DNA interactions.

.. figure:: _static/rna-chrom-processing-pipeline.png
   :width: 100%
   :alt: The diagram of a typical processing pipeline for RNA-DNA interactome data
   :align: center

   The typical `RNAchrom` pipeline involves stages such as data trimming, alignment, and annotation. The interactions between RNA and DNA components are extracted, filtered, and then visualized or further analyzed.

Key Features
------------

- **Input Management**: Automatically validates input files and stages them for further processing.
- **Modular Design**: Includes flexible modules for deduplication, trimming, alignment, and more.
- **Comprehensive Analysis**: Handles both one-to-all (e.g., RAP) and all-to-all (e.g., GRID-seq) experimental designs.
- **Quality Control**: Integrates with FastQC and MultiQC for comprehensive quality assessments.
- **Flexible Customization**: Users can skip stages or customize configurations specific to their experimental setup.

Pipeline Modules
----------------

`RNAchrom` includes a variety of modules that can be utilized as needed in the workflow --
see :doc:`modules_reference` for the full per-module reference:

==================== ==============================================
Module               Description
==================== ==============================================
SAMPLESHEET_CHECK    Validates and stages input samples.
DEDUP                Deduplicates sequencing reads (fastq-dupaway / fastuniq / clumpify / seqkit_rmdup).
FASTP / TRIM         Trims adapters and low-quality bases (fastp / TrimGalore / Trimmomatic).
SMARTSEQ_FILTER      GGG-filter for Red-C raw reads.
PEAR / BBMERGE       Merges overlapping paired-end reads before bridge search.
BITAP_DEBRIDGE       Separates RNA and DNA parts by bridge sequence (also: chartools debridger).
RSITES               Applies restriction-site trimming pattern to RNA/DNA parts.
ALIGN                Aligns RNA and DNA parts to the reference genome (HISAT2 / STAR / BWA-MEM / Bowtie2).
BAM_TO_CONTACTS      Joins aligned RNA and DNA BAMs into raw contact pairs.
FILTER_CONTACTS      Filters contacts by edit distance and CIGAR pattern.
DETECT_STRAND        Votes and corrects RNA-part strand orientation.
ANNOTATION           Annotation voting for RNA parts of contacts (two variants: old / smart).
BLACKLIST            Removes contacts overlapping a blacklist BED.
BARDIC / MACS2       Identifies significant peaks of chromatin-interacting RNAs (ATA / OTA).
NORMALISATION        Background-model normalisation of contacts (ATA).
CHROMATIN_POTENTIAL  Computes chromatin potential per RNA from UU/UM contacts + RNA-seq.
UCARNA_ASSEMBLY      Assembles unannotated chromatin-associated RNAs (ucaRNAs) with StringTie.
MULTIQC              Aggregates results across multiple samples for comparative analysis.
==================== ==============================================

Contents:

.. toctree::
   :hidden:

   self

.. toctree::
   :maxdepth: 3

   installation
   input
   configuration
   stages
   modules_reference
   experiment_types
   results

.. toctree::
  :maxdepth: 3
  :caption: Stages
  :titlesonly:

  ./examples/RNAchrom_pipeline_walkthrough.ipynb

.. toctree::
  :maxdepth: 3
  :caption: Optional Stages
  :titlesonly:

  ./examples/ucaRNA_identififcation.ipynb
  ./examples/Chromatin_potential.ipynb
  
Setup and Configuration
-----------------------

The pipeline can be configured using custom parameter settings to fit the needs of different experimental designs. Key configuration settings include:

- **Genome and Annotation**: Provides support for multiple reference genomes and annotation files.
- **Toolchain Configuration**: Choice of tools for alignment (e.g., HISAT2, STAR) and trimming (e.g., Trimmomatic, FastP).
- **Output Management**: Options for generating summarized reports and logs for comprehensive analysis.

User Support and Community
--------------------------

- **Documentation**: Detailed installation and execution instructions available.
- **Community Support**: Engage with the community through forums and GitHub issues.
- **Contributions**: Open to contributions from the research community to enhance features.

For additional help and support, please check our community forums and our GitHub repository.

* :ref:`genindex`