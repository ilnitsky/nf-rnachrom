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

`RNAchrom` includes a variety of modules that can be utilized as needed in the workflow:

============ ==============================================
Module       Description
============ ==============================================
INPUT_CHECK  Validates and stages input samples.
DEDUP        Deduplicates sequencing reads.
TRIM         Trims sequencing reads to remove adaptors and low-quality bases.
ALIGN        Aligns reads to the reference genome.
ANNOTATE     Conducts annotation voting for the DNA segments.
MACS2        Identifies significant peaks in RNA-DNA interaction data.
MULTIQC      Aggregates results across multiple samples for comparative analysis.
============ ==============================================

Contents:

.. toctree::
   :hidden:

   self

.. toctree::
   :maxdepth: 3

   quickstart
   installation
   configuration

.. toctree::
  :maxdepth: 3
  :caption: Tutorials
  :titlesonly:

  ./examples/RNAchrom_pipeline_walkthrough.ipynb
  
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