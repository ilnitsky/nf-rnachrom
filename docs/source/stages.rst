Stages of  Data Analysis
========================

This pipeline outlines the key stages for analyzing all-to-all interactome sequencing data, focusing on the steps that lead to the final results of RNA-DNA contact pairs, RNA annotation, and significant peaks of chromatin-interacting RNAs.

1. Input and Preprocessing
--------------------------
* Read input data (FASTQ files) and validate sample information
* Perform quality control with FASTQC

2. Deduplication (optional)
---------------------------
* Remove duplicate reads using tools like:

  - fastq-dupaway
  - fastuniq
  - clumpify

3. Trimming
-----------
* Trim low-quality bases and adapters using tools like:

  - fastp
  - Trimmomatic
  - BBduk
  - cutadapt

4. Bridge Processing (for specific experiment types)
----------------------------------------------------

While creating description sequence for experiments that have non-separated raw reads, we
include bridge processing step. 

Based on the bridge sequence, included in the configuration file single-end or paired-end reads are
split into files RNA.fastq and DNA.fastq. Bridge/Linker Sequence in the configuration file
always has the following orientation 5'-{RNA}-{Forward Bridge Sequence}-{DNA}-3'

Bridge processing tools include our own tool based on bitap-search and debridge.jl program based on fuzzy
search from charseq https://github.com/straightlab/chartools/tree/main/Jchartools

* For experiments like GRID-seq, RADICL-seq, iMARGI, etc., process the bridge sequences that connect RNA and DNA parts
* Use tools like BBMerge or PEAR to merge paired-end reads
* Separate RNA and DNA parts based on the bridge sequence

5. Restriction sites filtering 
----------------------------------------------------

- DNA Sequence Start (\*): Begin your DNA sequence with the ``*`` symbol to indicate the start of the DNA part.

- Add Sequence (`+[CATG]`): Use the + operator followed by the sequence you want to add in square brackets. For example, +[CATG]* means you are adding the sequence "CATG" 
to the 5' of the DNA part.




- RNA Sequence Start (`.`): If you need to specify an RNA sequence, use the . symbol to indicate the start of the RNA part. In this case, it seems to be used as an endpoint.

6. Alignment
------------
* Align RNA and DNA reads to the reference genome using tools like:

  - HISAT2
  - STAR
  - BWA-MEM
  - Bowtie2

6. Post-alignment Processing
----------------------------
* Filter aligned reads for uniqueness and mismatches
* Convert BAM files to BED format

7. Contact Generation
---------------------
* Join RNA and DNA parts to create raw contacts
* Perform strand detection and correction

8. CIGAR Filtering (optional)
-----------------------------
* Filter contacts based on CIGAR strings to improve quality

9. Merging Replicates
---------------------
* Combine data from replicate experiments

10. Chromosome Splitting (optional)
-----------------------------------
* Split data by chromosomes for parallel processing

11. Annotation and Voting
-------------------------
* Annotate RNA parts of contacts using reference annotation
* Perform voting to resolve conflicting annotations

12. Background Model Generation
-------------------------------
* Create a background model for normalization

13. Normalization
-----------------
* Normalize raw contacts using the background model
* Perform additional normalization steps (N2, scaling)

14. Peak Calling (for One-to-All experiments)
---------------------------------------------
* Use MACS2 to call significant peaks of chromatin-interacting RNAs

15. Statistics and Visualization
--------------------------------
* Generate statistics at various stages of the pipeline
* Create plots and visualizations of the results

16. MultiQC Report
------------------
* Compile a comprehensive quality control report using MultiQC

Main Results
------------
The main results of this pipeline are:

1. Pairs of RNA and DNA contacts, stored in tab-separated files
2. Annotation of the RNA parts of the contacts
3. Significant peaks of chromatin-interacting RNAs (for One-to-All experiments)
4. Various statistics and quality control metrics throughout the process

Note: This pipeline is flexible and can handle different types of all-to-all interactome sequencing data, with options to customize the workflow based on the specific experiment type and analysis requirements.