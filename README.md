# ![nf-core/rnachrom](docs/images/nf-core-rnachrom_logo_light.png#gh-light-mode-only) ![nf-core/rnachrom](docs/images/nf-core-rnachrom_logo_dark.png#gh-dark-mode-only)


[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)

## Introduction


Full documentation is available at ([`ReadTheDocs`](https://nf-rnachrom.readthedocs.io/en/latest/))


**nf-core/rnachrom** is a comprehensive and flexible bioinformatics pipeline designed to process RNA-DNA interactome sequencing data. It efficiently handles large-scale data from various experimental methods including all-to-all approaches (GRID-seq, RADICL-seq, iMARGI, Red-C) and one-to-all approaches (ChART, RAP, CHIRP). The pipeline provides a streamlined workflow for analyzing RNA-DNA interactions, from raw sequencing data to annotated contacts and statistical analyses.

The pipeline supports:
- Both all-to-all (ATA) and one-to-all (OTA) RNA-chromatin interaction protocols
- Bridge sequence processing for ATA methods like Red-C, ChAR-seq, GRID-seq and RADICL-seq
- Multiple alignment strategies for both RNA and DNA components
- Comprehensive annotation of interaction sites with genomic features
- Statistical analysis and normalization of interaction data
- Integration with matching RNA-seq data for enhanced analysis (termed "chromatin potential")

## Pipeline Workflow

The pipeline implements a multi-stage workflow that handles different types of RNA-chromatin interaction data:

1. **Data Preprocessing**
   - Quality control with FastQC
   - Adapter trimming with FastP
   - Optional deduplication for PCR duplicates

2. **Protocol-Specific Processing**
   - **One-to-all methods** (ChART, RAP, CHIRP): Direct trimming and alignment
   - **All-to-all methods**:
     - **Separate RNA/DNA  reads**: Sorting and processing RNA and DNA parts separately
     - **Chimera reads with linker**: Bridge sequence identification and processing
     - **iMARGI**: BWA alignment and specialized processing

3. **Read Processing and Alignment**
   - Short read trimming and filtering
   - PEAR merging for paired-end reads
   - Bridge splitting with custom tools
   - Alignment using appropriate tools (STAR, BWA, Bowtie2, HISAT2)

4. **Contact Extraction and Processing**
   - BAM to contacts conversion
   - Edit distance filtering and CIGAR filtering
   - Integration of RNA-DNA parts into contact tables
   - Contact normalization and deduplication

5. **Annotation and Analysis**
   - RNA-DNA contact annotation with genomic features
   - RNA part annotation with transcriptome data
   - Strand detection for RNA components
   - Optional chromatin potential calculation with RNA-seq data

6. **Downstream Analysis**
   - Peak calling with MACS2 or BARDIC
   - Statistical analysis and visualization
   - Final contact table generation with header information
   - Comprehensive reporting and data visualization

The workflow adaptively handles different experimental protocols and data formats, providing a complete solution from raw sequencing data to biologically interpretable results.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->

<!-- ## Default Steps

The pipeline includes the following major steps:

1. Read QC and preprocessing ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`FastP`](https://github.com/OpenGene/fastp))
2. Bridge detection and read splitting (protocol-specific)
3. Alignment of RNA and DNA components ([`STAR`](https://github.com/alexdobin/STAR), [`BWA`](https://github.com/lh3/bwa), [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/), [`HISAT2`](http://daehwankimlab.github.io/hisat2/))
4. Contact extraction and filtering
5. Annotation and genomic feature assignment
6. Statistical analysis and visualization
7. Comprehensive report generation ([`MultiQC`](http://multiqc.info/)) -->

## Usage

:information_source: **Note:**
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.

Nextflow Windows installation: [this page](https://www.nextflow.io/blog/2021/setup-nextflow-on-windows.html)

### Running the pipeline

First, prepare a samplesheet with your input data according to your experiment type. The pipeline supports different types of RNA-chromatin interaction experiments:

#### All-to-all (ATA) methods (e.g., MARGI, RADICL-seq)

A simple ATA samplesheet should include the RNA and DNA parts of the experiment:

```csv
sample,rna,dna
HFFc6_imargi_1,SRR8206679_1.fastq.gz,SRR8206679_2.fastq.gz
HFFc6_imargi_2,SRR8206680_1.fastq.gz,SRR8206680_2.fastq.gz
```

#### ATA with RNA-seq control

For ATA experiments with matching RNA-seq data, use the following format:

```csv
sample,rna,dna,description
HFFc6_imargi_1,SRR8206679_1.fastq.gz,SRR8206679_2.fastq.gz,tissue:"HFFc6";rnaseq:"rnaseq_HFFc6"
HFFc6_imargi_2,SRR8206680_1.fastq.gz,SRR8206680_2.fastq.gz,tissue:"HFFc6";rnaseq:"rnaseq_HFFc6"
rnaseq_HFFc6,SRR8206681_1.fastq.gz,SRR8206681_2.fastq.gz,tissue:"HFFc6"
```

Note that RNA-seq samples should always start with `rnaseq_` prefix.

#### One-to-all (OTA) methods (e.g., ChART, RAP, CHIRP)

For one-to-all methods, use the format:

```csv
sample,fastq_1,fastq_2,description
SCIRT_chart_39,SRR10044362.fastq.gz,,tissue:"K562"
```

If you have single-end data, leave fastq_2 empty. For paired-end data, include both files.

#### OTA with input controls

To specify input controls for OTA methods:

```csv
sample,fastq_1,fastq_2,control,description
SCIRT_chart_39,SRR10044362.fastq.gz,,SCIRT_INPUT,tissue:"K562"
SCIRT_INPUT,SRR10044359.fastq.gz,,
```

Use the `control` column to specify which sample serves as the input control.

### Running the pipeline


Check if you have Apptainer, Conda, or Docker installed:
- Apptainer (formerly Singularity): [Installation Guide](https://apptainer.org/docs/admin/main/installation.html)
- Conda: [Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- Docker: [Installation Guide](https://docs.docker.com/get-docker/)



Clone the repository with Nextflow project:
```bash
git clone https://github.com/ilnitsky/nf-rnachrom.git
```


Now, you can run the pipeline using:

```bash
nextflow run ./nf-rnachrom \
   -profile <docker/apptainer/conda/...> \
   --input samplesheet.csv \
   --exp_type <redc/imargi/radicl/chart/rap/chirp/...> \
   --outdir <OUTDIR>
```

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::

## Nextflow Modules implemented in this pipeline


This pipeline includes modules from both local and nf-core sources:

- **Local modules** (`./modules/local`): Custom implementations specific to RNA-chromatin interaction analysis
- **nf-core modules** (`./modules/nf-core`): Standard bioinformatics modules from the nf-core library


### Local Modules



- **Quality Control**
  - `fastqc`: Quality control of FASTQ files
  - `multiqc`: Aggregates bioinformatics analyses into a single report
  - `fastp`: All-in-one FASTQ preprocessor

- **Deduplication**
  - `fastq-dupaway`: 
  - `fastuniq`: 
  - `clumpify`: Removes duplicate reads from FASTQ files
  - `seqkit_rmdup`:


- **Trimming**
  - `trimgalore`: Trim adapter sequences and low quality regions
  - `trimmomatic`: Flexible read trimming tool
  - `bbduk`: BBMap suite

- **Managing bridge sequence**
  - `bitap`: 
  - `chartools`: 
  - `tagdust`: 

- **Preprocessing**
  - `dedup`: Removes duplicates from sequencing data
  - `debridge`: Processes bridge sequences in MARGI/RADICL-seq data
  - `rsites`: Analyzes restriction enzyme sites
  - `bam_to_contacts`: Converts aligned BAM files to RNA-DNA contact files
  - `filter_contacts`: Filters RNA-DNA contact files based on various criteria

- **Alignment**
  - `star`: RNA-seq read alignment
  - `hisat2`: Hierarchical indexing for spliced alignment of transcripts
  - `bwa`: Burrows-Wheeler Aligner for DNA sequences
  - `bowtie2`: Ultrafast short read aligner

- **Analysis**
  - `annotation`: Annotates RNA-DNA contacts with genomic features
  - `normalisation`: Performs normalization of contact data
  - `chromatin_potential`: Computes chromatin interaction potential
  - `bardic`: Implements BARDIC algorithm for interaction analysis
  - `xrna_assembly`: Assembles novel RNA transcripts from RNA parts of contacts
  - `detect_strand`: Determines strand orientation for RNA-DNA contacts


- **Peak Calling**
  - `macs2`: Model-based Analysis of ChIP-Seq data
  - `bardic`: Implements BARDIC algorithm for peak calling 

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/rnachrom/usage) and the [parameter documentation](https://nf-co.re/rnachrom/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/rnachrom/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/rnachrom/output).

The pipeline produces the following outputs:


###  Analysis Results
- Raw RNA-DNA contact tables in TSV format
- Annotated and normalized interaction data (contacts) with genomic features
- Peak calls
- Chromatin potential calculations (if RNA-seq provided)

- MultiQC reports aggregating pipeline statistics
- Statistical summaries and visualization data
- Log files for troubleshooting



## Test Run with External Data

Download test dataset:

```bash
wget http://bioinf.fbb.msu.ru/ken/nextflow/test_data_nf-rnachrom.tar.gz && \
tar -xzvf test_data_nf-rnachrom.tar.gz -C test_data_nf-rnachrom 
```

To test the pipeline functionality using pre-configured test data, you can run:

```bash
nextflow run ./nf-rnachrom \
  -profile test_full,apptainer \
  --outdir test_results
```

This test run will download example GRID-seq data from bioinf.fbb.msu.ru/ken and 
process the data through all pipeline stages. It will prepare GRCh38 index for genome.

The `test_full` profile automatically configures all necessary parameters including reference genomes, annotation files, and processing options optimized for this test dataset.

> **Note:** The full test requires approximately 8GB of RAM and 4 CPU cores. At least 10 GB of free disk space is required: the download size is around 5000MB, the apptainer image size is 4000MB. The test should complete in about 30-45 minutes on a standard workstation.

## Credits

nf-core/rnachrom was originally written by Ivan Ilnitskiy.

We thank the following people for their extensive assistance in the development of this pipeline:
**Ivan Markov** (Moscow State University), **Arina Nikolskaya** (Moscow State University), **Anastasia Zharikova**

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#rnachrom` channel](https://nfcore.slack.com/channels/rnachrom) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/rnachrom for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
