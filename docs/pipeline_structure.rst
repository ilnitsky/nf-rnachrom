Pipeline Structure of  nf-rnachrom 
=============================================================



The workflow is modularized for flexibility. Here's an overview of the main components:

**Modules**

- **Input Check**: Validates and processes the input samplesheet.
- **FastQC**: Performs quality control checks on raw sequence data.
- **Deduplication (DEDUP)**: Removes duplicate reads to ensure data accuracy.
- **Trimming (TRIM)**: Trims adaptors and low-quality bases from reads.
- **Alignment**: Aligns reads to the reference genome using tools like HISAT2.
- **Post-Alignment Processing**: Includes deduplication, BAM sorting, and conversion to BED format.
- **Experiment-Specific Workflows**:
  - **OTA/ATA**: Processes one-to-all or all-to-all experiments.
  - **Peak Calling**: Uses MACS2 for peak detection in RAP, CHIRP, and similar experiments.
  - **Contact Mapping**: Analyzes interactions in GRID-seq and related methods.

**Workflow Diagram**

.. code-block:: text

   Input Check --> FastQC --> [Dedup] --> [Trim] --> Alignment --> Post-Alignment --> Analysis

