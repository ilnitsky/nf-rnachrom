Input file preparation
======================

Please follow through our guidelines for preparing gene annotation in a format compatible with various bioinformatics tools.

Gene Annotation
---------------

We require the gene annotation in ``GTF/GFF`` format suitable for the aligner and a ``"bed-like`` format with the following columns:

:: 

   chr     start     end     name     strand     biotype     source

Requirements for Bed-like Gene Annotation File
----------------------------------------------

- **1-based Indexing**: Use 1-based indexing for specifying coordinates.
- **Unique "name" Field**: Ensure that values in the "name" field are unique.
- **No Duplicate Coordinates**: Ensure that there are no genes with identical coordinates (absence of duplicates).
- **Chromosome Naming Consistency**: Chromosome names in the annotation file must match those in the genome assembly.

  .. note::
     For instance, when canonical chromosomes in the genome align with those in the GTF, but non-canonical ones do not, careful attention is required.

- **Single GTF File Usage**: If utilizing a single GTF file, you can convert it to a bed-like file using our script. 
If employing multiple annotations for a shared gene annotation, manual processing is required (refer to provided examples).

Shared Gene Annotation Issues
-----------------------------

When using multiple annotations, be aware of potential issues related to shared gene annotation: