.. _setup_and_installation:

Setup and Installation
======================


This guide will walk you through the process of setting up and installing the RNAchrom Nextflow pipeline.



Prerequisites
-------------
Before you begin, ensure you have the following installed on your system::

    Java 8 or later
    Nextflow (version 20.04.0 or later)
    Anaconda

The pipeline uses Conda environments specified in ``envs/`` directory . Ensure you have Conda installed on your system.

Installing Nextflow
-------------------
    To install Nextflow, run the following commands:

.. code-block:: bash

    $ curl -s https://get.nextflow.io | bash
    $ mv nextflow ~/bin/

Make sure ~/bin is in your PATH.

Cloning the Repository
Clone the RNAchrom repository:

.. code-block:: bash

    $ git clone https://github.com/your-repo/rnachrom.git
    $ cd rnachrom

Pipeline Dependencies
---------------------
The pipeline automatically installs several tools during execution. These include::

    RnaChromATA
    Bitap
    StereoGene
    fastq-dupaway

Installation Process
--------------------


The PrepareSoftware process in the pipeline handles the installation of these tools:

**RnaChromATA**


Installed from the bin/RnaChromATA/setup.py file in the project directory.

**Bitap**


Cloned from the GitHub repository and compiled:

.. code-block:: bash

    $ git clone https://github.com/ikm4rkov/RawReadsProcessor.git


**StereoGene**


Cloned from GitHub and compiled:

.. code-block:: bash

    $ git clone https://github.com/favorov/stereogene.git
    $ cd stereogene/src
    $ make

**fastq-dupaway**


Cloned from GitHub and compiled:

.. code-block:: bash

    $ git clone https://github.com/AndrewSigorskih/fastq-dupaway.git
    $ cd fastq-dupaway
    $ make


Running the Pipeline
--------------------


Before running pipeline with real data, try running the pipeline in test mode. Thus,
all required tools and conda environments will be installed. 

To run the pipeline:

.. code-block:: bash

    $ nextflow run main.nf --input samples.csv --outdir results

Troubleshooting
---------------
If you encounter issues during installation or execution:

Ensure all prerequisites are correctly installed.
Check the log files in the work directory for error messages.
Verify that all required input files are present and correctly formatted.
For more detailed information, refer to the pipeline documentation or open an issue on the GitHub repository.