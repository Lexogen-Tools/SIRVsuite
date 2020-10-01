=========
SIRVsuite
=========

SIRVsuite is a command line tool to analyze performance of SIRV set 3
and 4 spike-ins based on alignment and transcript count data.

SIRVsuite is permitted under the following licence xxx.

**General usage**

::

    python SIRVsuite.py [-h] -i |  SAMPLE_SHEET -o OUTPUT_DIR [ [-a|--all-modules] | --coverage | --ERCC-correlation | --SIRV-concentration ] [--experiment-name EXPERIMENT_NAME]

Getting started
---------------

1. Installation
---------------

To install SIRVsuite, an environment for all depedent packages needs to
be created. To this end, install/sirvsuite\_env.yml can be used via
conda to create one via

::

    conda env create -p PATH_TO_CONDA -f install/sirvsuite_env.yml

Another option would be to install python packages directly, which is
not recommended due to possible dependency conflicts.

2. Preparing sample sheet
=========================

The SIRVsuite consists of the following modules: coverage, SIRV
concentration and ERCC correlation. The modules have different
requirements in terms of input files type and necessary information for processing. Therefore, a .csv file comprised of such information needs
to be created. We call this type of file a sample sheet. An example of a
valid sample sheet:

::

    sample_name;alignment_path;counting_path;read_orientation;counting_method;counting_feature;library_prep_type;replication_group
    sample_name_1;/home/user/alignment_data/sample_name_1.bam;/home/user/transcipt_count_data/sample_name_1.tsv;FWD;mix2;transcript;whole;1
    sample_name_2;/home/user/alignment_data/sample_name_2.bam;/home/user/transcipt_count_data/sample_name_2.tsv;FWD;mix2;transcript;whole;2

The SIRVsuite tool will automatically check whether specified module can
be processed based on the sample sheets.

Sample sheet content description
--------------------------------

Sample sheet is required to have a following format:

- used ";" separator
- trailing whitespace or tab is allowed (trimmed during reading process)
- UTF-8 encoding
- every column must have a name which is predefined (see more info below)
- column values are case insensitive
- the order of columns can be arbitrary

**Allowed columns**:

General columns: 1. sample\_name: any set of characters to identify
samples (this will be printed in the final graphics).

2. library\_prep\_type: whole (whole transcriptome library prep) or qs
   (QuantSeq library preparation)

Concentration specific columns:

3. counting\_path: valid path to the count file.
4. counting\_method: mix2, cufflinks or htseq. Defines, how the count
   file should be read.
5. counting\_feature: gene or transcript. Defines, whether ERCC
   correlation plots and tables are quantified (gene counts) or ERCC
   correlation + SIRV heatmap and boxplot are quantified (transcript
   counts).
6. replicate\_group(optional): any set of characters to identify
   replicate groups (this will be printed in the final graphics). An
   optional column, if the same replicate group assigned to multiple
   samples, their mean value will be used for quantification

Coverage specific columns:

7. alignment\_path: valid paths to a .bam file.
8. read\_orientation: fwd, rev or none. Use "fwd" or "rev" for
   strand-specific libraries, "none" for non-strand specific libraries.

Any other columns will be ignored.

3. Running SIRVsuite
--------------------

SIRVsuite accepts the following arguments:
::

    required arguments:
      -i, --sample-sheet       path to the sample sheet
      -o, --output-dir         directory for output files

    selectively required arguments* (at least one is required):
      -a, --all-modules        triggers all available modules
      --ERCC-correlation        triggers processing of ERCCs ratios
      --SIRV-concentration      triggers processing of SIRVs relative concentration
      --coverage                triggers coverage processing

      *Note that using "-a" is the same as "--ERCC-correlation --SIRV-concentration --coverage".
       A valid usage is to specify at least one of the modules or to use -a argument.

    optional arguments:
      --experiment-name         name of the experiment to displayed in the final graphics (if empty, general title will be used)
      -h, --help                show help message and exit

Example commands:
::

  # run SIRVsuite to generate coverage plots, data, relative abundance with heatmap and boxplot from sample information in sample_sheet.csv
  python SIRVsuite.py -i sample_sheet.csv -o /home/user/SIRVsuite_output/ --experiment-name "sequencing-run-1" --coverage --SIRV-concentration

  # run SIRVsuite to perform whole analysis from sample information in sample_sheet.csv
  python SIRVsuite.py -i sample_sheet.csv -o /home/user/SIRVsuite_output/ --experiment-name "sequencing-run-1" -a


Output data:
------------

The pipeline will create subfolder for every specified module with module-specific output data.

**Coverage module**

The module processes .bam files + SIRV-set 3-4 annotation, calculates coverage (expected + measured) and creates 3 types of output:

- CoD table,
- coverage data in bigwig format,
- coverage plot.

The CoD metrics allows to measure the resemblence between expected (theoretical) and measured (real) coverage. The theoretical coverage is calculated based on annotated distribution of exons, whilst the measured coverage is quantified from the reads obtained from the sequencer. For CoD applies
CoD >= 0

Measured coverage in bigwig (.bw) format can be used, for example, in a IGV browser to inspect spike-in coverage interactively. See more info about bigwig: http://genome.ucsc.edu/goldenPath/help/bigWig.html.

**ERCC-Correlation module**

The module processes transcript or gene counts and input concentration of ERCCs. It creates two types of output:

- correlation table,
- correlation plot.

The correlation table consist of R^2 values for each sample, the correlation plot displays an overview of distribution of ERCC gene concentration ratios.

**SIRV-concentration module**

The module processes transcript FPKM values for SIRVs and creates 3 types of output:

- boxplot,
- relative concentration table,
- heatmap.
