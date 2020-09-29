SIRVsuite is a command line tool to analyze performance of SIRV set 3
and 4 spike-ins based on alignment and transcript count data.

SIRVsuite is permitted under the following licence xxx.

**General usage**

::

    python SIRVsuite.py [-h] -i SAMPLE_SHEET -o OUTPUT_DIR [-a|--all-modules] [--coverage|--ERCC-correlation|--SIRV-concentration] [--experiment-name EXPERIMENT_NAME]

Getting started
===============

1. Installation
===============

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
requirements in terms of input files and necessary parameters for ther
processing. Therefore, a .csv file comprised of such information needs
to be created. We call this type of file a sample sheet. An example of a
valid sample sheet:

::

    sample_name;alignment_path;counting_path;read_orientation;counting_method;counting_feature;library_prep_type;replication_group
    sample_name_1;/home/user/alignment_data/sample_name_1.bam;/home/user/transcipt_count_data/sample_name_1.tsv;FWD;mix2;transcript;whole
    sample_name_2;/home/user/alignment_data/sample_name_2.bam;/home/user/transcipt_count_data/sample_name_2.tsv;FWD;mix2;transcript;whole

The SIRVsuite tool will automatically check whether specified module can
be processed based on the sample sheets.

Sample sheet content description
--------------------------------

Sample sheet is required to have a following format: - used ";"
separator - trailing whitespace or tab is allowed (trimmed during
reading process) - UTF-8 encoding - every column must have a name which
is predefined (see more info below) - column values are case insensitive
- the order of columns can be arbitrary

**Allowed columns**:

General columns: 1. sample\_name: any set of characters to identify
samples (this will be printed in the final graphics).

2. library\_prep\_type: whole (whole transcriptome library prep) or qs
   (QuantSeq library preparation)

Concentration specific columns:

3. counting\_path: valid path to count files.
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
====================

.. math::
   FPKM_{expected} = \frac{\sum_{i=1}^{n}FPKM_i}{n}

   FPKM_{i,rel} = \frac {FPKM_i} {FPKM_{expected}} = \frac {FPKM_i} {\sum_{i=1}^{n}FPKM_i} / \frac {1} {n}
