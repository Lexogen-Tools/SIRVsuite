# SIRVsuite

SIRVsuite is a command line tool to analyze performance of SIRV set 3 and 4 spike-ins based on alignment and transcript count data. For more specific details about SIRV sets visit: https://www.lexogen.com/sirvs/sirv-sets/.

SIRVsuite is permitted under the following licence xxx.

**General usage**
```
SIRVsuite [-h] -i SAMPLE_SHEET -o OUTPUT_DIR [-a|--all-modules] [--coverage|--ERCC-correlation|--SIRV-concentration] [--experiment-name EXPERIMENT_NAME]
```

## Getting started

## 1. Installation

### Requirements
Non-python requirements:
- cairo >= 1.15.10
- zlib
- libcurl (and the curl-config config)

### a) Docker

<!--
To install SIRVsuite, an environment for all depedent packages needs to be created. Thus, install/sirvsuite_env.yml can be used via conda command
```
conda env create -p CONDA_PATH/envs/sirvsuite -f install/sirvsuite_env.yml
```

to create a virtual conda environment, from which SIRVsuite.py can run. Conda is installed to the home directory by default. In this case CONDA_PATH would refer to /home/user_name/anaconda3 or /home/user_name/miniconda3. See more info about conda enironments: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html. 

Another option is to install python packages directly, which is not recommended due to possible dependency conflicts.
-->

You can build a docker image by running from a project directory:

```
docker build . -t 'sirvsuite:latest'
```

Then, you can run the SIRVsuite via command

```
docker run -v DATA_DIR_PATH:/data SIRVsuite [-h] -i SAMPLE_SHEET -o OUTPUT_DIR [-a|--all-modules] [--coverage|--ERCC-correlation|--SIRV-concentration] [--experiment-name EXPERIMENT_NAME]
```

Please note that yous need to fill into the sample sheet alignment and counting paths which correspond to the mapped directory. 

Example usage:

First, locate to the cloned repository root location. You can then invoke SIRVsuite via:

```
docker run -it -v $(pwd):/data sirvsuite SIRVsuite -i /data/examples/sample_sheet_test_SIRVset4_docker.tsv -o /data/out -a
```

To make local files visible inside of docker container, the current working directory is mapped into /data directory inside of docker container. Therefore, -i and -o arguments and sample sheet path information need to be changed accordingly. The command will create out/ folder in the project root directory containing output data of all modules.

### b) PyPI

Currently being implemented...

Soon available via: 
```
pip install SIRVsuite
```

### c) github

Valid after publishing on github

### d) gitlab

You can install sirv-suite using gitlab repo for internal testing purposes.

It is recommended to create a virtual python environment or conda environment with python version 3.6-3.8 first. Next, 
after env activation you can run:

```
pip3 install git+http://my_token:sZtLBXmrwFFzmvLiyp-c@10.90.1.56:10080/Bioinfo/sirv-suite.git@SIRVsuite_v0.1
```

## 2. Preparing sample sheet
The SIRVsuite consists of the following modules: coverage, SIRV concentration and ERCC correlation. The modules have different requirements in terms of input files and necessary parameters for ther processing. Therefore, a .csv file comprised of such information needs to be created. We call this type of file a sample sheet. An example of a valid sample sheet:

```
sample_name;alignment_path;counting_path;read_orientation;counting_method;counting_feature;library_prep_type;replication_group
sample_name_1;/home/user/alignment_data/sample_name_1.bam;/home/user/transcipt_count_data/sample_name_1.tsv;FWD;mix2;transcript;whole
sample_name_2;/home/user/alignment_data/sample_name_2.bam;/home/user/transcipt_count_data/sample_name_2.tsv;FWD;mix2;transcript;whole
```

The SIRVsuite tool will automatically check whether specified module can be invoked based on the sample sheet specification.

### Sample sheet content description

Sample sheet is required to have a following format:
- used ";" separator
- trailing whitespace or tab is allowed (trimmed during reading process)
- UTF-8 encoding
- every column must have a name which is predefined (see more info below)
- column values are case insensitive
- the order of columns can be arbitrary

**Allowed columns**:

Columns in the sample sheet can be divided into different categories. General columns are always required, the other columns relate to the module of interest.

General columns:
- sample_name: any set of characters to identify samples (this will be printed in the final graphics).
- library_prep_type: whole (whole transcriptome library prep) or qs (QuantSeq library preparation)

SIRV-concentration & ERCC-correlation:

- counting_path: valid path to count files.
- counting_method: mix2, cufflinks or htseq. Defines, how the count file should be read.
- counting_feature: gene or transcript. Defines, whether ERCC correlation plots and tables are quantified (gene counts) or ERCC correlation + SIRV heatmap and boxplot are quantified (transcript counts).
- replicate_group(optional): replicate groups definition, the same value assigned to multiple samples, their mean value will be used for quantification and will be used instead of sample names if the final graphics visualization. If replicate group is to be defined for a subset of samples, use "none" value to treat samples separately and use sample_names in the graphics instead.

Coverage:

- alignment_path: valid path to a .bam file.
- read_orientation: fwd, rev or none. Use "fwd" or "rev" for strand-specific libraries, "none" for non-strand specific libraries.

Any other column will be ignored.

## 3. Running SIRVsuite

SIRVsuite accepts the following arguments:

    required arguments:
      -i, --sample-sheet       path to the sample sheet
      -o, --output-dir         directory for output files

    selectively required arguments* (at least one is required):
      -a, --all-modules        triggers all available modules
      --ERCC-correlation        triggers processing of ERCCs ratios
      --SIRV-concentration      triggers processing of SIRVs relative concentration
      --coverage                triggers coverage processing

      * Note that using "-a" is the same as "--ERCC-correlation --SIRV-concentration --coverage".
       A valid usage is to specify at least one of the modules or to use -a argument.

    optional arguments:
      --experiment-name         name of the experiment to displayed in the final graphics (if empty, general title will be used)
      -h, --help                show help message and exit

Example commands:
```
  # run SIRVsuite to generate coverage plots, data, relative abundance with heatmap and boxplot from sample information in sample_sheet.csv
  python SIRVsuite.py -i sample_sheet.csv -o /home/user/SIRVsuite_output/ --experiment-name "sequencing-run-1" --coverage --SIRV-concentration

  # run SIRVsuite to perform whole analysis from sample information in sample_sheet.csv
  python SIRVsuite.py -i sample_sheet.csv -o /home/user/SIRVsuite_output/ --experiment-name "sequencing-run-1" -a
```

## 4. Output data
The pipeline will create subfolder for every specified module with module-specific output data.

**Coverage module**

The module processes .bam files + SIRV-set 3-4 annotation, calculates coverage (expected + measured) and creates 3 types of output:

- CoD table,
- coverage data in bigwig format,
- coverage plot.

The CoD metrics allows to measure the resemblence between expected (theoretical) and measured (real) coverage. The theoretical coverage is calculated based on annotated distribution of exons, whilst the measured coverage is quantified from the reads obtained from the sequencer.

CoD value is given as follows:

<p align="center"><img src="https://latex.codecogs.com/svg.latex?\Large&space;CoD=\frac{\sum_{i=1}^{n}(cov_{theory}-cov_{real,scaled})^2}{\sum_{i=1}^{n}cov_{theory,i}}"\></p>

 For CoD applies:
 CoD >= 0
 , where values around 0 indicate an ideal match between expected and measured coverage.

Measured coverage in bigwig (.bw) format can be used, for example, in a IGV browser to inspect spike-in coverage interactively. See more info about bigwig: http://genome.ucsc.edu/goldenPath/help/bigWig.html and IGV: http://software.broadinstitute.org/software/igv/UserGuide.

Coverage plot can serve as an overview of exon distribution for different transcript variants, the corresponding expected coverage based on annotation and the measured read distribution fits into these regions. In addition, it provides a basic statistics along with a CoD value.

An example of a coverage plot:

<p align="center"><img src="./docs/output_preview/coverage_preview.png" width=1280></p>

**ERCC-Correlation module**

The module processes transcript or gene counts and input concentration of ERCCs. It creates two types of output:

- correlation table,
- correlation plot.

The correlation table consist of Pearson correlation R values, the correlation plot is a scatterplot revealing dependency between theoretical and measured concentration for ERCC genes based on concentration table: https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt.

An example of a correlation plot:

<p align="center"><img src="./docs/output_preview/ERCC_correlation_preview.png" width=800></p>

**SIRV-concentration module**

The module processes transcript FPKM values for SIRVs and creates 3 types of output:

- relative concentration table (concentration/relative_concentration.tsv),
- boxplot of SIRV relative transcript concentration (concentration/SIRV_boxplot.png),
- heatmap of Log<sub>2</sub> Fold Change SIRV transcript relative concentrations (concentration/SIRV_heatmap.png).

The SIRV E0 mix (present in SIRV set 3 and 4) is comprised of equimolar transcripts. This enables the folowing calculation procedure. Consider SIRV transcript i of n SIRV transcripts present in a sample. Given a measured FPKM concentration FPKM<sub>i</sub>,  FPKM<sub>expected</sub> for each transcript is quantified as follows:

<p align="center"><img src="https://latex.codecogs.com/svg.latex?\Large&space;FPKM_{i,expected}=\frac{\sum_{i=1}^{n}FPKM_i}{n}"\></p>

and we can define relative FPKM value as a ratio of estimated and expected relative abundance using formula

<p align="center"><img src="https://latex.codecogs.com/svg.latex?\Large&space;FPKM_{i,rel}=\frac{FPKM_i}{FPKM_{expected}}=\frac{FPKM_i}{\sum_{i=1}^{n}FPKM_i} \ \frac{1}{n}"\></p>

The FPKM<sub>i,rel</sub> values can be found in relative concentration/relative_concentration.tsv in the output directory.

Log<sub>2</sub> Fold Change of FPKM<sub>i,rel</sub> is displayed in a heatmap, showing the difference between expected and calculated values. The green signalizes match between concentrations, blue transcript underexpression and red overexpression.

The distribution of FPKM<sub>i,rel</sub> is summarized into a boxplot.

An example of a SIRV concentration boxplot (mean - read circle, median - white strip):

<p align="center"><img src="./docs/output_preview/SIRV_boxplot_preview.png" width=800></p>

An example of a SIRV concentration heatmap for SIRV1 and SIRV2:
<p align="center"><img src="./docs/output_preview/SIRV_heatmap_preview_SIRV12.png" width=350></p>
