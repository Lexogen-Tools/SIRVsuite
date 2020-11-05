# SIRVsuite

SIRVsuite is a command line tool to QC an RNA-Seq workflow using Lexogen's SIRV spike-in controls. For more specific details about Lexogen's SIRV spike-in controls visit: https://www.lexogen.com/sirvs/sirv-sets/.

SIRVsuite is published under the following [licence](https://github.com/Lexogen-Tools/SIRVsuite/blob/master/LICENCE).

**General usage**
```
SIRVsuite [-h] -i SAMPLE_SHEET -o OUTPUT_DIR [-a|--all-modules] [--coverage|--ERCC-correlation|--SIRV-concentration] [--experiment-name EXPERIMENT_NAME]
```

## Getting started

To get started with SIRVsuite analysis, you need to:

1. Install SIRVsuite
2. Prepare data + sample sheet
3. Run SIRVsuite

## 1. Installation

### Requirements

The following requirements have to be fulfilled prior to installation of SIRVsuite. When installing via docker these requirements are taken care of in Dockerfile.

Non-python requirements:
- cairo >= 1.15.10 (e.g., in ubuntu 16.04 and higher: libcairo2-dev)
- zlib
- libcurl (and the curl-config config) (e.g., in ubuntu: libcurl4-openssl-dev (this package combines both libcurl and curl-config))

Python requirements:
- numpy = 1.19.2 (recommended & tested)

### Installation method

#### a) PyPI

You can install SIRVsuite from GitHub by executing:

```
pip3 install git+https://github.com/Lexogen-Tools/SIRVsuite
```

It is recommended to create a python virtual environment or conda environment with python version 3.6-3.8 first. To create a python virtual environment, you can execute:

```
python3.6 -m venv sirvsuite
source sirvsuite/bin/activate
```

Once the sirvsuite venv has been created and activated, run the following commands to prepare the virtual env for installation:

```
pip3 install --upgrade pip
pip3 install numpy==1.19.2
```

Now, you can install SIRVsuite via:

```
pip3 install SIRVsuite
```

#### b) GitHub

First, all requirements mentioned above have to be installed or an environment needs to be created as in the previous example. Next, SIRVsuite can be installed directly from GitHub:

```
git clone https://github.com/Lexogen-Tools/SIRVsuite.git
python setup.py install
```

#### c) Docker

You can build a docker image from the cloned or downloaded SIRVsuite GitHub project. From inside the root directory of your local SIRVsuite project workspace execute the following command:

```
docker build . -t 'sirvsuite:latest'
```

Now, you can run SIRVsuite via the following command line:

```
docker run -v DATA_DIR_PATH:/data sirvsuite SIRVsuite [-h] -i SAMPLE_SHEET -o OUTPUT_DIR [-a|--all-modules] [--coverage|--ERCC-correlation|--SIRV-concentration] [--experiment-name EXPERIMENT_NAME]
```

Please note that you need to fill into the sample sheet alignment and counting paths which correspond to the mapped directory. 

Example usage:

From inside the root directory of your lcoal SIRVsuite project workspace you can invoke SIRVsuite on the included example data via:

```
docker run -it -v $(pwd):/data sirvsuite SIRVsuite -i /data/examples/sample_sheet_test_SIRVset4_docker.tsv -o /data/out -a
```

To make local files visible inside the docker container, the current working directory is mapped into the /data directory inside the docker container. Therefore, the -i and -o arguments and sample sheet path information need to be changed accordingly. The command will create out/ folder in the project root directory containing output data of all modules.

## 2. Preparing data + the sample sheet
The SIRVsuite consists of the following modules:

* SIRV coverage
* SIRV concentration
* ERCC correlation

The input for these modules is specified via a sample sheet, which is a csv file with the following format:

```
sample_name;alignment_path;counting_path;read_orientation;counting_feature;replication_group
sample_name_1;/home/user/alignment_data/sample_name_1.bam;/home/user/transcipt_count_data/sample_name_1.tsv;FWD;transcript;none
sample_name_2;/home/user/alignment_data/sample_name_2.bam;/home/user/transcipt_count_data/sample_name_2.tsv;FWD;transcript;rgroup_1
sample_name_3;/home/user/alignment_data/sample_name_3.bam;/home/user/transcipt_count_data/sample_name_3.tsv;FWD;transcript;rgroup_1
sample_name_4;/home/user/alignment_data/sample_name_4.bam;/home/user/transcipt_count_data/sample_name_4.tsv;FWD;transcript;none
```

### Alignment input

A sorted BAM file is required as input for the SIRVsuite coverage module. The SIRVsuite has been tested on BAM files computed by <a href=https://github.com/alexdobin/STAR>STAR aligner</a>. 

### Count input

A transcript quantification file is required as input for the SIRV concentration and ERCC correlation module. Since SIRVsuite was tested with <a href=https://www.lexogen.com/store/mix2-analysis-software>Mix<sup>2</sup></a> for transcript quantification, we assume that three columns of the  <a href=https://www.lexogen.com/store/mix2-analysis-software>Mix<sup>2</sup></a>  output format are present. These columns are given as follows:


| Field name | Example |                          Description                          |
| :--------- | :------ | :------------------------------------------------------------ |
| gene_ID | SIRV1 | A globally unique identifier for the genomic locus of the transcript. |
| tracking_ID | SIRV101 | A unique identifier for the transcript. |
| FPKM_CHN | 4158.42 | FPKM compatible hits norm. FPKM_CHN is calculated counting only the fragments,<br> which are compatible with a reference transcript. |

To use SIRVsuite in conjunction with another transcript quantification method than <a href=https://www.lexogen.com/store/mix2-analysis-software>Mix<sup>2</sup></a>, the relevant output of this transcript quantification has to be extracted into the format specified in the table above.

### Sample sheet content description

Sample sheets have the following format:
- ";" is the field separator
- trailing whitespace or tab is allowed
- UTF-8 encoding
- column names must be chosen from a predefined set (see info below)
- column values are case insensitive
- columns can be specified in arbitrary order

**Allowed columns**

Columns in the sample sheet can be divided into different categories. General columns are always required, the other columns relate to the module of interest.

General columns:
- sample_name: any set of characters to identify samples (this name will be printed in the final graphics).

SIRV-concentration & ERCC-correlation:

- counting_path: valid path to transcript quantification files.
- counting_feature: gene or transcript. Please note the transcript quantification files are required to run all the modules. Specifying "gene" allows to run only the coverage and ERCC correlation module. 
- replicate_group(optional): replicate groups definition, the same value assigned to multiple samples, their mean value will be used for quantification and will be used instead of sample names in the final graphics visualization. If a replicate group is to be defined for a subset of samples, use "none" to treat samples separately and use sample_names in the graphics instead.

Coverage:

- alignment_path: valid path to a .bam file.
- read_orientation: fwd or rev. It is necessary to specify library strandeness.

Any other column will be ignored.

## 3. Running SIRVsuite

SIRVsuite accepts the following arguments:

    required arguments:
      -i, --sample-sheet       path to the sample sheet
      -o, --output-dir         directory for output files
    
    selectively required arguments* (at least one is required):
      -a, --all-modules        triggers all available modules
      --ERCC-correlation       triggers processing of ERCCs ratios
      --SIRV-concentration     triggers processing of SIRVs relative concentration
      --coverage               triggers coverage processing
    
      * Note that using "-a" is the same as "--ERCC-correlation --SIRV-concentration --coverage".
       A valid usage is to specify at least one of the modules or to use -a argument.
    
    optional arguments:
      --experiment-name         name of the experiment to displayed in the final graphics (if empty, general title will be used)
      -h, --help                show help message and exit
      -v, --version             show version and exit

Example commands:
```
  # run SIRVsuite to generate coverage plots, data, relative abundance with heatmap and boxplot from sample information in sample_sheet_test_SIRVset4.tsv
  SIRVsuite -i CLONED_PROJECT_PATH/examples/sample_sheet_test_SIRVset4.tsv -o /home/user/SIRVsuite_output/ --experiment-name "sequencing-run-1" --coverage --SIRV-concentration

  # run SIRVsuite to perform whole analysis from sample information in sample_sheet.csv
  SIRVsuite -i CLONED_PROJECT_PATH/examples/sample_sheet_test_SIRVset4.tsv -o /home/user/SIRVsuite_output/ --experiment-name "sequencing-run-1" -a
```

## 4. Output data
The pipeline creates a subfolder for each specified module, which contains the output data of the module.

**Coverage module**

The module processes .bam files using the SIRV-set 3-4 annotation, calculates coverage (expected + measured) and creates 3 types of output:

- coverage/CoD_table_whole.tsv: Coefficient of Deviation (CoD) table,
- coverage/bigwig/SAMPLE_NAME_STRAND.bw coverage data in bigwig format,
- coverage/coverage_plots/SAMPLE_NAME/GENE_ID/coverage.png: coverage plot.

The CoD metric measures the similarity between expected and measured coverage. For the exon of a SIRV transcript the expected coverage is assumed to be a constant value. For overlapping exons from different SIRVs the expected coverage is the constant value times the number of overlapping exons. The CoD value is calculated from the expected coverage and the measured coverage multiplied by a scaling factor s as follows.

<p align="center"><img src="./docs/formulas/CoD.svg"></p>

Here, variable i ranges from the start to the end of the SIRV gene and s is chosen as follows.

<p align="center"><img src="./docs/formulas/scale_factor.svg"></p>

The CoD value is greater or equal to zero with values close to 0 indicating a good match between expected and measured coverage. 

The coverage plot generated by the coverage module visualizes the annotation of the SIRV gene in the upper part of the graphic and prints CoD values for the plus and minus strand. The lower part of the graphic shows the scaled measured coverage (grey) and expected coverage (blue and green for minus and plus strand). The x-axis of this graphic shows the coordinates of the SIRV exons, with introns plotted at a constant width. Where the scaled measured coverage has been truncated, this is visualized by a red arrow and a number corresponding to the maximal number of reads in the truncated region.

An example of a coverage plot is given below.

<p align="center"><img src="./docs/output_preview/coverage_preview.png" width=1280></p>

The spike-in coverage can also be inspected interactively, by loading the measured coverage (stored in <a href=http://genome.ucsc.edu/goldenPath/help/bigWig.html>bigwig</a> (.bw) format) in a genome browser (e.g. the <a href=http://software.broadinstitute.org/software/igv/UserGuide>IGV</a> browser).

**ERCC correlation module**

The module processes transcript or gene counts and input concentration of ERCCs. It creates two types of output:

- correlation/ERCC_correlation.tsv: correlation table
- correlation/SAMPLE_NAME/ERCC_correlation.png: correlation plot

The correlation table consist of Pearson correlation R values between measured and true ERCC concentrations. The true ERCC concentrations are obtained from <a href=https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt>this table</a>. The correlation plot, generated by the ERCC correlation module, is a scatterplot showing the true concentration of the ERCCs on the x-axis and the measured concentration on the y-axis. An example of such a correlation plot is given below.

<p align="center"><img src="./docs/output_preview/ERCC_correlation_preview.png" width=800></p>

**SIRV concentration module**

SIRV transcripts in the SIRV E0 mix (present in SIRV set 3 and 4) are equimolar. Hence, assuming the sum of the FPKM values of all SIRV transcripts correctly quantifies the overall concentration of the SIRVs, the expected FPKM value of each SIRV is given by:

<p align="center"><img src="./docs/formulas/FPKM_exp.svg"></p>

The observed FPKM value relative to the expected FPKM value can then be calculated as follows:

<p align="center"><img src="./docs/formulas/FPKM_rel.svg"></p>

The SIRV concentration module calculates these relative FPKM values and generates the following 3 output files in the output directory. The latter two files visualize the data in the first:

- concentration/relative_concentration.tsv: table of relative SIRV concentrations
- concentration/SIRV_boxplot.png: boxplot of relative SIRV concentrations
- concentration/SIRV_heatmap.png: heatmap of Log<sub>2</sub> relative SIRV concentrations

The FPKM<sub>i,rel</sub> values can be found in concentration/relative_concentration.tsv. The relative FPKM values of all SIRV transcripts are visualized as a boxplot in concentration/SIRV_boxplot.png. An example of such a boxplot is given below. The mean and median are visualized by the red circle and white strip, respectively.

<p align="center"><img src="./docs/output_preview/SIRV_boxplot_preview.png" width=800></p>

The output file concentration/SIRV_heatmap.png displays a heatmap of the log<sub>2</sub> relative concentrations for each SIRV transcript grouped by SIRV gene. The values are visualized on a red/green/blue color map, with green indicating log<sub>2</sub> values close to 0, red and blue indicating values close to 1 and -1, respectively. Values above 1 and below -1 are visualized by red and green. An example of a SIRV concentration heatmap for SIRV1 and SIRV2 is given below.

<p align="center"><img src="./docs/output_preview/SIRV_heatmap_preview_SIRV12.png" width=350></p>
