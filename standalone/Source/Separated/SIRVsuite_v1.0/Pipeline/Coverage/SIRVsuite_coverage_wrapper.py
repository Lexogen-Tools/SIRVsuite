from SIRVsuite_coverage_tool import SIRVsuite_coverage as cpv
import numpy as np
import os
import argparse as ap
# "ENSG00000096384"

parser = ap.ArgumentParser()
parser.add_argument('-g','--gene-id', action = 'store', help = "Specify genes IDs for visualisation. If no gene is specified, the tool will process all genes available in annotation.", required = True, nargs = "+")
parser.add_argument('-o','--output-path', action = 'store', help = "Specify output folder path.", required = True, nargs = "+")
parser.add_argument('-f','--alignment-file', action = 'store', help = "Specify path to alignment BAM file.", required = False, nargs = "+")
parser.add_argument('-b','--bigwig-file', action = 'store', help = "Specify path to bigwig file.", required = False, nargs = "+")
parser.add_argument('-s','--strandeness', action = 'store', help = "Specify library strandeness in the same order as aligment-path arguments.", required = False, nargs = "+", choices = ["FWD","REV"])
parser.add_argument('-a','--annotation', action = 'store', help = "Specify path to the annotation file.", required = True, nargs = 1)
parser.add_argument('--UTRmode', action = 'store_true', help = "If specified, only coverage of an artifical UTR regions will be generated. Default value is 200 nt before 3' end and can be modified by --UTR-length option.")
parser.add_argument('--UTR-length', action = 'store', help = "Parameter which defines an artificial UTR length in number of nt. Default value is 200.", default = 200, required = False)
parser.add_argument('--transition-ends', action = 'store_true', help = "If specified, transition ends will be calculated for theoretical coverage.")
parser.add_argument('--sample-name', action = 'store', help = "Specify names of the samples, which will be visualized as a title in coverage plots.", default = "undefined", nargs = '+')
parser.add_argument('--experiment-name', action = 'store', default = "undefined", help = "Specify the name of an experiment, which will be visualized as a title in coverage plots.", nargs = 1)


args = parser.parse_args()
genes = args.gene_id

alignment_path = args.alignment_file
bigwig_path = args.bigwig_file

if (alignment_path == None and bigwig_path != None):
        input_mode = "bigwig"
elif (alignment_path != None and bigwig_path == None):
        input_mode = "alignment"
else:
        raise Exception("Error while loading input. Please choose either bigwig or BAM file input.")

strandeness = args.strandeness
annotation = args.annotation[0]
out_path = args.output_path[0]
transitionLengths = args.transition_ends
UTRmode = args.UTRmode
sample_name = args.sample_name
experiment_name = args.experiment_name[0]

if input_mode == 'alignment':
        if (len(strandeness) != len(alignment_path) ):
                raise Exception('different number of alignment files and library prep strandeness specified')
        else:
                bamFileList = np.transpose(np.array([alignment_path,strandeness]))

if (sample_name == 'undefined'):
        sample_name = []
        for path in alignment_path:
                file_name = os.path.basename(path)
                dir_names = os.path.split(os.path.dirname(path))
                sample_name.append(dir_names[-1] + "/" + file_name)
else:
        if (input_mode == 'alignment'):
                if (len(sample_name) != len(alignment_path) ):
                        raise Exception('different number of sample names and alignment files specified') 
                sample_name = sample_name
        elif (input_mode == 'bigwig'):
                if (len(sample_name) != len(bigwig_path) ):
                        raise Exception('different number of sample names and alignment files specified') 
                sample_name = sample_name

if (experiment_name == 'undefined'):
        experiment_name = ""

x = cpv(annotationFile = annotation,
        geneList = genes,
        sampleNames = sample_name,
        outputPath=out_path,
        UTRmode = UTRmode,
        UTRlength = args.UTR_length,
        transitionLengths = transitionLengths,
        experimentName = experiment_name,
        createBigWig = True)

x.process_annotation()
x.calculate_theoretical_coverage()

if input_mode == 'alignment':
        x.calculate_real_coverage(bamFileList = bamFileList)
elif input_mode == 'bigwig':
        x.calculate_real_coverage_BigWig(bigwig_path = bigwig_path)

x.calculate_stats()
x.visualize2png(1)
