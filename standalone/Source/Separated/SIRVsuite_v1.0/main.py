import argparse as ap

from Pipeline.Coverage.SIRVsuiteCoverage import SIRVsuiteCoverage
from Pipeline.Concentration.SIRVsuite_concentration import SIRVsuiteConcentration
from Pipeline.Correlation.ERCC_correlation import ERCCcorrelation
from Pipeline.helper import *


parser = ap.ArgumentParser()
parser.add_argument('-i','--sample-sheet', action = 'store', help = "Specify path to the sample sheet.", required = True, nargs = 1)
parser.add_argument('-o','--output-dir', action = 'store', help = "Specify output directory.", required = True, nargs = 1)
parser.add_argument('-a','--all-modules', action = 'store_true', help = "This options has the same functionality as --ERCC-correlation --SIRV-concentration --coverage and executes all 3 modules.", required=False)
parser.add_argument('--ERCC-correlation', action = 'store_true', help = "Specify to calculate correlation between expected and measured concentration for ERCC genes.")
parser.add_argument('--SIRV-concentration', action = 'store_true', help = "Specify to create heatmap and boxplots for SIRV E0 relative concentration.")
parser.add_argument('--coverage', action = 'store_true', help = "Specify to create coverages for all spike-in genes.")
parser.add_argument('--experiment-name',action = 'store', default="", help = "Specify name of an experiment.", required = False, nargs = 1)
args = parser.parse_args()

modules_to_execute = []

if (args.all_modules or args.ERCC_correlation or args.SIRV_concentration):
    modules_to_execute.append("concentration")

if (args.coverage or args.all_modules):
    modules_to_execute.append("coverage")

out = args.output_dir[0]

input_dict = read_sample_sheet("/home/tdrozd/development/sirv-suite/test/input/sample_sheet_merged.tsv", modules_to_execute = modules_to_execute)

if "concentration" in modules_to_execute:
    if (args.ERCC_correlation):
        a = ERCCcorrelation(sample_sheet=input_dict["concentration"], output_dir=out)
    
    if (args.SIRV_concentration):
        b = SIRVsuiteConcentration(sample_sheet=input_dict["concentration"], output_dir=out)
        b.create_sirvsuite_boxplot(c.data)
        b.create_sirvsuite_heatmap(c.data)

if "coverage" in modules_to_execute:
    c = SIRVsuiteCoverage(sample_sheet=input_dict["coverage"], output_dir=out, experiment_name = "")
    c.expected_coverage()
    c.bam_to_coverage()
    c.calc_statistics()
    c.coverage_plot()