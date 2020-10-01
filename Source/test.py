"""
import argparse as ap

from Pipeline.Coverage.SIRVsuiteCoverage import SIRVsuiteCoverage
from Pipeline.Concentration.SIRVsuite_concentration import SIRVsuiteConcentration
from Pipeline.Correlation.ERCC_correlation import ERCCcorrelation
from Pipeline.helper import *


parser = ap.ArgumentParser()
parser.add_argument('-s','--sample-sheet', action = 'store', help = "Specify path to the sample sheet.", required = True, nargs = 1)
parser.add_argument('-o','--output-dir', action = 'store', help = "Specify output directory.", required = True, nargs = 1)
parser.add_argument('-a','--all-modules', action = 'store_true', help = "This options has the same functionality as --ERCC-correlation --SIRV-concentration --coverage and executes all 3 modules.", required=False)
parser.add_argument('--ERCC-correlation', action = 'store_true', help = "Specify to calculate correlation between expected and measured concentration for ERCC genes.")
parser.add_argument('--SIRV-concentration', action = 'store_true', help = "Specify to create heatmap and boxplots for SIRV E0 relative concentration.")
parser.add_argument('--coverage', action = 'store_true', help = "Specify to create coverages for all spike-in genes.")
parser.add_argument('--experiment-name',action = 'store', default="", help = "Specify name of an experiment.", required = False, nargs = 1)
args = parser.parse_args()


modules_to_execute = ["coverage", "concentration"]

input_dict = read_sample_sheet("/home/tdrozd/development/sirv-suite/test/input/sample_sheet_merged.tsv", modules_to_execute = modules_to_execute)
out = "/home/tdrozd/development/sirv-suite/standalone/output/"

if "concentration" in modules_to_execute:
    c = ERCCcorrelation(sample_sheet=input_dict["concentration"], output_dir=out)

    c = SIRVsuiteConcentration(sample_sheet=input_dict["concentration"], output_dir=out)
    c.create_sirvsuite_boxplot(c.data)
    c.create_sirvsuite_heatmap(c.data)

if "coverage" in modules_to_execute:
    k = SIRVsuiteCoverage(sample_sheet=input_dict["coverage"], output_dir=out, experiment_name = "")
    k.expected_coverage()
    k.bam_to_coverage()
    k.calc_statistics()
    k.coverage_plot()

"""

import pyBigWig as bw

bb = bw.open("/home/tdrozd/test_data/out/coverage/bigwig/NGS2.54_0002_antisense.bw")

print (bb.chroms())
print (bb.header())

