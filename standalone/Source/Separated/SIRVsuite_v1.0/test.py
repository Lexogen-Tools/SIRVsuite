from Pipeline.Coverage.SIRVsuiteCoverage import SIRVsuiteCoverage
from Pipeline.Concentration.SIRVsuite_concentration import SIRVsuiteConcentration
from Pipeline.Correlation.ERCC_correlation import 
from Pipeline.helper import *
import pyranges
import numpy as np

def create_testing_dataset(annotation=None,output_path="./",mode='cufflinks'):
    annot_table = pyranges.read_gtf(annotation, as_df= True)
    annot_table["FPKM"] = np.random.normal(size=len(annot_table))
    annot_ranges = pyranges.PyRanges(annot_table)
    annot_ranges.to_gtf(path = output_path)

input_dict = read_sample_sheet("/home/tdrozd/development/sirv-suite/test/input/sample_sheet_merged.tsv")

#k = SIRVsuiteCoverage(sample_sheet=input_dict["coverage"], output_dir="/home/tdrozd/development/sirv-suite/", experiment_name = "")
#k.expected_coverage(transition_lengths=(25,30))
#k.bam_to_coverage(test_dict)
#k.calc_statistics()
#print ()


a = get_relative_abundance(input_dict["concentration"])
b = create_sirvsuite_heatmap(a, output_path = "/home/tdrozd/development/sirv-suite/test")