from Pipeline.Coverage.SIRVsuiteCoverage import SIRVsuiteCoverage
from Pipeline.Concentration.SIRVsuite_concentration import *
import pyranges
import numpy as np

def create_testing_dataset(annotation=None,output_path="./",mode='cufflinks'):
    annot_table = pyranges.read_gtf(annotation, as_df= True)
    annot_table["FPKM"] = np.random.normal(size=len(annot_table))
    annot_ranges = pyranges.PyRanges(annot_table)
    annot_ranges.to_gtf(path = output_path)
"""
test_dict = {"test_sample":{"bam":"/home/data/slurm_scratch/tdrozd/DAP/DEV/67/output/Aligned.sortedByCoord.out.bam",
                                "lib_strand": "rev"}}

samples = ["sample1","sample2","sample3"]
genes = ["gene1","gene2","gene3","gene4"]
strands = ["+","-"]
stat_attrs = ["read_cnts","CoD","scaling_factor"]

l = [samples,genes,strands]

k = SIRVsuiteCoverage()
k.bam_to_coverage(test_dict)

"""

"""
def create_nested_dict(input_list, init = 0):

    input = input_list[0]
    input_list.pop(0)
    l = len(input_list)

    if (l > 0):
        temp = createNestedDict(input_list)
        temp = [copy.deepcopy(temp) for i in range(len(input))]
        return dict(zip(input, temp))

    elif (l == 0):
        n = len(input)
        init_list = [init] * n
        return dict(zip(input, init_list))
"""

create_new_dataset = False

if (create_new_dataset == True):
    create_testing_dataset(
        annotation="/home/tdrozd/data/SIRVome_SIRVset4/SIRV_isoforms_multi-fasta-annotation_190418a.gtf",
        output_path = "/home/tdrozd/development/sirv-suite/test/test_sirv_cufflinks.gtf")

file_dict = load_sample_sheet("/home/tdrozd/development/sirv-suite/test/input/sheet_cufflinks.tsv")
relative_abundance = get_relative_abundance(file_dict)

export_data(relative_abundance, output_path = "/home/tdrozd/development/sirv-suite/test/")
