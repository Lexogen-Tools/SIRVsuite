from Pipeline.Coverage.SIRVsuiteCoverage import SIRVsuiteCoverage
import copy

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

#final = createNestedDict(l)
