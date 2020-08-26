import pysam as ps
import pandas as pd
import numpy as np
import pyBigWig as bwig
import pyranges
from scipy.stats import norm
import os
import re
import copy

class SIRVsuiteCoverage():

    def __init__(self, gene_list = ["SIRV1","SIRV2","SIRV3","SIRV4","SIRV5","SIRV6","SIRV7"]):
        self.verbose = "DEBUG"
        self.target_gene_id = gene_list
        self._strands = ["+","-"]
        self.output_path = "/home/tdrozd/development/sirv-suite/test/"
        self.annotation_path = "/home/tdrozd/data/SIRVome_SIRVset4/SIRV_isoforms_multi-fasta-annotation_190418a.gtf"
        print ("SIRVsuite coverage creator initialized")

    def __values2intervals__(self, data, start_pos = 0):
        tmp = np.logical_and(data[1:] != data[:-1], np.logical_or(np.isnan(data[1:]) == False, np.isnan(data[:-1]) == False))
        starts = np.append(np.array([start_pos]), np.where(tmp)[0] + start_pos + 1)
        ends = np.append(starts[1:], start_pos + len(data))
        values = data[starts - start_pos]
        return (starts, ends, values)

    def __strand2text__(self, text):
        if (text == "+"):
            return "sense"
        elif (text == "-"):
            return "antisense"

    def create_nested_dict(input_list, init = 0):
        # a recursive approach to create and initialize nested dictionary
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

    def __cov_data_export__(self, output_path, output_type = "bigwig"):
        # function used for exporting coverage data to "bigwig" or "cov" file
        if (hasattr(self, "bam_coverage")):
            cov_dict = self.bam_coverage
            
            for strand in self._strands:
                for sample in cov_dict.keys():
                    if (output_type == "bigwig"):
                        with bwig.open(os.path.join(output_path, sample+"_"+self.__strand2text__(strand)+".bw"), "w") as bw:
                            header = [(gene, len(cov_dict[sample][gene][strand])) for gene in sorted(cov_dict[sample].keys())] 
                            bw.addHeader(header)
                            for gene in sorted(cov_dict[sample].keys()):
                                (starts, ends, values) = self.__values2intervals__(cov_dict[sample][gene][strand])
                                chroms = np.array([gene] * len(values))
                                bw.addEntries(chroms, starts, ends = ends, values = values)

                    elif (output_type == "cov"):
                        for gene in sorted(cov_dict[sample].keys()):
                            with open(os.path.join(output_path, sample+"_"+gene+"_"+self.__strand2text__(strand)), "w") as cov:
                                cov.write("position"+"\t"+"coverage"+"\n")
                                coverage = cov_dict[sample][gene][strand]
                                for idx in range(len(coverage)):
                                    cov.write(str(idx+1)+"\t"+str(coverage[idx])+"\n") 
        else:
            print ("gfdgfdgdfg")
    
    def bam_to_coverage(self, sample_dict = None, output_type = "bigwig"):

        if (sample_dict == None):
            raise Warning("sample_dict must be specified")
            return 

        if self.verbose == "DEBUG":
            print ("Calculating real coverage from input files")

        stats = {"num_reads":0,
                "CoD":0,
                "scaling_factor":0}

        bam_coverage = {} 
        stat_dict = {}
        
        for sample in sample_dict.keys():

            strandeness = sample_dict[sample]["lib_strand"]
            bam_path = sample_dict[sample]["bam"]

            if (not os.path.exists(bam_path)):
                print (bam_path+"does not exist.. skipping to the next file")

            bamFile = ps.AlignmentFile(bam_path,"rb")

            bam_coverage[sample] = {}
            stat_dict[sample] = {}

            for gene in self.target_gene_id:

                length = [bamFile.lengths[i] for i in range(len(bamFile.references)) if bamFile.references[i] == gene][0] 
                bam_coverage[sample][gene] = {}
                stat_dict[sample][gene] = {}

                if (strandeness == "none"):
                    bam_coverage[sample][gene]["none"] = np.zeros(length)
                    stat_dict[sample][gene]["none"] = copy.deepcopy(stats)

                elif (strandeness in ["rev","fwd"]):
                    bam_coverage[sample][gene]["+"] = np.zeros(length)
                    bam_coverage[sample][gene]["-"] = np.zeros(length)
                    stat_dict[sample][gene]["+"] = copy.deepcopy(stats)
                    stat_dict[sample][gene]["-"] = copy.deepcopy(stats)

                for fragmentRead in bamFile.fetch(gene):

                    if (strandeness == "none"):
                        strand = strandeness

                    elif(strandeness in ["rev","fwd"]):
                        # inferring strandeness of a read
                        if ((fragmentRead.is_paired and fragmentRead.is_read1 and fragmentRead.is_reverse ) or 
                            (fragmentRead.is_paired and fragmentRead.is_read2 and not fragmentRead.is_reverse) or 
                            (not fragmentRead.is_paired and not fragmentRead.is_reverse)):

                            strand_test = "+"
                        else:
                            strand_test = "-"

                        # swap strands if reverse library prep
                        if (strandeness == "rev" and strand_test == "+"):
                            strand_test = "-"
                        elif (strandeness == "rev" and strand_test == "-"):
                            strand_test = "+"

                    stat_dict[sample][gene][strand_test]["num_reads"] += 1
                    bam_coverage[sample][gene][strand_test][fragmentRead.positions] += 1

        self.cov_stats = stat_dict
        self.bam_coverage = bam_coverage

        self.__cov_data_export__(self.output_path)

    def _path_features(self, path):
        matching_pattern = "((?:(?:.*\\/)*(.*)\\/))*(?:(.*)\\.(.*))*"
        feature_match = re.match(matching_pattern, filePath)
        out_dict = {"parent_dir": feature_match.group(2),
                    "file_name": feature_match.group(3),
                    "extension": feature_match.group(4),
                    "path": feature_match.group(1)}
        return out_dict


    def load_annotation(self, annotation_path):
        annotation_type = self._path_features(annotation_path)["extension"]
        
        if (annotation_type == "gtf"):
            pyranges.read_gtf(annotation_path, as_df=True)
        elif (annotation_type == "bed"):
            pyranges.read_bed(annotation_path, as_df=True)



test_dict = {"test_sample":{"bam":"/home/data/slurm_scratch/tdrozd/DAP/DEV/67/output/Aligned.sortedByCoord.out.bam",
                                "lib_strand": "rev"}}

samples = ["sample1","sample2","sample3"]
genes = ["gene1","gene2","gene3","gene4"]
strands = ["+","-"]
stat_attrs = ["read_cnts","CoD","scaling_factor"]

l = [samples,genes,strands]

k = SIRVsuiteCoverage()
k.bam_to_coverage(test_dict)
