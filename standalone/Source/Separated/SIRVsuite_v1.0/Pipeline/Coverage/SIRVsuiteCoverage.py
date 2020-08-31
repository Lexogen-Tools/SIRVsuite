import pysam as ps
import pandas as pd
import numpy as np
import pyBigWig as bwig
import pyranges
from scipy.stats import norm
import os
import re
import copy
import warnings

class SIRVsuiteCoverage():
    """
    SIRVsuiteCoverage class for creating coverage plots & statistics

    The class is concepted to be used in following steps:
        1. load annotation file, which contains**
        2. calculate expected coverage based on the annotation
        3. calculate coverage from an alignment file
        4. calculate statistics
        5. create coverage plots

    ** by default, the annotation is loaded from a Resources folder in SIRVsuite project 
    (annotation is identical to the one which can be found https://www.lexogen.com/sirvs/download/ under section SIRV-Set 4)

    """
    def __init__(self, gene_list = ["SIRV1","SIRV2","SIRV3","SIRV4","SIRV5","SIRV6","SIRV7"]):
        self.verbose = "DEBUG"
        self.target_gene_id = gene_list
        self.output_path = "/home/tdrozd/Lexogen/ssuite/out/"
        self.annotation_path = "/home/tdrozd/data/SIRVome_SIRVset4/SIRV_isoforms_multi-fasta-annotation_190418a.gtf"
        
        ## constants
        self._strands = ["-","+"]
        #self._strandNames = {"+":"sense","-":"antisense"}
        
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

    def __create_nested_dict__(self, input_list, init = 0):
        # a recursive approach to create and initialize nested dictionary
        input = input_list[0]
        input_list.pop(0)
        l = len(input_list)

        if (l > 0):
            temp = __create_nested_dict__(input_list)
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
            warnings.warn("No coverage detected to export..")
    
    def bam_to_coverage(self, sample_dict = None, output_type = "bigwig"):

        if (sample_dict == None):
            raise ValueError("Sample_dict must be specified")
            return 

        if self.verbose == "DEBUG":
            print ("Calculating real coverage from input files")

        stats = {"num_reads": 0}

        bam_coverage = {} 
        stat_dict = {}
        
        for sample in sample_dict.keys():

            strandeness = sample_dict[sample]["lib_strand"]
            bam_path = sample_dict[sample]["bam"]

            if (not os.path.exists(bam_path)):
                if self.verbose == "debug":
                    print (bam_path+"does not exist.. skipping to the next file")
            else:
                index_path = bam_path + ".bai"

                if (not os.path.exists(index_path)):
                    ps.index(index_path)
                    if self.verbose == "debug":
                        print (".bai index file not detected.. creating index file..")

            bamFile = ps.AlignmentFile(bam_path,"rb")

            bam_coverage[sample] = {}
            stat_dict[sample] = {}

            annotation_df = self.annotation_df["whole"]

            for gene in self.target_gene_id:

                if (hasattr(self, "annotation_df")):
                    contig_id = annotation_df["Chromosome"][annotation_df["gene_id"] == gene].values[0]
                else:
                    contig_id = gene

                if (gene == "ERCC-00002"):
                    print ()
         
                bam_coverage[sample][gene] = {}
                stat_dict[sample][gene] = {}

                # check if annotation processed
                if (hasattr(self, "gene_coords")):
                    bam_gene_start = self.gene_coords[gene][0]
                    bam_gene_end = self.gene_coords[gene][1]
                    length = bam_gene_end - bam_gene_start + 1
                else:
                    length = [bamFile.lengths[i] for i in range(len(bamFile.references)) if bamFile.references[i] == gene][0]
                    bam_gene_start = 1
                    bam_gene_end = length

                if (strandeness == "none"):
                    bam_coverage[sample][gene]["none"] = np.zeros(length)
                    stat_dict[sample][gene]["none"] = copy.deepcopy(stats)

                elif (strandeness in ["rev","fwd"]):
                    bam_coverage[sample][gene]["+"] = np.zeros(length)
                    bam_coverage[sample][gene]["-"] = np.zeros(length)
                    stat_dict[sample][gene]["+"] = copy.deepcopy(stats)
                    stat_dict[sample][gene]["-"] = copy.deepcopy(stats)

                for fragmentRead in bamFile.fetch(contig_id, bam_gene_start, bam_gene_end):

                    if (strandeness == "none"):
                        strand = strandeness

                    elif(strandeness in ["rev","fwd"]):
                        # inferring strandeness of a read
                        if ((fragmentRead.is_paired and fragmentRead.is_read1 and fragmentRead.is_reverse ) or 
                            (fragmentRead.is_paired and fragmentRead.is_read2 and not fragmentRead.is_reverse) or 
                            (not fragmentRead.is_paired and not fragmentRead.is_reverse)):

                            strand = "+"
                        else:
                            strand = "-"

                        # swap strands if reverse library prep
                        if (strandeness == "rev" and strand == "+"):
                            strand = "-"
                        elif (strandeness == "rev" and strand == "-"):
                            strand = "+"

                    stat_dict[sample][gene][strand]["num_reads"] += 1
                    positions = np.array(fragmentRead.positions)
                    positions = positions[positions<bam_gene_end] - bam_gene_start
                    bam_coverage[sample][gene][strand][positions] += 1
                
                """
                except:
                    print ("there was an error with extracting fragments for specified contig.. skipping")
                    continue
                """

        self.cov_stats = stat_dict
        self.bam_coverage = bam_coverage

        self.__cov_data_export__(self.output_path)
        

    def _path_features(self, path):
        matching_pattern = "((?:(?:.*\\/)*(.*)\\/))*(?:(.*)\\.(.*))*"
        feature_match = re.match(matching_pattern, path)
        out_dict = {"parent_dir": feature_match.group(2),
                    "file_name": feature_match.group(3),
                    "extension": feature_match.group(4),
                    "path": feature_match.group(1)}
        return out_dict


    def load_annotation(self, annotation_path):
        annotation_type = self._path_features(annotation_path)["extension"]
        
        try:
            if self.verbose == "DEBUG":
                print ("loading "+annotation_path)
            if (annotation_type == "gtf"):
                annotation_df = pyranges.read_gtf(annotation_path, as_df=True)
            elif (annotation_type == "bed"):
                annotation_df = pyranges.read_bed(annotation_path, as_df=True)

            # limit only to a feature of interest
            annotation_df[annotation_df["Feature"]=="exon"] 
        except:
            raise ValueError("An error occured while loading annotation..")
        
        self.annotation_df = {"whole": annotation_df}

    def get_continous_coverage_ends(self, starts, ends):
        ## This function returnes continous segment of oveerlapping features (exons) from list of starts and ends of the features 
        sorted_idx = np.argsort(starts)

        contig_start = [starts[sorted_idx[0]]]
        contig_end = []

        for i in range(0,len(sorted_idx)-1):

            prev_idx = ends[sorted_idx[i]]
            next_idx = starts[sorted_idx[i+1]]

            if (next_idx > prev_idx):
                contig_end.append(prev_idx)
                contig_start.append(next_idx)
            
        contig_end.append(ends[sorted_idx[-1]])

        return contig_start, contig_end

    def calc_start_transition(self,transcriptLength,mean=200,std=80):
        ## This function estimates probability distribution of a fragment according to its length, mean and std of the model distribution
        
        lengthProbability = norm.pdf(np.arange(1,transcriptLength+1),mean,std)
        
        lengthProbability /= np.sum(lengthProbability)
        
        fragmentProbability = np.zeros((transcriptLength,transcriptLength))
        
        for length in np.arange(1,transcriptLength):
            fragmentProbability[length,0:transcriptLength-length+1] = lengthProbability[length]/(transcriptLength-length+1)
        
        startDistance = np.sum(fragmentProbability,axis=0)
        
        adjustedLength = np.sum(lengthProbability*(transcriptLength-np.arange(1,transcriptLength+1)+1))
        
        startDistance /= np.max(startDistance)
        
        return startDistance, lengthProbability, adjustedLength

    def calc_statistics(self):
        if (not hasattr(self, "cov_stats") or not hasattr(self, "bam_coverage")):
            raise ValueError("Measured coverage not detected")
        elif (not hasattr(self, "expected_coverage")):
            raise ValueError("Expected coverage not detected")

        for sample in self.bam_coverage.keys():
            for gene in self.bam_coverage[sample].keys():
                for strand in self.expected_coverage[gene].keys():
                    expected_cov = self.expected_coverage[gene][strand]
                    measured_cov = self.bam_coverage[sample][gene][strand]
                    self.cov_stat
                    self.cov_stats[sample][gene][strand]["CoD"] = self.CoD(measured_cov, expected_cov)

    def CoD(self,real_cov,expect_cov):
        ## This function calculate CoD from 2 vectors (real and expected coverage)
        ## and returns coefficient for scaling and CoD

        if (len(real_cov) != len(expect_cov)):
            raise ValueError("real and expected coverage must have the same length")

        
        sumBAC_theory = np.sum(expect_cov) # sumBAC (sum of base annotation count), sum of expected coverage
        
        sumCoverage_real = np.sum(real_cov) # sum of real coverage

        scaling_coefficient = sumCoverage_real/sumBAC_theory # calculation of scaling factor for real coverage (for comparison)
        
        ## Correction of 0 factor
        if ( scaling_coefficient == 0 ): 
            scaling_coefficient = 1
            
        ## Scaling
        coverageReal_scaled = real_cov/scaling_coefficient
            
        ## COD computation 
        difference = expect_cov - coverageReal_scaled
        COD = np.sum((difference)**2)/sumBAC_theory
        COD = np.round(COD,4)
        
        return COD, scaling_coefficient

    def annot_2_UTR(self, region_length = 200):
        ## This function modifies annotation data frame and trims the transcripts to a length close to 3' (where QuantSeq is mostly expressed)

        if (not hasattr(self, "annotation_df")):
            warnings.warn("There was no annotation detected.. skipping UTR mode..")
            return
        else:
            annotation_df = self.annotation_df["whole"]
        
        UTR_annotation = annotation_df[0:0] # Creating an empty data frame with the same structure as annotation data frame

        for strand in self._strands:
            transcripts = np.unique(annotation_df[annotation_df["Strand"]  ==  strand]["transcript_id"])

            for transcript in transcripts:

                annot_transcript = annotation_df[annotation_df["transcript_id"]  ==  transcript]

                if ( strand  ==  '+'):
                    annot_transcript = annot_transcript.sort_values(by=["Start"], axis = 0, ascending = False).reset_index(drop=True)
                    exon_lengths = annot_transcript["End"] - annot_transcript["Start"]

                    remainder = region_length
                    UTR_data_frame = annot_transcript

                    for index, length in enumerate(exon_lengths):
                        remainder = remainder - length
                        if (remainder < 0):
                            UTR_data_frame.loc[index,"Start"] = UTR_data_frame.loc[index,"Start"] - remainder
                            UTR_data_frame = UTR_data_frame.loc[0:index,:]
                            break

                elif ( strand  ==  '-'):
                    annot_transcript = annot_transcript.sort_values(by=["Start"], axis = 0, ascending = True).reset_index(drop=True)
                    exon_lengths = annot_transcript["End"] - annot_transcript["Start"]

                    remainder = region_length
                    UTR_data_frame = annot_transcript

                    for index, length in enumerate(exon_lengths):
                        remainder = remainder - length
                        if ( remainder < 0 ):
                            UTR_data_frame.loc[index,"End"] = UTR_data_frame.loc[index,"End"] + remainder
                            UTR_data_frame = UTR_data_frame.loc[0:index,:]
                            break
                
                UTR_annotation = UTR_annotation.append(UTR_data_frame)

        UTR_annotation = UTR_annotation.reset_index(drop=True)

        self.annotation_df["UTR"] = UTR_annotation


    def expected_coverage(self, transition_lengths = None):
        if (not hasattr(self, "annotation_df")):
            raise ValueError("No annotation detected for theoretical coverage.. Please load annotation in .gtf or .bed format.")
        
        expected_coverage = dict()

        self.gene_coords = {}

        for mode in self.annotation_df.keys():

            annotation_df = self.annotation_df[mode]

            self.target_gene_id = np.unique(annotation_df["gene_id"])

            for gene in self.target_gene_id:
                gene_annot = annotation_df[annotation_df["gene_id"]==gene]

                end_coord = np.max(gene_annot["End"])
                start_coord = np.min(gene_annot["Start"])

                self.gene_coords[gene] = (start_coord, end_coord) 
                
                expected_coverage[gene] = {}
                
                for strand in self._strands:
                    
                    stranded_gene_annot = gene_annot[gene_annot["Strand"] == strand]

                    # skip when no record for gene and strand
                    if (len(stranded_gene_annot) == 0):
                        continue

                    # declare base for coverage
                    expected_coverage[gene][strand] = np.zeros(end_coord-start_coord+1)

                    for transcript in np.unique(stranded_gene_annot["transcript_id"]):
                        
                        transcript_annot = stranded_gene_annot[stranded_gene_annot["transcript_id"] == transcript].copy(deep=True).reset_index(drop = True)                    

                        exon_lengths = np.array(transcript_annot["End"]) - np.array(transcript_annot["Start"]) + 1

                        if (transition_lengths is not None):
                            weights,_,_ = self.calc_start_transition(np.sum(exon_lengths), mean = transition_lengths[0], std = transition_lengths[1])
                        else:
                            weights = np.ones(np.sum(exon_lengths))

                        weight_idx = 0

                        for row_idx,row in transcript_annot.iterrows():

                            shifted_start = row["Start"] - start_coord
                            shifted_end = row["End"] - start_coord + 1

                            # apply weights for transition lengths for transcript start and end
                            weights = weights * weights[::-1]
                            expected_coverage[gene][strand][shifted_start:shifted_end] += weights[weight_idx:(weight_idx + exon_lengths[row_idx])]

                            weight_idx += exon_lengths[row_idx]

        self.expected_coverage = expected_coverage


test_dict = {"test_sample":{"bam":"/home/tdrozd/Lexogen/ssuite/NGS2.59/NGS2.59_0045/Aligned.sortedByCoord.out.bam",
                                "lib_strand": "rev"}}

samples = ["sample1","sample2","sample3"]
genes = ["gene1","gene2","gene3","gene4"]
strands = ["+","-"]
stat_attrs = ["read_cnts","CoD","scaling_factor"]

l = [samples,genes,strands]

k = SIRVsuiteCoverage()
k.load_annotation("./Resources/SIRV_isoforms_multi-fasta-annotation_190418a.gtf")
#k.annot_2_UTR(200)
k.expected_coverage(transition_lengths=(25,30))
k.bam_to_coverage(test_dict)
k.calc_statistics()
print ()
