import pysam as ps
import pandas as pd
import numpy as np
import pyBigWig as bwig
import pyranges
from scipy.stats import norm
import os
import copy
import warnings
from ..helper import *
import sys
from Pipeline.Coverage.CairoDrawer import CairoDrawer

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
    def __init__(self, sample_sheet = None, output_dir = "./", gene_list = ["SIRV1","SIRV2","SIRV3","SIRV4","SIRV5","SIRV6","SIRV7"], experiment_name = ""):

        self.verbose = "DEBUG"
        self.target_gene_id = gene_list

        if output_dir == None:
            self.output_path = "/".join(__file__.split("/")[:-6]) + "/coverage"
        else:
            self.output_path = output_dir

        self.annotation_path = "/".join(__file__.split("/")[:-3]) + "/Resources/SIRVsuite_annotation.gtf"
        self.load_annotation(self.annotation_path)

        self.experiment_name = experiment_name

        if sample_sheet == None:
            raise ValueError("please specify path to the sample sheet..")
        else:
            self.sample_dict = sample_sheet
        
        ## constants
        self._strands = ["-","+"]
        
        print ("SIRVsuite coverage creator initialized")

    def __values2intervals__(self, data, start_pos = 0):
        """
        A method to convert data array (coverage) into intervals specified by starts, ends and values
        """
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

    def __cov_data_export__(self, output_path, output_type = "bigwig"):
        # function used for exporting coverage data to "bigwig" or "cov" file
        if (hasattr(self, "bam_coverage")):
            cov_dict = self.bam_coverage
            
            for strand in self._strands:
                for sample in cov_dict.keys():

                    if (output_type.lower() == "bigwig"):
                        with bwig.open(os.path.join(output_path, sample+"_"+self.__strand2text__(strand)+".bw"), "w") as bw:
                            header = [(gene, len(cov_dict[sample][gene][strand])) for gene in sorted(cov_dict[sample].keys()) if strand in cov_dict[sample][gene].keys()] 
                            bw.addHeader(header)
                            for gene in sorted(cov_dict[sample].keys()):
                                (starts, ends, values) = self.__values2intervals__(cov_dict[sample][gene][strand])
                                chroms = np.array([gene] * len(values))
                                bw.addEntries(chroms, starts, ends = ends, values = values)

                    elif (output_type.lower() == "cov"):
                        for gene in sorted(cov_dict[sample].keys()):
                            with open(os.path.join(output_path, sample+"_"+gene+"_"+self.__strand2text__(strand)), "w") as cov:
                                cov.write("position"+"\t"+"coverage"+"\n")
                                coverage = cov_dict[sample][gene][strand]
                                for idx in range(len(coverage)):
                                    cov.write(str(idx+1)+"\t"+str(coverage[idx])+"\n") 
        else:
            warnings.warn("No coverage detected to export..")

    def set_output_dir(self, output_dir, create = False):

        if ((not create) & (not os.path.exists(output_dir))):
            raise ValueError("Provided folder does not exist.. please use create=True argument to make a new directory or change to valid path..")
            return 
        elif ((create) & (not os.path.exists(output_dir))):
            os.makedirs(output_dir)
            print ("New directories created.")

        self.output_path = output_dir
        print ("output directory changed to: "+self.output_path)
    
    def bam_to_coverage(self, output_type = "bigwig"):
        """
        A method for loading coverage from bam files
        """
        
        sample_dict = self.sample_dict

        if self.verbose == "DEBUG":
            print ("Calculating real coverage from input files")

        stats = {"num_reads": 0}

        bam_coverage = dict()
        stat_dict = dict()
        
        for sample in sample_dict.keys():

            strandeness = sample_dict[sample]["read_orientation"]
            bam_path = sample_dict[sample]["alignment_path"]

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

            bam_coverage[sample] = dict()
            stat_dict[sample] = dict()

            annotation_df = self.annotation_df["whole"]

            for gene in self.target_gene_id:

                if (hasattr(self, "annotation_df")):
                    contig_id = annotation_df["Chromosome"][annotation_df["gene_id"] == gene].values[0]
                else:
                    contig_id = gene
         
                bam_coverage[sample][gene] = dict()
                stat_dict[sample][gene] = dict()

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
                    positions = positions[(positions<bam_gene_end) & (positions>bam_gene_start)] - bam_gene_start
                    bam_coverage[sample][gene][strand][positions] += 1 

        self.cov_stats = stat_dict
        self.bam_coverage = bam_coverage

        # export coverage to .bw format
        bigwig_path = os.path.join(self.output_path, "coverage/bigwig")
        if (not os.path.exists(bigwig_path)):
            os.makedirs(bigwig_path)

        self.__cov_data_export__(bigwig_path, output_type="bigWig")        

    def load_annotation(self, annotation_path):
        """
        A method for loading annotation using pyranges library
        """
        annotation_type = path_features(annotation_path)["extension"]
        
        try:
            if self.verbose == "DEBUG":
                print ("loading "+annotation_path)
            if (annotation_type == "gtf"):
                annotation_df = pyranges.read_gtf(annotation_path, as_df=True)
            elif (annotation_type == "bed"):
                annotation_df = pyranges.read_bed(annotation_path, as_df=True)

            # limit only to a feature of interest
            annotation_df[annotation_df["Feature"]=="exon"]
            # trim to required cols only 
            annotation_df = annotation_df[["Chromosome", "Start", "End", "Strand", "gene_id", "transcript_id"]]
        except:
            raise ValueError("An error occured while loading annotation..")
        
        self.annotation_df = {"whole": annotation_df}

    def get_continous_coverage_ends(self, starts, ends):
        """
        The method returnes continous segment of oveerlapping features (exons) from list of starts and ends of the features
        """
        
        sorted_idx = np.argsort(starts)

        contig_start = [starts[sorted_idx[0]]]
        contig_end = []
        thre = 0
        
        for i in range(0,len(sorted_idx)-1):

            prev_idx = ends[sorted_idx[i]]
            next_idx = starts[sorted_idx[i+1]]

            if prev_idx > thre:
                thre = prev_idx

            if (next_idx > thre):
                contig_end.append(prev_idx)
                contig_start.append(next_idx)
            
        contig_end.append(ends[sorted_idx[-1]])

        return contig_start, contig_end

    def calc_start_transition(self,transcriptLength,mean=200,std=80):
        """
        The method estimates probability distribution of a fragment according to its length, mean and std of the model distribution
        """
        
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
        """
        The method checks for bam coverages and expected coverages and then calculates CoD and scaling factor for particular samples, genes and strands 
        """

        # check if all attributes present
        if (not hasattr(self, "cov_stats") or not hasattr(self, "bam_coverage")):
            raise ValueError("Measured coverage not detected")
        elif (not hasattr(self, "expected_coverage")):
            raise ValueError("Expected coverage not detected")
        
        CoD_table_path = os.path.join(self.output_path, "coverage/") 
        if not os.path.exists(CoD_table_path):
            os.makedirs(CoD_table_path)
        
        experiment_name = self.experiment_name
        if experiment_name != "":
            experiment_name = "_" + experiment_name

        for mode in self.annotation_df.keys():
            with open(CoD_table_path+"CoD_table"+experiment_name+"_"+mode+".tsv","w") as CoD_table:
                header = ";" + ";".join(list(sorted(self.bam_coverage.keys()))) + "\n"
                CoD_table.write(header)

                rows = dict()
                
                for sample in sorted(self.bam_coverage.keys()):
                    for gene in self.bam_coverage[sample].keys():
                        for strand in self.expected_coverage[mode][gene].keys():
                            expected_cov = self.expected_coverage[mode][gene][strand]
                            measured_cov = self.bam_coverage[sample][gene][strand]

                            # calculate CoD and scaling factor and assign them to the cov_stats dictionary  
                            self.cov_stats[sample][gene][strand]["CoD"], self.cov_stats[sample][gene][strand]["scale"] = self.CoD(measured_cov, expected_cov) 
                            
                            row_key = gene+"("+strand+")"
                            if row_key not in rows.keys():
                                rows[row_key] = str(self.cov_stats[sample][gene][strand]["CoD"])
                            else:
                                rows[row_key] += ";" + str(self.cov_stats[sample][gene][strand]["CoD"])
            
                for row_idx in rows.keys():
                    CoD_table.write(row_idx + ";" + rows[row_idx]+"\n") 
                   

    def CoD(self,real_cov,expect_cov):
        """
        The method calculate CoD from 2 vectors (real and expected coverage) and returns coefficient for scaling and CoD

        """

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
        """
        This method modifies annotation data frame and trims the transcripts to a length close to 3' (where QuantSeq is mostly expressed)
        
        """

        # check if annotation loaded
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
                
                # initiliaze parameters for a particular strand
                if ( strand  ==  '+'):
                    ascending = False
                    const = -1
                    UTR_df_col = "Start"
                elif ( strand  ==  '-'):
                    ascending = True
                    const = 1
                    UTR_df_col = "End"

                # sort annotation dataframe by Start column
                annot_transcript = annot_transcript.sort_values(by=["Start"], axis = 0, ascending = ascending).reset_index(drop=True)
                exon_lengths = annot_transcript["End"] - annot_transcript["Start"]
                remainder = region_length
                UTR_data_frame = annot_transcript
                
                # 
                for index, length in enumerate(exon_lengths):
                    remainder = remainder - length
                    if (remainder < 0):
                        UTR_data_frame.loc[index,UTR_df_col] = UTR_data_frame.loc[index,UTR_df_col] + remainder*const
                        UTR_data_frame = UTR_data_frame.loc[0:index,:]
                        break
                
                UTR_annotation = UTR_annotation.append(UTR_data_frame)


        UTR_annotation = UTR_annotation.reset_index(drop=True)


        self.annotation_df["UTR_test"] = UTR_annotation


    def expected_coverage(self, transition_lengths = False, tlen_model_param = (0,0), mode = "normal", region_length = 200):
        """
        Method to create expected coverage based on provided annotation.
        
        Experimental feature: transition lengths modelling
            input: list, ndarray or tuple of a length 2, specifying model parameters in the following order [mean, std], e.g. [100, 50] 

        Returns: 3-level dictionary gene_id -> strand -> expected_coverage
        """

        # check if annotation has been loaded
        if (not hasattr(self, "annotation_df")):
            raise ValueError("No annotation detected for theoretical coverage.. Please load annotation in .gtf or .bed format.")
        
        if (len(tlen_model_param) != 2):
            raise ValueError("The length of passed array 'transition_lengths' does no correspond model requirement..")
    
        if mode == "UTR":
            self.annot_2_UTR(region_length)

        expected_coverage = dict()

        self.gene_coords = dict()

        for mode in self.annotation_df.keys():

            annotation_df = self.annotation_df[mode]

            #self.target_gene_id = np.unique(annotation_df["gene_id"])
            expected_coverage[mode] = dict()

            for gene in self.target_gene_id:
                gene_annot = annotation_df[annotation_df["gene_id"]==gene]

                if (len(gene_annot) == 0):
                    warnings.warn("No record for gene "+gene+" has been found in the provided annotation file.. skipping gene..")

                end_coord = np.max(gene_annot["End"])
                start_coord = np.min(gene_annot["Start"])

                self.gene_coords[gene] = (start_coord, end_coord) 
                
                expected_coverage[mode][gene] = dict()
                
                for strand in self._strands:
                    
                    stranded_gene_annot = gene_annot[gene_annot["Strand"] == strand]

                    # skip when no record for gene and strand
                    if (len(stranded_gene_annot) == 0):
                        continue

                    # initialize base for coverage with 0s 
                    expected_coverage[mode][gene][strand] = np.zeros(end_coord-start_coord+1)

                    for transcript in np.unique(stranded_gene_annot["transcript_id"]):
                        
                        transcript_annot = stranded_gene_annot[stranded_gene_annot["transcript_id"] == transcript].copy(deep=True).reset_index(drop = True)                    

                        exon_lengths = np.array(transcript_annot["End"]) - np.array(transcript_annot["Start"]) + 1

                        # if parameters for transition transcript lengths are defined, use probabilistic model to calculate multiplicative weights for starts and ends,
                        # otherwise set weights to 1
                        if ( transition_lengths ):
                            if (tlen_model_param == (0,0)):
                                tlen_model_param = (25,30)

                            weights,_,_ = self.calc_start_transition(np.sum(exon_lengths), mean = tlen_model_param[0], std = tlen_model_param[1])
                        else:
                            weights = np.ones(np.sum(exon_lengths))

                        weight_idx = 0

                        for row_idx,row in transcript_annot.iterrows():

                            shifted_start = row["Start"] - start_coord
                            shifted_end = row["End"] - start_coord + 1

                            # apply weights for transition lengths for transcript start and end
                            weights = weights * weights[::-1]
                            expected_coverage[mode][gene][strand][shifted_start:shifted_end] += weights[weight_idx:(weight_idx + exon_lengths[row_idx])]

                            weight_idx += exon_lengths[row_idx]

        self.expected_coverage = expected_coverage

    def coverage_plot(self):
        
        # define base sizes, coordinates
        page_width = 2000
        page_height = 1324

        title_y = 40
        title_size = 40

        panel_gap_y = title_size 

        header_y = title_size + panel_gap_y
        header_height = page_height / 6

        exon_panel_x = page_width / 28
        exon_panel_width = page_width - 2*exon_panel_x
        exon_panel_height = page_height / 3

        transcript_line_offset_x = exon_panel_x + exon_panel_width * 0.1
        
        transcript_line_width = exon_panel_width * 0.8

        transcript_text_x = exon_panel_x + (transcript_line_offset_x - exon_panel_x)/2
        coverage_panel_x = exon_panel_x
        
        coverage_panel_width = exon_panel_width
        coverage_panel_height = exon_panel_height 
        
        if "UTR" in self.annotation_df.keys():
            mode = "UTR"
        else:
            mode = "whole"

        if self.experiment_name != "":
            tab = {"Experiment":self.experiment_name}
        else:
            tab = {}
            header_height = header_height*3/4

        for gene in self.target_gene_id:
            
            max_expect_covs = max([max(self.expected_coverage["whole"][gene][s]) for s in self.expected_coverage["whole"][gene].keys()])
            
            gene_annot = self.annotation_df["whole"][self.annotation_df["whole"]["gene_id"] == gene]
            transcripts = sorted(set(gene_annot["transcript_id"]))
            start_pos, end_pos = self.get_continous_coverage_ends(list(gene_annot["Start"]), list(gene_annot["End"]))
            segment_lengths = [end_pos[i] - start_pos[i] for i in range(len(start_pos))]
            
            total_segment_length = sum(segment_lengths)
            num_segments = len(start_pos)
            
            exon_height = exon_panel_height/len(transcripts)*0.6
            transcript_row_gap = exon_panel_height/(len(transcripts)+1)

            intersegment_gap = 40
            draw_length = transcript_line_width - (num_segments+1)*intersegment_gap
            
            tab["Gene"] = gene

            basic_scaling = (coverage_panel_height)/(max_expect_covs*exon_height)

            for sample in self.bam_coverage.keys():

                scale_coef = max([self.cov_stats[sample][gene][s]["scale"] for s in self.cov_stats[sample][gene].keys() if "CoD" in self.cov_stats[sample][gene][s].keys()])  
                
                output_path = os.path.join(self.output_path, "coverage/coverage_plots/"+sample+"/"+gene+"/")
                
                d = CairoDrawer(output_path+"coverage.png", width = page_width, height = page_height)

                ## CREATE HEADER
                d.draw_text(text="SIRVsuite coverage", y = title_y, x = page_width/2, h_align="center", font_size=title_size)
                
                tab["Sample"] = sample
                tab["Mode"] = mode
                
                exon_panel_y = header_y + header_height + panel_gap_y
                transcript_line_offset_y = exon_panel_y + exon_panel_height * 0.05
                transcript_text_y = transcript_line_offset_y
                coverage_panel_y = exon_panel_y + exon_panel_height + panel_gap_y
                
                d.draw_table(table_dict=tab,x=exon_panel_x,y=header_y,width=exon_panel_width/3,height=header_height)

                d.draw_rectangle(x = exon_panel_x,
                    y = exon_panel_y,
                    width = exon_panel_width,
                    height = exon_panel_height,
                    color_fill=(.5,.5,.5),
                    alpha = 0.3)

                t_y = exon_panel_y + transcript_row_gap*2/3 

                for transcript in transcripts:

                    transcript_annot = gene_annot[gene_annot["transcript_id"] == transcript]

                    if "+" in set(transcript_annot["Strand"]):
                        color_exon = (60/255,140/255,80/255)
                    elif "-" in set(transcript_annot["Strand"]):
                        color_exon = (0/255,120/255,180/255)
                    
                    d.draw_line(x = transcript_line_offset_x, 
                        y = t_y + exon_height/2, 
                        width = transcript_line_width, line_width = 2)

                    segment_start = transcript_line_offset_x + intersegment_gap
                    for segment_idx in range(len(start_pos)):
                        exons_in_segment = transcript_annot[(transcript_annot["Start"] >= start_pos[segment_idx]) & (transcript_annot["End"] <= end_pos[segment_idx])]
                        lengths = exons_in_segment["End"] - exons_in_segment["Start"] + 1
                    
                        for exon_idx in range(len(exons_in_segment)):
                            relative_pos_x = segment_start + (int(exons_in_segment.iloc[exon_idx]["Start"]) - start_pos[segment_idx] + 1) / total_segment_length * draw_length
                            relative_width = (int(exons_in_segment.iloc[exon_idx]["End"]) - int(exons_in_segment.iloc[exon_idx]["Start"]) + 1) / total_segment_length * draw_length

                            d.draw_rectangle(x=relative_pos_x, y=t_y, width=relative_width, height=exon_height, alpha=1, color_fill=color_exon)

                        segment_start += segment_lengths[segment_idx]/total_segment_length*draw_length + intersegment_gap

                    t_y += transcript_row_gap

                    d.draw_text(text = transcript, 
                        x = transcript_text_x, 
                        y = t_y - 34, 
                        font_size = 28,
                        v_align="bottom")
                
                for strand in self.expected_coverage["whole"][gene].keys():
                    
                    segment_start = transcript_line_offset_x + intersegment_gap

                    for segment_idx in range(len(start_pos)):

                        gene_pos = self.gene_coords[gene]

                        start = start_pos[segment_idx] - gene_pos[0]
                        end = end_pos[segment_idx] - gene_pos[0]

                        expected_cov = self.expected_coverage["whole"][gene][strand]
                        max_e = max(expected_cov)
                        expected_cov = expected_cov[start:end]

                        if (strand == "+"):
                            upside_down = True
                            color_fill = (60/255,140/255,80/255)
                        elif (strand == "-"):
                            upside_down = False
                            color_fill = (0/255,120/255,180/255)

                        real_cov_scaled = self.bam_coverage[sample][gene][strand][start:end] / scale_coef 
                        
                        expected_coverage_height_total = coverage_panel_height/2 * max_e / max_expect_covs

                        d.draw_signal(signal = expected_cov, x = segment_start, y = coverage_panel_y + coverage_panel_height/2,
                        width = segment_lengths[segment_idx]/total_segment_length*draw_length, height = expected_coverage_height_total,
                        mode = "segment", upside_down = upside_down, line_width = 1, color_fill = color_fill, alpha_fill = 0.5, alpha_line = 0.8, y_max = max_expect_covs)

                        d.draw_signal(signal = real_cov_scaled, x = segment_start, y = coverage_panel_y + coverage_panel_height/2,
                        width = segment_lengths[segment_idx]/total_segment_length*draw_length, height = expected_coverage_height_total,
                        mode = "normal", upside_down = upside_down, 
                        y_max=max_expect_covs,
                        line_width = 1, color_fill = color_fill, alpha_fill = 0.5, alpha_line = 0.8)
                         
                        segment_start += segment_lengths[segment_idx]/total_segment_length*draw_length + intersegment_gap
                d.finish()
