import pysam as ps
import numpy as np
import pyBigWig as bwig
import gtfparse
import pandas as pd
from scipy.stats import norm
import os
import copy
import warnings
from ..helper import path_features
from . import cairoDrawer as draw
import logging


log = logging.getLogger(__name__.split(".")[-1])


class SIRVsuiteCoverage():
    """
    SIRVsuiteCoverage class for creating coverage plots & statistics

    The class is concepted to be used in following steps:
        1. load annotation file, which contains**
        2. calculate expected coverage based on the annotation
        3. calculate coverage from an alignment file
        4. calculate statistics
        5. create coverage plots

    ** by default, the annotation is loaded from a Resources folder
    in SIRVsuite project (annotation is identical to the one which can
    be found https://www.lexogen.com/sirvs/download/ under section SIRV-Set 4)

    """

    def __init__(self, sample_sheet=None, output_dir="./", experiment_name=""):
        """
        """

        self.output_path = output_dir

        log.info("calculating spike-in gene coverage")
        log.debug("SIRVsuite coverage creator initialized")
        
        self.annotation_path = "/".join(__file__.split("/")[:-3]) + "/Resources/SIRV_isoforms_multi-fasta-annotation_190418a.gtf"
        self.load_annotation(self.annotation_path)

        self.experiment_name = experiment_name

        if sample_sheet is None:
            raise ValueError("please specify path to the sample sheet..")
        else:
            self.sample_dict = sample_sheet

        # constants
        self._strands = ["-", "+"]

    def __values2intervals__(self, data, start_pos=0):
        """
        A method to convert data array (coverage) into intervals specified by starts, ends and values
        """
        tmp = np.logical_and(data[1:] != data[:-1], np.logical_or(np.isnan(data[1:]) == False, np.isnan(data[:-1]) == False))
        starts = np.append(np.array([start_pos]), np.where(tmp)[0] + start_pos + 1)
        ends = np.append(starts[1:], start_pos + len(data))
        values = data[starts - start_pos]
        
        return (starts, ends, values)

    def __strand2text__(self, text):
        """
        A helper method to convert
        """
        if (text == "+"):
            return "sense"
        elif (text == "-"):
            return "antisense"

    def __cov_data_export__(self, output_path, output_type="bigwig"):
        """

        """

        log.debug("Exporting coverage data")

        # function used for exporting coverage data to "bigwig" or "cov" file
        if (hasattr(self, "bam_coverage")):
            cov_dict = self.bam_coverage

            for strand in self._strands:
                for sample in cov_dict.keys():

                    if (output_type.lower() == "bigwig"):
                        with bwig.open(os.path.join(output_path, sample+"_"+self.__strand2text__(strand)+".bw"), "w") as bw:
                            header = [(gene, len(cov_dict[sample][gene][strand])) for gene in sorted(cov_dict[sample].keys()) if strand in cov_dict[sample][gene].keys()]
                            if (len(header) == 0):
                                continue

                            bw.addHeader(header)
                            for gene in sorted(cov_dict[sample].keys()):
                                (starts, ends, values) = self.__values2intervals__(cov_dict[sample][gene][strand])
                                chroms = np.array([gene] * len(values))

                                # if trimmed, need to use + self.gene_coords[gene][0] to start & end
                                bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=values.tolist())

                    elif (output_type.lower() == "cov"):
                        for gene in sorted(cov_dict[sample].keys()):
                            with open(os.path.join(output_path, sample+"_"+gene+"_"+self.__strand2text__(strand)), "w") as cov:
                                cov.write("position"+"\t"+"coverage"+"\n")
                                coverage = cov_dict[sample][gene][strand]
                                for idx in range(len(coverage)):
                                    cov.write(str(idx+1)+"\t"+str(coverage[idx])+"\n")
        else:
            warnings.warn("No coverage detected to export..")

    def bam_to_coverage(self, output_type="bigwig"):
        """
        A method for loading coverage from bam files
        """

        sample_dict = self.sample_dict

        log.debug("Calculating real coverage from input files")

        stats = {"num_reads": 0}

        bam_coverage = dict()
        stat_dict = dict()

        index_created = False

        for sample in sample_dict.keys():

            strandeness = sample_dict[sample]["read_orientation"]
            bam_path = sample_dict[sample]["alignment_path"]

            if (not os.path.exists(bam_path)):
                log.debug(bam_path+"does not exist.. skipping to the next file")
            else:
                index_path = bam_path + ".bai"

                if (not os.path.exists(index_path)):
                    ps.index(bam_path)
                    index_created = True
                    log.info(".bai index file not detected.. creating index file..")

            bamFile = ps.AlignmentFile(bam_path, "rb")

            bam_coverage[sample] = dict()
            stat_dict[sample] = dict()

            annotation_df = self.annotation_df["whole"]

            for gene in self.target_gene_id:

                if (hasattr(self, "annotation_df")):
                    contig_id = annotation_df["seqname"][annotation_df["gene_id"] == gene].values[0]
                else:
                    contig_id = gene

                if contig_id not in bamFile.references:
                    #del self._expected_coverage["whole"][contig_id]
                    #del self.cov_stats[contig_id]
                    continue

                bam_coverage[sample][gene] = dict()
                stat_dict[sample][gene] = dict()

                # check if annotation processed
                if (hasattr(self, "gene_coords")):
                    """
                    if (False):
                        bam_gene_start = self.gene_coords[gene][0]
                        bam_gene_end = self.gene_coords[gene][1]
                        length = bam_gene_end - bam_gene_start + 1
                        length = [bamFile.lengths[i] for i in range(len(bamFile.references)) if bamFile.references[i] == gene][0]
                    """
                    bam_gene_start = 1
                    bam_gene_end = self.gene_coords[gene][1]
                    length = bam_gene_end - bam_gene_start + 1
                else:
                    raise ValueError("expected coverage attribute not detected")

                if strandeness in ["rev", "fwd"]:
                    bam_coverage[sample][gene]["+"] = np.zeros(length)
                    bam_coverage[sample][gene]["-"] = np.zeros(length)
                    stat_dict[sample][gene]["+"] = copy.deepcopy(stats)
                    stat_dict[sample][gene]["-"] = copy.deepcopy(stats)

                for fragmentRead in bamFile.fetch(contig_id, bam_gene_start, bam_gene_end):

                    if(strandeness in ["rev", "fwd"]):
                        # inferring strandeness of a read
                        if ((fragmentRead.is_paired and fragmentRead.is_read1 and not fragmentRead.is_reverse) or
                            (fragmentRead.is_paired and fragmentRead.is_read2 and fragmentRead.is_reverse) or
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
                    positions = positions[(positions < bam_gene_end) & (positions > bam_gene_start)] - bam_gene_start
                    bam_coverage[sample][gene][strand][positions] += 1

            if index_created:
                os.remove(index_path)

        self.cov_stats = stat_dict
        self.bam_coverage = bam_coverage

        # export coverage to .bw format
        cov_track_path = os.path.join(self.output_path, "coverage", output_type)
        if (not os.path.exists(cov_track_path)):
            os.makedirs(cov_track_path)

        self.__cov_data_export__(cov_track_path, output_type=output_type)

    def load_annotation(self, annotation_path):
        """
        A method for loading annotation using gtfparse library
        """
        annotation_type = path_features(annotation_path)["extension"]

        try:
            log.debug("Loading "+annotation_path)

            if (annotation_type == "gtf"):
                # disable root log due to gtfparse logging interference
                log_root = logging.getLogger("")
                log_root.disabled = True
                # read gtf
                annotation_df = gtfparse.read_gtf(annotation_path)
                # enable root log
                log_root.disabled = False

            # limit only to a feature of interest
            annotation_df[annotation_df["feature"] == "exon"]
            # trim to required cols only
            annotation_df = annotation_df[["seqname", "start", "end", "strand", "gene_id", "transcript_id"]]

        except Exception as e:
            raise ValueError("An error occured while loading annotation: "+str(e))

        annotation_df["start"] -= 1
        self.annotation_df = {"whole": annotation_df}

    def get_continous_coverage_ends(self, starts, ends):
        """
        The method returnes continous segment of oveerlapping features (exons) from list of starts and ends of the features
        """

        sorted_idx = np.argsort(starts)

        contig_start = [starts[sorted_idx[0]]]
        contig_end = []
        thre = ends[sorted_idx[0]]

        for i in range(0, len(sorted_idx)-1):

            prev_idx = ends[sorted_idx[i]]
            next_idx = starts[sorted_idx[i+1]]

            if prev_idx > thre:
                thre = prev_idx

            if (next_idx > thre):
                contig_end.append(thre)
                contig_start.append(next_idx)

        contig_end.append(thre)

        return contig_start, contig_end

    def calc_start_transition(self, transcriptLength, mean=200, std=80):
        """
        The method estimates probability distribution of a fragment according to its length, mean and std of the model distribution
        """

        lengthProbability = norm.pdf(np.arange(1, transcriptLength+1), mean, std)

        lengthProbability /= np.sum(lengthProbability)

        fragmentProbability = np.zeros((transcriptLength, transcriptLength))

        for length in np.arange(1, transcriptLength):
            fragmentProbability[length, 0:transcriptLength-length+1] = lengthProbability[length]/(transcriptLength-length+1)

        startDistance = np.sum(fragmentProbability, axis=0)

        adjustedLength = np.sum(lengthProbability * (transcriptLength - np.arange(1, transcriptLength+1) + 1))

        startDistance /= np.max(startDistance)

        return startDistance, lengthProbability, adjustedLength

    def calc_statistics(self):
        """
        The method checks for bam coverages and expected coverages and then calculates CoD and scaling factor for particular samples, genes and strands
        """

        log.debug("Calculating statistics")

        if (not hasattr(self, "expected_coverage")):
            raise ValueError("Expected coverage not detected")

        # check if all attributes present
        if (not hasattr(self, "cov_stats") or not hasattr(self, "bam_coverage")):
            raise ValueError("Measured coverage not detected")

        CoD_table_path = os.path.join(self.output_path, "coverage/")
        if not os.path.exists(CoD_table_path):
            os.makedirs(CoD_table_path)

        experiment_name = self.experiment_name
        if experiment_name != "":
            experiment_name = "_" + experiment_name

        for mode in self.annotation_df.keys():
            with open(CoD_table_path+"CoD_table"+experiment_name+"_"+mode+".tsv", "w") as CoD_table:
                header = "gene_id\t" + "\t".join(list(sorted(self.bam_coverage.keys()))) + "\n"
                CoD_table.write(header)

                rows = dict()

                for sample in sorted(self.bam_coverage.keys()):
                    for gene in self.bam_coverage[sample].keys():
                        for strand in self._expected_coverage[mode][gene].keys():
                            expected_cov = self._expected_coverage[mode][gene][strand]
                            measured_cov = self.bam_coverage[sample][gene][strand]
                            
                            # calculate CoD and scaling factor and assign them to the cov_stats dictionary
                            self.cov_stats[sample][gene][strand]["CoD"], self.cov_stats[sample][gene][strand]["scale"] = self.CoD(measured_cov, expected_cov)

                            row_key = gene+"("+strand+")"
                            if row_key not in rows.keys():
                                rows[row_key] = str(self.cov_stats[sample][gene][strand]["CoD"])
                            else:
                                rows[row_key] += "\t" + str(self.cov_stats[sample][gene][strand]["CoD"])

                for row_idx in rows.keys():
                    CoD_table.write(row_idx + "\t" + rows[row_idx]+"\n")

    def CoD(self, real_cov, expect_cov):
        """
        The method calculate CoD from 2 vectors (real and expected coverage) and returns coefficient for scaling and CoD

        """

        if (len(real_cov) != len(expect_cov)):
            raise ValueError("real and expected coverage must have the same length")

        sumBAC_theory = np.sum(expect_cov)  # sumBAC (sum of base annotation count), sum of expected coverage

        sumCoverage_real = np.sum(real_cov)  # sum of real coverage

        scaling_coefficient = sumCoverage_real / sumBAC_theory  # calculation of scaling factor for real coverage (for comparison)

        # Correction of scale coef
        if (scaling_coefficient < 1):
            scaling_coefficient = 1

        # Scaling
        coverageReal_scaled = real_cov / scaling_coefficient

        # COD computation
        difference = expect_cov - coverageReal_scaled
        COD = np.sum((difference)**2) / sumBAC_theory
        COD = np.round(COD, 4)

        return COD, scaling_coefficient

    def annot_2_UTR(self, region_length=200):
        """
        This method modifies annotation data frame and trims the transcripts to a length close to 3' (where QuantSeq is mostly expressed)

        """

        # check if annotation loaded
        if (not hasattr(self, "annotation_df")):
            warnings.warn("There was no annotation detected.. skipping UTR mode..")
            return
        else:
            annotation_df = self.annotation_df["whole"]

        UTR_annotation = annotation_df[0:0]  # Creating an empty data frame with the same structure as annotation data frame

        for strand in self._strands:
            transcripts = np.unique(annotation_df[annotation_df["strand"] == strand]["transcript_id"])
            for transcript in transcripts:
                annot_transcript = annotation_df[annotation_df["transcript_id"] == transcript]

                # initiliaze parameters for a particular strand
                if (strand == '+'):
                    ascending = False
                    const = -1
                    UTR_df_col = "start"
                elif (strand == '-'):
                    ascending = True
                    const = 1
                    UTR_df_col = "end"

                # sort annotation dataframe by Start column
                annot_transcript = annot_transcript.sort_values(by=["start"], axis=0, ascending=ascending).reset_index(drop=True)
                exon_lengths = annot_transcript["end"] - annot_transcript["start"]
                remainder = region_length
                UTR_data_frame = annot_transcript

                for index, length in enumerate(exon_lengths):
                    remainder = remainder - length
                    if (remainder < 0):
                        UTR_data_frame.loc[index, UTR_df_col] = UTR_data_frame.loc[index, UTR_df_col] + remainder*const
                        UTR_data_frame = UTR_data_frame.loc[0:index, :]
                        break

                UTR_annotation = UTR_annotation.append(UTR_data_frame)

        UTR_annotation = UTR_annotation.reset_index(drop=True)

        self.annotation_df["UTR_test"] = UTR_annotation

    def expected_coverage(self, transition_lengths=False, tlen_model_param=(0, 0), mode="normal", region_length=200):
        """
        Method to create expected coverage based on provided annotation.

        Experimental feature: transition lengths modelling
            input: list, ndarray or tuple of a length 2, specifying model parameters in the following order [mean, std], e.g. [100, 50]

        Returns: 3-level dictionary gene_id -> strand -> expected_coverage
        """

        log.debug("Calculating expected coverage from the annotation")

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

            self.target_gene_id = np.unique(annotation_df["gene_id"])

            expected_coverage[mode] = dict()

            for gene in self.target_gene_id:
                gene_annot = annotation_df[annotation_df["gene_id"] == gene]

                if (len(gene_annot) == 0):
                    warnings.warn("No record for gene "+gene+" has been found in the provided annotation file.. skipping gene..")

                end_coord = np.max(gene_annot["end"])
                start_coord = 1

                self.gene_coords[gene] = (start_coord, end_coord)

                expected_coverage[mode][gene] = dict()

                for strand in self._strands:

                    stranded_gene_annot = gene_annot[gene_annot["strand"] == strand]

                    # skip when no record for gene and strand
                    if (len(stranded_gene_annot) == 0):
                        continue

                    # initialize base for coverage with 0s
                    expected_coverage[mode][gene][strand] = np.zeros(end_coord-start_coord+1)

                    for transcript in np.unique(stranded_gene_annot["transcript_id"]):

                        transcript_annot = stranded_gene_annot[stranded_gene_annot["transcript_id"] == transcript].copy(deep=True).reset_index(drop=True)
                        transcript_annot["start"] += 1

                        exon_lengths = np.array(transcript_annot["end"]) - np.array(transcript_annot["start"]) + 1

                        # if parameters for transition transcript lengths are defined, use probabilistic model to calculate multiplicative weights for starts and ends,
                        # otherwise set weights to 1
                        if (transition_lengths):
                            if (tlen_model_param == (0, 0)):
                                tlen_model_param = (25, 30)

                            weights, _, _ = self.calc_start_transition(np.sum(exon_lengths), mean=tlen_model_param[0], std=tlen_model_param[1])
                        else:
                            weights = np.ones(np.sum(exon_lengths))

                        weight_idx = 0

                        for row_idx, row in transcript_annot.iterrows():

                            shifted_start = row["start"] - start_coord
                            shifted_end = row["end"] - start_coord + 1

                            # apply weights for transition lengths for transcript start and end
                            weights = weights * weights[::-1]
                            expected_coverage[mode][gene][strand][shifted_start:shifted_end] += weights[weight_idx:(weight_idx + exon_lengths[row_idx])]

                            weight_idx += exon_lengths[row_idx]

        self._expected_coverage = expected_coverage

    def coverage_plot(self):

        """
        This method utilizes CairoDrawer class to draw objects onto surface of
        interest. In current version, .png format of the surface is hardcoded.
        Here, all dimensions, coordinates, colors and draw elements within 
        the surface are defined. It creates the graphical output separately 
        for each gene and sample.
        """

        log.debug("Creating coverage plots")

        # define base sizes, coordinates
        page_width = 2000
        page_height = 1324

        title_y = 40
        title_size = 40

        panel_gap_y = title_size

        header_y = title_size + panel_gap_y
        header_height = page_height/6

        exon_panel_x = page_width/28
        exon_panel_width = page_width - 2*exon_panel_x
        exon_panel_height = page_height/3

        transcript_line_offset_x = exon_panel_x + exon_panel_width*0.1

        transcript_line_width = exon_panel_width*0.8

        transcript_text_x = exon_panel_x + (transcript_line_offset_x - exon_panel_x)/2

        coverage_panel_height_orig = page_height/3

        if "UTR" in self.annotation_df.keys():
            mode = "UTR"
        else:
            mode = "whole"

        if self.experiment_name != "":
            tab = {"Experiment": self.experiment_name}
            page_height = 1400
        else:
            tab = {}
            header_height = header_height*3/4

        for gene in self._expected_coverage["whole"].keys():

            max_expect_covs = max([max(self._expected_coverage["whole"][gene][s])
                for s in self._expected_coverage["whole"][gene].keys()])

            gene_annot = self.annotation_df["whole"][self.annotation_df["whole"]["gene_id"] == gene]
            transcripts = sorted(set(gene_annot["transcript_id"]))
            start_pos, end_pos = self.get_continous_coverage_ends(list(gene_annot["start"]), list(gene_annot["end"]))
            segment_lengths = [end_pos[i] - start_pos[i] for i in range(len(start_pos))]

            total_segment_length = sum(segment_lengths)
            num_segments = len(start_pos)

            coverage_panel_height = coverage_panel_height_orig
            exon_panel_y = header_y + header_height + panel_gap_y

            existing_strands = list(self._expected_coverage["whole"][gene].keys())

            if (len(transcripts) < 3):
                exon_panel_height = 60*len(transcripts)
                exon_height = exon_panel_height / len(transcripts)*0.6
                transcript_row_gap = exon_panel_height / (len(transcripts)+1) * (1-0.6)
            else:
                real_cov_limit_factor = 1
                exon_panel_height = page_height/3
                exon_height = exon_panel_height / len(transcripts)*0.6
                transcript_row_gap = exon_panel_height / (len(transcripts)+1) * (1-0.6)

            
            coverage_panel_y = exon_panel_y + exon_panel_height + panel_gap_y
            coverage_panel_upper_y = coverage_panel_y 

            x_axis_y = coverage_panel_y + coverage_panel_height + panel_gap_y
            # fix, if only one strand, scale
            if len(existing_strands) == 1:
                if (existing_strands[0] == "+"):
                    coverage_panel_y += panel_gap_y/2 + coverage_panel_height + coverage_panel_height_orig/3
                    x_axis_y = coverage_panel_y + coverage_panel_height/2 + panel_gap_y
                elif (existing_strands[0] == "-"):
                    coverage_panel_y -= coverage_panel_height/2

            # depreciated factor
            # basic_scaling = (coverage_panel_height) / (max_expect_covs*exon_height)

            intersegment_gap = 40
            draw_length = transcript_line_width - (num_segments+1) * intersegment_gap

            tab["Gene"] = gene

            for sample in self.bam_coverage.keys():

                if gene not in self.cov_stats[sample].keys():
                    continue

                scale_coef = max([
                    self.cov_stats[sample][gene][s]["scale"] for s in self.cov_stats[sample][gene].keys()
                    if "CoD" in self.cov_stats[sample][gene][s].keys()
                    ])

                output_path = os.path.join(self.output_path, "coverage/coverage_plots/" + sample + "/" + gene + "/")

                d = draw.CairoDrawer(output_path + "coverage.png", width=page_width, height=page_height)

                d.draw_text(text="SIRVsuite coverage", y=title_y, x=page_width/2, h_align="center", font_size=title_size)

                tab["Sample"] = sample
                tab["Mode"] = mode

                d.draw_table(table_dict=tab, x=exon_panel_x, y=header_y, width=exon_panel_width/3, height=header_height)

                d.draw_rectangle(
                    x=exon_panel_x,
                    y=exon_panel_y,
                    width=exon_panel_width,
                    height=exon_panel_height,
                    color_fill=(.5, .5, .5),
                    alpha=0.15
                    )

                transcript_y = exon_panel_y + transcript_row_gap

                for transcript in transcripts:

                    transcript_annot = gene_annot[gene_annot["transcript_id"] == transcript]

                    if "+" in set(transcript_annot["strand"]):
                        color_exon = (60/255, 140/255, 80/255)
                    elif "-" in set(transcript_annot["strand"]):
                        color_exon = (0/255, 120/255, 180/255)

                    d.draw_line(
                        x=transcript_line_offset_x,
                        y=transcript_y + exon_height/2,
                        width=transcript_line_width, line_width=2
                        )

                    d.draw_text(
                        text=transcript,
                        x=transcript_text_x,
                        y=transcript_y + exon_height/2,
                        font_size=28,
                        v_align="center"
                        )

                    segment_start = transcript_line_offset_x + intersegment_gap
                    for segment_idx in range(len(start_pos)):

                        exons_in_segment = transcript_annot[(transcript_annot["start"] >= start_pos[segment_idx]) & (transcript_annot["end"] <= end_pos[segment_idx])]
                        lengths = list(exons_in_segment["end"] - exons_in_segment["start"] + 1)

                        for exon_idx in range(len(exons_in_segment)):
                            relative_pos_x = segment_start + (int(exons_in_segment.iloc[exon_idx]["start"]) - start_pos[segment_idx] + 1) / total_segment_length * draw_length
                            relative_width = lengths[exon_idx] / total_segment_length * draw_length
                            d.draw_rectangle(x=relative_pos_x, y=transcript_y, width=relative_width, height=exon_height, alpha=1, color_fill=color_exon)

                        segment_start += segment_lengths[segment_idx]/total_segment_length*draw_length + intersegment_gap

                    transcript_y += transcript_row_gap + exon_height

                stats_table = {}

                for strand in existing_strands:

                    stats_table["CoD("+strand+")"] = "%.4f" % (self.cov_stats[sample][gene][strand]["CoD"])
                    stats_table["total reads("+strand+")"] = "%d" % (self.cov_stats[sample][gene][strand]["num_reads"])

                    expected_cov = self._expected_coverage["whole"][gene][strand]
                    real_cov_scaled = self.bam_coverage[sample][gene][strand] / scale_coef
                    max_e = max(expected_cov)

                    if len(existing_strands) == 1:
                        real_cov_limit_factor = max([max(real_cov_scaled), max_e]) / max_e
                        if len(transcript) < 3:
                            max_cov_scale_limit = (coverage_panel_y - coverage_panel_upper_y)*2 / coverage_panel_height
                            real_cov_limit_factor = min([real_cov_limit_factor, max_cov_scale_limit])
                        else:
                            min([real_cov_limit_factor, 2.0])

                    segment_start = transcript_line_offset_x + intersegment_gap

                    axis_shift_y = 5

                    # draw coverage levels

                    if (strand == "+"):
                        upside_down = False
                        color_fill = (60/255, 140/255, 80/255)
                        rotate_axis = 270
                        factor = -1

                    elif (strand == "-"):
                        upside_down = True
                        color_fill = (0/255, 120/255, 180/255)
                        rotate_axis = 90
                        factor = 1

                    delta = coverage_panel_height / 2 / max_expect_covs*factor
                    for i in range(1, int(max_expect_covs)+1):
                        d.draw_line(
                            x=transcript_line_offset_x,
                            y=coverage_panel_y + coverage_panel_height / 2 + delta*i,
                            width=transcript_line_width, line_width=1, alpha=0.3
                            )
                        d.draw_line(
                            x=transcript_line_offset_x - 5,
                            y=coverage_panel_y + coverage_panel_height / 2 + delta*i,
                            width=10, color=(0.3, 0.3, 0.3)
                            )
                        d.draw_text(
                            text=str(int(i*scale_coef)),
                            x=transcript_line_offset_x - 10,
                            y=coverage_panel_y + coverage_panel_height/2 + delta*i,
                            h_align="right", v_align="center", font_size=16, color_rgb=(0.3, 0.3, 0.3)
                            )

                    d.draw_line(
                        x=transcript_line_offset_x,
                        y=coverage_panel_y + coverage_panel_height/2 + axis_shift_y*factor,
                        width=coverage_panel_height/2 * real_cov_limit_factor + panel_gap_y/2,
                        end_shape=("right", "arrow"), rotate=rotate_axis, color=(0.3, 0.3, 0.3)
                        )

                    d.draw_text(
                        text="reads" + strand, x=transcript_line_offset_x - panel_gap_y*2,
                        y=coverage_panel_y + coverage_panel_height/2 + factor*real_cov_limit_factor*coverage_panel_height/4,
                        rotate=90, font_size=24, rotation_point="center")
                    
                    d.draw_text(
                        text="( " + self.__strand2text__(strand) + " )",
                        x=transcript_line_offset_x - panel_gap_y*2 - 34,
                        y=coverage_panel_y + coverage_panel_height/2 + factor*real_cov_limit_factor*coverage_panel_height/4,
                        rotate=90, font_size=24, rotation_point="center")

                    for segment_idx in range(len(start_pos)):

                        gene_pos = self.gene_coords[gene]

                        start = start_pos[segment_idx] - gene_pos[0] + 1
                        end = end_pos[segment_idx] - gene_pos[0] + 1

                        expected_cov_slice = expected_cov[start:end]
                        real_cov_scaled_slice = real_cov_scaled[start:end]

                        expected_coverage_height_total = coverage_panel_height/2

                        d.draw_signal(
                            signal=expected_cov_slice, mode="segment", upside_down=upside_down,
                            x=segment_start, y=coverage_panel_y + coverage_panel_height/2,
                            width=segment_lengths[segment_idx]/total_segment_length*draw_length,
                            height=expected_coverage_height_total, y_max=max_expect_covs,
                            line_width=1, color_fill=color_fill, alpha_fill=0.7, alpha_line=0.9
                            )

                        d.draw_signal(
                            signal=real_cov_scaled_slice,
                            x=segment_start, y=coverage_panel_y + coverage_panel_height/2,
                            width=segment_lengths[segment_idx] / total_segment_length * draw_length,
                            height=expected_coverage_height_total * real_cov_limit_factor,
                            mode="normal", upside_down=upside_down,
                            y_max=max_expect_covs * real_cov_limit_factor, scale_factor=scale_coef,
                            color_fill=(0.2, 0.2, 0.2), color_line=(0.2, 0.2, 0.2),
                            alpha_fill=0.7, alpha_line=0.8, line_width=1
                            )

                        # draw coordinates
                        d.draw_line(
                            x=segment_start, y=x_axis_y,
                            width=segment_lengths[segment_idx]/total_segment_length*draw_length,
                            end_shape=("both", "line"), color=(0.3, 0.3, 0.3)
                            )
                        d.draw_text(
                            text=str(gene_pos[0] + start), x=segment_start, y=x_axis_y + 6,
                            font_size=16, v_align="center", h_align="left", rotate=90
                            )
                        d.draw_text(
                            text=str(gene_pos[0] + end - 1),
                            x=segment_start + segment_lengths[segment_idx]/total_segment_length*draw_length,
                            y=x_axis_y + 6, rotate=90, font_size=16, v_align="bottom", h_align="left")

                        segment_start += segment_lengths[segment_idx]/total_segment_length*draw_length + intersegment_gap

                d.draw_table(
                    table_dict=stats_table, x=page_width*2/3, y=header_y,
                    width=exon_panel_width/3, height=len(stats_table) * header_height/4
                    )

                d.finish()
