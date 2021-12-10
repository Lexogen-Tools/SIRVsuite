
import numpy as np
import os 
import re
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
from ..helper import *
from ..countReader import countReader
import logging

log = logging.getLogger(__name__.split(".")[-1])

## This module relevant to full-trancriptome libraries, i.e. CORALL

class SIRVsuiteConcentration():

    def __init__(self, sample_sheet = None, output_dir = "./", experiment_name = ""):
        """
        """

        self.output_dir = output_dir
        self.experiment_name = experiment_name

        if (sample_sheet == None):
            raise ValueError("Please load and pass sample sheet..")
        
        grouped_samples = dict()
        for sample in sample_sheet.keys():
            if "replicate_group" in sample_sheet[sample].keys():
                replicate_group = sample_sheet[sample]["replicate_group"]
                if replicate_group == "none":
                    replicate_group = sample
            else:
                replicate_group = sample
            
            if not replicate_group in grouped_samples.keys():
                grouped_samples[replicate_group] = dict()
                grouped_samples[replicate_group]["path"] = [sample_sheet[sample]["counting_path"]]
                method = "mix2"
            else:
                grouped_samples[replicate_group]["path"].append(sample_sheet[sample]["counting_path"])
        
        self.grouped_samples = grouped_samples

        cnts = dict()
        for group in grouped_samples.keys():
            c = countReader()
            cnts[group] = c.read_counting_file(files = grouped_samples[group]["path"], spike_in_type=["SIRV","ERCC"], counting_method=method)
        
        self.cnts = cnts
        self.data = self.get_relative_abundance(cnts)
        self.export_data(self.data)

    def heatmap_generator(self, data, row_labels, col_labels, title = "SIRVsuite heatmap", ax=None,
                cbar_kw={}, cbarlabel="", no_x_ticks = True, **kwargs):
        """
        This method creates a heatmap with a colorbar 
        """

        if not ax:
            ax = plt.gca()

        colorsList = [(0,0,1),(0,1,0),(1,0,0)]
        CustomCmap = matplotlib.colors.LinearSegmentedColormap.from_list("my_cmap",colorsList,N=256,gamma=1)
        
        im = ax.imshow(data, **kwargs, vmin = -1, vmax = 1, cmap = CustomCmap)
        ax.yaxis.set_label_coords(-0.1,10)
        
        fraction = 1

        cbar = ax.figure.colorbar(im, ax = ax, **cbar_kw, shrink = fraction*(6/len(row_labels)), pad = 0.05)
        cbar.ax.set_ylabel(cbarlabel, rotation = -90, va = "top", fontsize = 16, labelpad = 20)

        if not no_x_ticks:
            ax.set_xticks(np.arange(data.shape[1]))
            ax.set_xticklabels(col_labels)
        else:
            ax.set_xticks([])
            ax.set_xticklabels([])

        ax.set_yticks(np.arange(data.shape[0]))
        ax.set_yticklabels(row_labels)

        ax.tick_params(top=True, bottom=False,
                    labeltop=True, labelbottom=False, labelsize = 16)

        plt.setp(ax.get_xticklabels(), rotation=-90, ha="center",
                rotation_mode="default")

        for _, spine in ax.spines.items():
            spine.set_visible(False)

        ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
        ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
        ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
        ax.tick_params(which="minor", bottom=False, left=False)

        return im, cbar


    def annotate_heatmap(self, image, data=None, valfmt="{x:.2f}",
                        textcolors=("black", "white"),
                        threshold=None, **textkw):
        """
        Create annotation for heatmap to compare relative abundance values
        """

        if not isinstance(data, (list, np.ndarray)):
            data = image.get_array()

        if threshold is not None:
            threshold = image.norm(threshold)
        else:
            threshold = image.norm(data.max())/2.

        kw = dict(horizontalalignment="center",
                verticalalignment="center")
        kw.update(textkw)

        if isinstance(valfmt, str):
            valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

        texts = []
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                kw.update(color=textcolors[int(image.norm(data[i, j]) > threshold)])
                text = image.axes.text(j, i, valfmt(data[i, j], None), **kw)
                texts.append(text)

        return texts

    def __get_text_size__(self, fontsize,text):
        """
        
        """
        ffam = 'DejaVu Sans'
        fp = matplotlib.font_manager.FontProperties(
        family=ffam, style='normal', size=fontsize,
        weight='normal', stretch='normal')
        pth = matplotlib.textpath.TextPath((100,100), text, prop=fp)
        bb = pth.get_extents()

        return bb.width, bb.height

    def get_relative_abundance(self, cnts):
        """
        """
        norm_abund_dict = {}

        for group in cnts.keys():
            averaged_SIRV_abundance = {t:np.mean(cnts[group]["SIRV"][t]) for t in sorted(cnts[group]["SIRV"].keys())}
            ## NORMALIZATION OF TRANSCRIPT FPKM ##
            expected_quantity = np.mean(list(averaged_SIRV_abundance.values()))

            if (expected_quantity == 0):
                {t: "NaN" for t in sorted(averaged_SIRV_abundance.keys())}
            else:
                norm_abund_dict[group] = {t:( averaged_SIRV_abundance[t]/expected_quantity ) for t in sorted(averaged_SIRV_abundance.keys())}
        
        return (norm_abund_dict)

    def export_data(self, norm_abund_dict):
        """
        """
        groups = list(norm_abund_dict.keys())

        path = os.path.join(self.output_dir,"concentration/")

        if (not os.path.exists(path)):
            os.makedirs(path)

        with open(os.path.join(path,"relative_concentration.tsv"), "w") as out_file:
            out_file.write("transcript_id\t"+"\t".join(groups)+"\n")
            transcript_ids = sorted(norm_abund_dict[groups[0]].keys())
            for transcript in transcript_ids:
                out_file.write(transcript+"\t"+"\t".join([str(norm_abund_dict[group][transcript]) for group in groups])+"\n")

    def create_sirvsuite_boxplot(self, relative_abundance):
        """
        """
        log.info("creating SIRV E0 concentration boxplot")

        path = os.path.join(self.output_dir,"concentration/")

        if not os.path.exists(path):
            os.makedirs(path)
        
        fig, ax1 = plt.subplots(figsize=(12, 2 + len(relative_abundance)))
        ax1.set_xscale('log')

        if self.experiment_name != "":
            self.experiment_name += ": "

        output_name = os.path.join(path, "SIRVsuite_boxplot.png")
        ax1.set_title(self.experiment_name+"SIRV E0 relative concentration distribution", size = 16, loc = 'center',
                    fontdict = {'verticalalignment':'baseline'})

        heatmap_matrix = self.__count_dict_to_matrix__(relative_abundance)

        for i in range(np.shape(heatmap_matrix)[1]):
            relative_conc = heatmap_matrix[:,i]
            x = np.ones(len(relative_conc))*(i+1)
            ax1.scatter(relative_conc, x, alpha = 0.5, s = 10, color = 'gray')

        ax1.boxplot(heatmap_matrix,
                    0,
                    "rs",
                    0, 
                    widths = 0.5, 
                    boxprops = dict(facecolor='green', color='white', alpha = 0.7),
                    medianprops = dict(color='white', linewidth = 2),
                    whiskerprops = dict(linestyle = (0,(5,5)), alpha = 0.4, color = 'black'),
                    meanprops = dict(marker='o',markeredgecolor='red', markerfacecolor = 'none', markersize = 8 , alpha = 0.8),
                    capprops = dict(alpha = 0.4, color = 'black'),
                    showmeans = True,
                    patch_artist = True,
                    showfliers = False)
        
        ax1.plot([1,1],[0,len(relative_abundance) + 1], color = 'darkblue', alpha = .6)
        ax1.set_ylim([0,len(relative_abundance)+1])


        limit_x = heatmap_matrix.max()

        ax1.set_xlim([np.min(relative_conc)-np.min(relative_conc)*0.8,limit_x])
        
        ax1.set_yticklabels(relative_abundance.keys())
        plt.xlabel("relative SIRV transcript concentration")
        fig.tight_layout()
        fig.savefig(output_name, dpi = 300, pad_inches = 0.2)
    
    def __count_dict_to_matrix__(self, count_dict):
        """
        """
        conc_matrix = np.array([])
        groups = count_dict.keys()
        for group in groups:
            transcripts = sorted(count_dict[group].keys())
            conc_matrix = np.append(conc_matrix, [count_dict[group][t] for t in transcripts])

        return np.reshape(conc_matrix, (len(transcripts),len(groups)), order = 'F')
    

    def create_sirvsuite_heatmap(self, relative_abundance):
        """
        """
        log.info("creating SIRV E0 concentration heatmap")

        path = os.path.join(self.output_dir,"concentration/")

        if not os.path.exists(path):
            os.makedirs(path)

        groups = list(relative_abundance.keys())
        transcript_names = np.array((sorted(relative_abundance[groups[0]].keys())))

        genes = ["SIRV1","SIRV2","SIRV3","SIRV4","SIRV5","SIRV6","SIRV7"]

        fig = plt.figure(figsize=(2+len(relative_abundance),40))

        vectorized_len = np.vectorize(len)
        max_length_index =  np.argmax(vectorized_len(groups))
        text_width, text_height = self.__get_text_size__(16, groups[max_length_index])

        if self.experiment_name != "":
            self.experiment_name += ": "

        plt.title(self.experiment_name+"SIRV E0 relative concentrations", pad = text_width + 10 + 25, fontdict = {'fontsize': 24})
        plt.axis('off')

        transcript_num = np.zeros(len(genes)).astype(int)
        for i in range(len(genes)):
            for transcript in transcript_names:
                if re.match(genes[i], transcript):
                    transcript_num[i] += 1

        spec = gridspec.GridSpec(ncols=1, nrows = len(genes), height_ratios = transcript_num)

        for idx in range(len(genes)):

            filtered = [ i for i in range(len(transcript_names)) if bool(re.match(genes[idx],transcript_names[i])) ]

            heatmap_matrix = self.__count_dict_to_matrix__(relative_abundance) ## creating matrix for heatmap
            heatmap_matrix = heatmap_matrix[filtered,:]
            t_filtered = transcript_names[filtered]

            LFC = np.log2(heatmap_matrix)

            ax = fig.add_subplot(spec[idx])

            text_width, text_height = self.__get_text_size__(16, genes[idx])

            ax.text(-3.5, float(transcript_num[idx])/2 - 0.5," "*10+genes[idx]+" "*10,
            rotation = 90,
            verticalalignment = "center",
            horizontalalignment = "center",
            fontsize = 16, 
            bbox = {'facecolor':'silver'})

            if (idx == 0):
                draw_x_ticks = False
            else:
                draw_x_ticks = True

            im, cbar = self.heatmap_generator(LFC, t_filtered, groups, ax=ax, cbarlabel="Log2 Fold Change ($c_r$)", no_x_ticks=draw_x_ticks)
            texts = self.annotate_heatmap(im, data = LFC, valfmt="{x:.2f}", threshold=-10)

        plt.margins(1,1)
        fig.tight_layout()

        output_name = os.path.join(path, "SIRVsuite_heatmap.png") 
        fig.savefig(output_name, dpi = 300, pad_inches = 0.2, bbox_inches = 'tight')