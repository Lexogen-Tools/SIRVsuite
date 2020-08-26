
import numpy as np
import os 
import re
import csv 
import argparse as ap
import matplotlib.pyplot as plt
import matplotlib
import pyranges
from matplotlib import gridspec

## This module relevant to full-trancriptome libraries, i.e. CORALL

def read_mix2(file_dict, quantity_unit = 'fpkm_chn'):

    # function used for loading mix-square output
    
    SIRV_gene_set3 = ["SIRV1","SIRV2","SIRV3","SIRV4","SIRV5","SIRV6","SIRV7"]
    aliases = file_dict["alias"]
    files = file_dict["path"] 

    data_dict = {}

    for idx in range(len(files)):
        if os.path.exists(files[idx]):
            with open(files[idx],'r') as MIX2_file:
                MIX2_reader = csv.reader(MIX2_file,delimiter = "\t")
                init = True
                for row in MIX2_reader:
                    row = np.array(row)
                    if init:
                        columns = np.array(["tracking_ID","gene_ID","abundance","comp_frags_locus","FPKM_CHN","FPKM_THN"]) # trim to necessary info
                        indices = (row[:,None] == columns).argmax(axis=0) # find which indexes in row corresponds to desired columns

                        init = False
                    
                    if row[1] in SIRV_gene_set3:

                        selected_columns = row[indices]

                        if quantity_unit.lower() == 'fpkm_chn':
                            value = float(selected_columns[4])  
                        elif quantity_unit.lower() == 'counts':
                            value = float(selected_columns[2] * selected_columns[3]) # multiplying abundance and comp_frag_locus result in counts for transcripts  
                        elif quantity_unit.lower() == 'abundance':
                            value = float(selected_columns[2])
                        elif quantity_unit.lower() == 'fpkm_thn':
                            value = float(selected_columns[5])
                        else:
                            print ("unknown quantity unit specified.. supported options are fpkm_chn or abundance")
                            break

                        if selected_columns[0] not in data_dict.keys():
                            data_dict[selected_columns[0]] = [value]
                        else:
                            data_dict[selected_columns[0]].append(value)
                    else:
                        print ("Could not find file: "+files[idx]+". Skipping to the next file...")
                        continue
                
    return data_dict

def read_cufflinks(file_dict):

    SIRV_gene_set3 = ["SIRV1","SIRV2","SIRV3","SIRV4","SIRV5","SIRV6","SIRV7"]
    aliases = file_dict["alias"]
    files = file_dict["path"] 

    data_dict = {}

    for idx in range(len(files)):
        if os.path.exists(files[idx]):
            data = pyranges.read_gtf(files[idx], as_df=True)
            
            # check if format correct
            required_cols = ["FPKM","transcript_id","gene_id"]
            if any([col not in data for col in required_cols]):
                raise ValueError("could not find FPKM feature in the gtf..")

            filter_SIRV = data[data["gene_id"].isin(SIRV_gene_set3)]

            if (len(filter_SIRV) == 0 ):
                raise ValueError("No record of SIRVs found..")
            elif ( (len(filter_SIRV) > 0 and len(filter_SIRV) < 68 ) or len(filter_SIRV) > 68):
                raise ValueError("Unexpected number of SIRVs transcripts..")
            
            for index, row in data.iterrows():
                data_dict[row["transcript_id"]] = row["FPKM"]

        else:
            print ("Could not find file: "+files[idx]+". Skipping to the next file...")
            continue
    return

def path_features(filePath):
    # function returning features of a given path
    matching_pattern = "((?:(?:.*\\/)*(.*)\\/))*(?:(.*)\\.(.*))*"
    feature_match = re.match(matching_pattern, filePath)
    out_dict = {"parent_dir": feature_match.group(2),
                "file_name": feature_match.group(3),
                "extension": feature_match.group(4),
                "path": feature_match.group(1)}
    return out_dict

def heatmap_generator(data, row_labels, col_labels, title = "SIRVsuite heatmap", ax=None,
            cbar_kw={}, cbarlabel="", no_x_ticks = True, **kwargs):

    # This method creates a heatmap with a colorbar 

    if not ax:
        ax = plt.gca()

    colorsList = [(0,0,1),(0,1,0),(1,0,0)]
    CustomCmap = matplotlib.colors.LinearSegmentedColormap.from_list("my_cmap",colorsList,N=256,gamma=1)
    
    im = ax.imshow(data, **kwargs, vmin = -1, vmax = 1, cmap = CustomCmap)
    ax.yaxis.set_label_coords(-0.1,10)

    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    
    fraction = 1

    cbar = ax.figure.colorbar(im, ax = ax, **cbar_kw, norm = norm, cmap = CustomCmap, shrink = fraction*(6/len(row_labels)), pad = 0.05)
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
    # return also cbar
    return im, cbar


def annotate_heatmap(image, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    # Create annotation for heatmap to compare relative abundance values
    
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

def get_text_size(fontsize,text):
    ffam = 'DejaVu Sans'
    fp = matplotlib.font_manager.FontProperties(
    family=ffam, style='normal', size=fontsize,
    weight='normal', stretch='normal')
    pth = matplotlib.textpath.TextPath((100,100), text, prop=fp)
    bb = pth.get_extents()

    return bb.width, bb.height

def get_relative_abundance(sample_sheet):
    norm_abund_dict = {}

    extensions = np.array([sample_sheet[group]["path"] for group in sample_sheet.keys()]).flatten()
    extensions = np.unique([path_features(i)["extension"] for i in extensions])

    if (len(extensions) > 1):
        raise ValueError("It is forbidden to combine different transcript estimation tools.. Please specify input files using same method: Cufflinks or Mix2")

    for group in sample_sheet.keys():
        ## LOAD MIX2 ##
        if (extensions[0] in ("tsv","csv")):
            SIRV_abundance = read_mix2(sample_sheet[group], format = 'fpkm_chn')
        elif (extensions[0] == "gtf"):
            SIRV_abundance = read_cufflinks(sample_sheet[group])
        else:
            raise ValueError("Unsupported input format..")

        averaged_SIRV_abundance = {t:np.mean(SIRV_abundance[t]) for t in sorted(SIRV_abundance.keys())}
        
        ## NORMALIZATION OF TRANSCRIPT FPKM ##
        expected_quantity = np.mean(list(averaged_SIRV_abundance.values()))

        if (expected_quantity == 0):
            {t: "NaN" for t in sorted(averaged_SIRV_abundance.keys())}
        else:
            norm_abund_dict[group] = {t:( averaged_SIRV_abundance[t]/expected_quantity ) for t in sorted(averaged_SIRV_abundance.keys())}
        
    return norm_abund_dict

def export_data(norm_abund_dict, output_path):
    groups = list(norm_abund_dict.keys())
    with open(os.path.join(output_path+"relative_concentration.tsv"), "w") as out_file:
        out_file.write("\t"+"\t".join(groups)+"\n")
        transcript_ids = sorted(norm_abund_dict[groups[0]].keys())
        for transcript in transcript_ids:
            out_file.write(transcript+"\t"+"\t".join([str(norm_abund_dict[group][transcript]) for group in groups])+"\n")

def create_sirvsuite_boxplot(relative_abundance, output_path, experiment_name = ""):

    fig, ax1 = plt.subplots(figsize=(12,5))
    #ax1.set_xscale('log')
    output_name = os.path.join(output_path, "SIRVsuite_boxplot.png")
    ax1.set_title(experiment_name+": SIRV E0 relative concentration distribution", size = 20, pad = 10, loc = 'center',
                fontdict = {'verticalalignment':'baseline'})

    heatmap_matrix = count_dict_to_matrix(relative_abundance)

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

    # meanprops = dict(marker='o',markeredgecolor='red', markerfacecolor = 'green', markersize = 8 , alpha = 0.7),
    ax1.plot([1,1],[0,100], color = 'darkblue', alpha = .6)

    ax1.set_yticklabels(relative_abundance.keys())
    plt.xlabel("relative SIRV transcript concentration")
    fig.savefig(output_name, dpi = 300, quality = 100, pad_inches = 0.2)

def count_dict_to_matrix(count_dict):
    conc_matrix = np.array([])
    groups = count_dict.keys()
    for group in groups:
        transcripts = sorted(count_dict[group].keys())
        conc_matrix = np.append(conc_matrix, [count_dict[group][t] for t in transcripts])

    return np.reshape(conc_matrix, (len(transcripts),len(groups)), order = 'F')

def create_sirvsuite_heatmap(relative_abundance, output_path, experiment_name = ""):

    groups = list(relative_abundance.keys())
    transcript_names = np.array((sorted(relative_abundance[groups[0]].keys())))

    genes = ["SIRV1","SIRV2","SIRV3","SIRV4","SIRV5","SIRV6","SIRV7"]

    fig = plt.figure(figsize=(6,40))

    vectorized_len = np.vectorize(len)
    max_length_index =  np.argmax(vectorized_len(groups))
    text_width, text_height = get_text_size(16, groups[max_length_index])
    plt.title(experiment_name+": SIRV E0 relative concentrations", pad = text_width + 10 + 25, fontdict = {'fontsize': 24})
    plt.axis('off')

    transcript_num = np.zeros(len(genes)).astype(int)
    for i in range(len(genes)):
        for transcript in transcript_names:
            if re.match(genes[i], transcript):
                transcript_num[i] += 1

    spec = gridspec.GridSpec(ncols=1, nrows = len(genes), height_ratios = transcript_num)

    for idx in range(len(genes)):

        filtered = [ i for i in range(len(transcript_names)) if bool(re.match(genes[idx],transcript_names[i])) ]

        heatmap_matrix = count_dict_to_matrix(relative_abundance) ## creating matrix for heatmap
        heatmap_matrix = heatmap_matrix[filtered,:]
        t_filtered = transcript_names[filtered]

        LFC = np.log2(heatmap_matrix)

        ax = fig.add_subplot(spec[idx])

        text_width, text_height = get_text_size(16, genes[idx])

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

        im, cbar = heatmap_generator(LFC, t_filtered, groups, ax=ax, cbarlabel="Log2 Fold Change ($c_r$)", no_x_ticks=draw_x_ticks)
        texts = annotate_heatmap(im, data = LFC, valfmt="{x:.2f}", threshold=-10)

    plt.margins(1,1)
    fig.tight_layout()

    output_name = os.path.join(output_path, "SIRVsuite_heatmap.png") 
    fig.savefig(output_name, dpi = 300, quality = 100, pad_inches = 0.2, bbox_inches = 'tight')


def load_sample_sheet(sample_sheet_path):
    if (os.path.exists(sample_sheet_path)):
        file_dict = {}
        with open(sample_sheet_path) as sheet_file:
            sheet_reader = csv.reader(sheet_file,delimiter = "\t")
            for row in sheet_reader:
                if row[2] not in file_dict.keys():
                    file_dict[row[2]] = {"path": [row[0]], "alias": [row[1]]}
                else:
                    file_dict[row[2]]["path"].append(row[0])
                    file_dict[row[2]]["alias"].append(row[1])
    else:
        raise ValueError("wrong input sample sheet")
    return file_dict




## PARSING ARUGMENTS ##
"""
parser = ap.ArgumentParser()
parser.add_argument('-f','--file-names', action = 'store', help = "", required = True, nargs = "+")
parser.add_argument('-o','--output-path', action = 'store', help = "", required = True, nargs = 1)
parser.add_argument('-n','--sample-names', action = 'store', help = "", required = False, nargs = "+")
parser.add_argument('-g','--sample-groups', action = 'store', help = "", required = False, nargs = "+")
parser.add_argument('--experiment-name',action = 'store', help = "", required = False, default = "", nargs = 1)
args = parser.parse_args()



## READING MIX2 FILES

sample_names = np.array(args.sample_names)
files = args.file_names
groups = np.array(args.sample_groups)

## Checking arguments ##
#if (len(files) != len(sample_names)):
#    raise Exception("Error while loading input..")
    
if (np.all(groups == None)):
    groups = sample_names
else:
    if (len(files) != len(groups)):
        raise Exception("Error while loading input..")

groups_unique = np.unique(groups)

"""

experiment_name = "FFPE_k2"

file_dict = load_sample_sheet("/home/tdrozd/development/sirv-suite/test/input/sheet_cufflinks.tsv")
relative_abundance = get_relative_abundance(file_dict)
#export_data(relative_abundance, "/home/tdrozd/development/sirv-suite/test/")
#create_sirvsuite_heatmap(relative_abundance, output_path = "/home/tdrozd/development/sirv-suite/test/", experiment_name = experiment_name)
#create_sirvsuite_boxplot(relative_abundance, output_path = "/home/tdrozd/development/sirv-suite/test/", experiment_name = experiment_name)