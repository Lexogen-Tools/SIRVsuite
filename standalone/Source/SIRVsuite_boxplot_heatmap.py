
import numpy as np
from os import path
import re
import csv 
import argparse as ap
import matplotlib.pyplot as plt
import matplotlib

def getAbundanceFromMIX2Table(file_path, format = 'fpkm'):
    
    SIRV_gene_set3 = ["SIRV1","SIRV2","SIRV3","SIRV4","SIRV5","SIRV6","SIRV7"]
    gene_selection = SIRV_gene_set3
    
    if path.exists(file_path):
        with open(file_path,'r') as MIX2_file:
            MIX2_reader = csv.reader(MIX2_file,delimiter = "\t")
            init = True
            for row in MIX2_reader:
                row = np.array(row)
                if init:
                    columns = np.array(["tracking_ID","gene_ID","abundance","comp_frags_locus","FPKM_CHN"]) # trim to necessary info
                    indices = (row[:,None] == columns).argmax(axis=0) # find which indexes in row corresponds to desired columns
                    data_dict = {}
                    init = False
                
                if row[1] in gene_selection:
                    selected_columns = row[indices]
                    
                    if format.lower() == 'fpkm':
                        data_dict[selected_columns[0]] = float(selected_columns[4])  
                    elif format.lower() == 'counts':
                        data_dict[selected_columns[0]] = float(selected_columns[2] * selected_columns[3]) # multiplying abundance and comp_frag_locus result in counts for transcripts  
                    elif format.lower() == 'abundance':
                        data_dict[selected_columns[0]] = float(selected_columns[2])
                    else:
                        print ("unknown format specified.. supported options are FPKM or abundance")
                        break
                else:
                    continue
                
            return data_dict

def returnFeaturesFromPath(filePath):
    matching_pattern = "((?:(?:.*\\/)*(.*)\\/))*(?:(.*)\\.(.*))*"
    feature_match = re.match(matching_pattern, filePath)
    out_dict = {"parent_dir": feature_match.group(2),
                "file_name": feature_match.group(3),
                "extension": feature_match.group(4),
                "path": feature_match.group(1)}
    return out_dict

def is_valid_file(parser, arg):
    if not path.exists(arg):
        parser.error("The file %s does not exists!")
    else:
        return arg

def heatmap(data, row_labels, col_labels, title = "SIRVsuite heatmap", ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):

    if not ax:
        ax = plt.gca()

    colorsList = [(0,0,1),(0,1,0),(1,0,0)]
    CustomCmap = matplotlib.colors.LinearSegmentedColormap.from_list("my_cmap",colorsList,N=256,gamma=1)
    
    # interpolation='nearest' extent=[1,2,1,2]
    # Plot the heatmap
    im = ax.imshow(data, **kwargs, vmin = -1, vmax = 1, cmap = CustomCmap)
    ax.yaxis.set_label_coords(-0.1,10) # -0.1,10

    #ax.set_title(title, pad = 200, fontdict = {'fontsize': 24}, loc = 'right')

    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    # Create colorbar
    cbar = ax.figure.colorbar(im, aspect = 30, ax = ax, **cbar_kw, norm = norm, shrink = 0.3, cmap = CustomCmap)
    cbar.ax.tick_params(labelsize=14) 
    cbar.ax.set_ylabel(cbarlabel, rotation = -90, va = "top", fontsize = 24, labelpad = 30)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False, labelsize = 16)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-90, ha="center",
             rotation_mode="default")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

def getTextSize(fontsize,text):
    ffam = 'DejaVu Sans'
    fp = matplotlib.font_manager.FontProperties(
    family=ffam, style='normal', size=fontsize,
    weight='normal', stretch='normal')
    pth = matplotlib.textpath.TextPath((100,100), text, prop=fp)
    bb = pth.get_extents()

    return bb.width, bb.height

## PARSING ARUGMENTS ##

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

norm_abund_dict = {}

it = 0

for group in groups_unique:
    group_indices = np.where(groups == group)[0]
    
    for i in range(len(group_indices)):
        ## LOAD MIX2 ##
    
        SIRV_abundance = getAbundanceFromMIX2Table(files[group_indices[i]], format = 'FPKM')
        
        ## NORMALIZATION OF TRANSCRIPT FPKM ##

        if (it > 0):
            if SIRV_transcripts != sorted(SIRV_abundance):
                raise Exception("An error occured while loading abundance files. Probably SIRV transcript missing..")
        
        SIRV_transcripts = sorted(SIRV_abundance)

        SIRV_abund_list = np.array([SIRV_abundance[i] for i in SIRV_transcripts])
        expected_quantity = np.sum(SIRV_abund_list) / len(SIRV_abundance)
        if (i==0):
            SIRV_norm_abund = SIRV_abund_list / expected_quantity
        else:
            SIRV_norm_abund += SIRV_abund_list / expected_quantity
            if (i==len(group_indices) - 1):
                SIRV_norm_abund /= i+1
        
    norm_abund_dict[group] = SIRV_norm_abund

    it += 1

## HEATMAP

heatmap_list = list(norm_abund_dict.values())
heatmap_matrix = np.transpose(np.array(heatmap_list)) ## creating matrix for heatmap

LFC = np.log2(heatmap_matrix)

params = {'mathtext.default': 'regular' } 
plt.rcParams.update(params)

fig, ax = plt.subplots(figsize=(6,25))
im, cbar = heatmap(LFC, SIRV_transcripts, groups_unique, ax=ax, cbarlabel=r'$Log_{2}\/Fold\/Change$')
texts = annotate_heatmap(im, valfmt="{x:.2f}", threshold=-200, fontsize=14)


plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
plt.margins(0,0)
vectorized_len = np.vectorize(len)
max_length_index =  np.argmax(vectorized_len(groups))
text_width, text_height = getTextSize(16, groups[max_length_index])
plt.title(args.experiment_name[0], pad = text_width + 15 + 20, fontdict = {'fontsize': 30})
output_name = path.join(args.output_path[0], "SIRVsuite_heatmap.png") 

fig.savefig(output_name, dpi = 300, quality = 100, pad_inches = 0.2, bbox_inches = 'tight')

## BOXPLOTS

fig, ax1 = plt.subplots(figsize=(12,5))
ax1.set_xscale('log')
output_name = path.join(args.output_path[0], "SIRVsuite_boxplot.png")
ax1.set_title(args.experiment_name[0], size = 20, pad = 10, loc = 'center',
              fontdict = {'verticalalignment':'baseline'})

for i in range(0,len(heatmap_list)):
    x = np.ones(len(heatmap_list[i]))*(i+1)
    ax1.scatter(heatmap_list[i], x, alpha = 0.5, s = 10, color = 'gray')

ax1.boxplot(heatmap_list,
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

ax1.set_yticklabels(groups_unique)
plt.xlabel("relative SIRV transcript concentration")
fig.show()
fig.savefig(output_name, dpi = 300, quality = 100, pad_inches = 0.2)