import numpy as np
import csv 
import argparse as ap
import matplotlib.pyplot as plt

parser = ap.ArgumentParser()
parser.add_argument('-f','--file-names', action = 'store', help = "", required = True, nargs = "+")
parser.add_argument('-o','--output-name', action = 'store', help = "", required = True)
parser.add_argument('--title-name', action = 'store', help = "", required = False, default = "SIRV concentration")

args = parser.parse_args()

value_list = []
key_names = []
in_dict = dict()

for file_it,file in enumerate(args.file_names):
    k = open(file, 'r')
    r = csv.reader(k, delimiter=";")

    for row_it, row in enumerate(r):
        if row_it == 0:
            key = row[1]
            in_dict[key] = np.array([])
            key_names.append(key)
        else:
            in_dict[key] = np.append(in_dict[key], float(row[1]))
    
plot_list = []
for key in key_names:
    plot_list.append(in_dict[key])
    
tog = plot_list

fig, ax1 = plt.subplots(figsize=(12,5))
ax1.set_xscale('log')
ax1.set_title(args.title_name, size = 20, pad = 10, loc = 'center',
              fontdict = {'verticalalignment':'baseline'})

for i in range(0,len(tog)):
    x = np.ones(len(tog[i]))*(i+1)
    ax1.scatter(tog[i], x, alpha = 0.5, s = 10, color = 'gray')

ax1.boxplot(tog,
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

ax1.set_yticklabels(key_names)
plt.xlabel("relative SIRV transcript concentration")
fig.show()
fig.savefig(args.output_name, dpi = 300, quality = 100, pad_inches = 0.2)