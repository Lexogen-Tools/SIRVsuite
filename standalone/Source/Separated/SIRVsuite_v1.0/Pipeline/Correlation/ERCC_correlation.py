import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import csv
import os
import numpy as np
import math

def read_ERCC_concentration_table():
    # Reading of ERCC concentration table
    conc = {}
    ERCC_concentration_file = "Resources/ERCC92_Concentration.tsv"
    with open(ERCC_concentration_file, "r") as f:
        reader = csv.DictReader(f, delimiter = "\t")
        for row in reader:
            conc[row["ERCC ID"]] = row
    return (conc, reader.fieldnames[2:])

def get_data(file_list, file_type, sample_names):
    # TODO: This function needs to be ideally generalized due to usage for ERCC correlation + SIRV concentration  
    # The function returns a dictionary with following levels: sample_names -> gene or tracking_ID -> counts

    cnts_total = {}
    for index, path in enumerate(file_list):
        cnts=dict()
        if(file_type[index] == "COUNTING"):
            with open(path, "r") as f:
                for line in f:
                    gene, count = line.strip().split("\t")
                    if(gene[0:2] != "__"):
                        cnts[gene] = float(count)
        elif(file_type[index] == "MIX2"):
            with open(path, "r") as f:
                for line in f:
                    if(line.strip().split("\t")[0]=="tracking_ID"):
                        continue
                    # extract all the values, for now
                    _,_,_,_,_,_,_,_,_,FPKM_THN,_,_,_,_,_,_,_,_,_,_ = line.strip().split("\t")
                    
                    if("ERCC" in gene_ID):
                        cnts[gene_ID]=float(FPKM_THN)
                    else:
                        cnts[tracking_ID]=float(FPKM_THN)
        else:
            raise "Uknown type"
        cnts_total[sample_names[index]] = cnts
    return cnts_total

        
def ERCC_correlation(file_list, types, sample_names, ercc_spike_in, experiment_name = "unnamed"):
    # This functions loads list of input files (gene or transcript counts), type of quantification, names of samples, type of ERCC spike mix (Mix1 or Mix2) 
    # and experiment name, which is optional    

    ERCC_conc,_ = read_ERCC_concentration_table()

    output_dir = "test_output/Correlation/"

    cnt_data_total = get_data(file_list, types, sample_names)

    data_path = os.path.join(output_dir, "spikein_stats_%s.tsv"%(experiment_name))

    corr_per_sample = {}

    for sample_name in sample_names:

        cnt_data = cnt_data_total[sample_name]

        y = [cnt_data[ctrl] for ctrl in sorted(ERCC_conc.keys())]
        cnt = float(sum(y))
        y = [y_i/cnt if cnt>0 else 1e-100 for y_i in y]

        x = [float(ERCC_conc[ctrl][ercc_spike_in]) for ctrl in sorted(ERCC_conc.keys())]
        cnt = float(sum(x))
        x = [x_i/cnt for x_i in x]
        fig = plt.figure()
        fig.gca().loglog(x,y, "rx")
        x = np.array(x)
        y = np.array(y)
        ymin = 10**math.floor(math.log10(min(list(x[y>0])+list(y[y>0])))) # calculating axis limit 
        ymax = 10**math.ceil(math.log10(max(list(x[y>0])+list(y[y>0]))))
        fig.gca().set_xlabel("Theoretical concentration")
        fig.gca().set_ylabel("Measured concentration")
        fig.gca().set_xlim(ymin, ymax) 
        fig.gca().set_ylim(ymin, ymax)
        corr = pearsonr([math.log10(x[i]) for i in range(len(x)) if x[i]!=0 and y[i]!=0],
                                [math.log10(y[i]) for i in range(len(y)) if x[i]!=0 and y[i]!=0])[0]
        corr_per_sample[sample_name] = corr
        fig.gca().title.set_text("Sample %s (R=%.3f)"%(sample_name,corr))

        nested_dir = os.path.join(output_dir, sample_name)

        if ( not os.path.exists(nested_dir)):
            os.makedirs(nested_dir)

        plot_path = os.path.join(nested_dir, "ERCC_%s.png"%(sample_name))
        fig.savefig(plot_path)
        plt.close(fig)
    
    # write summary file for all given samples
    with open(data_path,"w") as summary_file:
        summary_file.write("Alias\tERCC\n")
        for sample_name in sample_names:
            summary_file.write("%s\t%f\n"%(sample_name, corr))
