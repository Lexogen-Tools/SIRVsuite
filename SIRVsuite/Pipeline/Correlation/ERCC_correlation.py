import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import csv
import os
import numpy as np
import math
from ..countReader import countReader
import logging

log = logging.getLogger(__name__.split(".")[-1])

class ERCCcorrelation():

    def __init__(self, sample_sheet = None, output_dir = "./", experiment_name = ""):
        """
        """
        
        if sample_sheet is not None:
            c = countReader()
            cnts = dict()
            for sample in sample_sheet.keys():    
                cnts[sample] = c.read_counting_file(files=[sample_sheet[sample]["counting_path"]],counting_method="mix2",spike_in_type=["ERCC"])
        
            self.ERCC_correlation(cnts, output_dir=os.path.join(output_dir,"correlation/"), experiment_name = experiment_name)

    def read_ERCC_concentration_table(self, ERCC_table_path):
        """
        Reading of ERCC concentration table
        """
        
        conc = {}

        with open(ERCC_table_path, "r") as f:
            reader = csv.DictReader(f, delimiter = "\t")
            for row in reader:
                conc[row["ERCC ID"]] = row
        return (conc, reader.fieldnames[2:])
            
    def ERCC_correlation(self, data, experiment_name = "", output_dir = "./"):
        """
        This functions loads list of input files (gene or transcript counts), type of quantification, names of samples, type of ERCC spike mix (Mix1 or Mix2) 
        and experiment name, which is optional   
        """

        log.info("calculating correlation of ERCC spike-in ratios")

        ERCC_table_path = os.path.dirname(__file__)
        ERCC_table_path = ERCC_table_path.replace("Pipeline/Correlation","Resources/ERCC92_Concentration.tsv")

        ERCC_conc,_ = self.read_ERCC_concentration_table(ERCC_table_path)
        ercc_spike_in = "Mix1"

        #cnt_data_total = get_data(file_list, types, sample_names)
        if experiment_name != "":
            experiment_name = experiment_name + ": "

        data_path = os.path.join(output_dir, "ERCC_correlation.tsv")

        corr_per_sample = {}

        for sample_name in data.keys():

            cnt_data = data[sample_name]["ERCC"]

            y = [cnt_data[ctrl][0] for ctrl in sorted(ERCC_conc.keys())]
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
            fig.gca().title.set_text("%s %s (R=%.3f)"%(experiment_name,sample_name,corr))

            nested_dir = os.path.join(output_dir, sample_name)

            if ( not os.path.exists(nested_dir)):
                os.makedirs(nested_dir)

            plot_path = os.path.join(nested_dir, "ERCC_correlation.png")
            fig.savefig(plot_path, dpi = 200)
            plt.close(fig)
        
        # write summary file for all given samples
        with open(data_path,"w") as summary_file:
            summary_file.write("Alias\tERCC Pearson R\n")
            for sample_name in data.keys():
                summary_file.write("%s\t%f\n"%(sample_name, corr_per_sample[sample_name]))
