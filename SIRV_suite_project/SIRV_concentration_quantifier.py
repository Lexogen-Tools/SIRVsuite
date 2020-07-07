import numpy as np
import argparse as ap
from os import path
import csv
import re
from pyexpat import features
from _csv import writer

def getAbundanceFromMIX2Table(file_path, format = 'fpkm'):
    if path.exists(file_path):
        with open(file_path,'r') as MIX2_file:
            MIX2_reader = csv.reader(MIX2_file,delimiter = "\t")
            it = 0
            init = True
            
            for row in MIX2_reader:
                
                row = np.array(row)
                
                if init:
                    columns = np.array(["tracking_ID","abundance","comp_frags_locus","FPKM_CHN"]) # trim to necessary info
                    indices = (row[:,None] == columns).argmax(axis=0) # find which indexes in row corresponds to desired columns
                    data_dict = {}
                    data_dict[columns[0]] = np.array([])
                    data_dict["quantity"] = np.array([])
                    init = False
                
                if row[0].startswith(("SIRV","ERCC")):
                    selected_columns = row[indices]
                    data_dict[columns[0]] = np.append(data_dict[columns[0]], selected_columns[0])
                    
                    if format.lower() == 'fpkm':
                        data_dict["quantity"] = np.append(data_dict["quantity"], float(selected_columns[3]))  
                    elif format.lower() == 'counts':
                        data_dict["quantity"] = np.append(data_dict["quantity"], float(selected_columns[1] * selected_columns[2])) # multiplying abundance and comp_frag_locus result in counts for transcripts  
                    elif format.lower() == 'abundance':
                        data_dict["quantity"] = np.append(data_dict["quantity"], float(selected_columns[1])) 
                    else:
                        print ("unknown format specified.. supported options are FPKM or abundance")
                        break
                    it += 1
                    
                else:
                    continue
                
            return data_dict

def returnFeaturesFromPath(filePath):
    matching_pattern = "(?:(?:.*\\/)*(.*)\\/)*(.*)\\.(.*)"
    feature_match = re.match(matching_pattern, filePath)
    out_dict = {"parent_dir": feature_match.group(1),
                "file_name": feature_match.group(2),
                "extension": feature_match.group(3)}
    return out_dict

def is_valid_file(parser, arg):
    if not path.exists(arg):
        parser.error("The file %s does not exists!")
    else:
        return arg

parser = ap.ArgumentParser()
parser.add_argument('-n','--normalization', nargs = '?', default = "nm1", help = "")
# parser.add_argument('-i','--sample-list', nargs = '?', required = False, type = lambda x: is_valid_file(parser, x), action = "store", metavar = "FILE")
parser.add_argument('-a','--average-replicates', action = 'store_true', help = "")
parser.add_argument('-o','--output-name', action = 'store', help = "", required = True)
parser.add_argument('-f','--file-name', action = 'store', help = "", required = True)
parser.add_argument('-s','--sample-name', action = 'store', help = "", required = True)

args = parser.parse_args()

# with open(args.one_file_only,'r') as csvfile_input:
#csv_file_input = open(args.one_file_only,'r')
#reader = csv.reader(csv_file_input, delimiter = ';')
SIRV_abundance_dict = {}
#conditions = np.array([])

#for row in reader:
    
#    if len(row) == 0:
#        continue
    
#    if row[0].startswith(('#','%',"-")) or len(row) == 0:
#        continue
#    else:
#        conditions = np.append(conditions, row[1])
# features = returnFeaturesFromPath(args.one_file_only)
sample_name = args.sample_name
SIRV_abundance_data = getAbundanceFromMIX2Table(args.file_name, format = 'fpkm')
SIRV_abundance_dict[sample_name] = SIRV_abundance_data
        
        # HERE GOES NORMALIZATION
        
expected_quantity = np.sum(SIRV_abundance_data["quantity"]) / len(SIRV_abundance_dict[sample_name]["tracking_ID"])
SIRV_abundance_dict[sample_name]["norm_abund"] = SIRV_abundance_dict[sample_name]["quantity"] / expected_quantity

with open(args.output_name, 'w', newline='') as csv_file_output:
    writer = csv.writer(csv_file_output, delimiter = ";")
    header_row = [""]
    header_row.extend(list(SIRV_abundance_dict))
    writer.writerow(header_row)
    
    for row in range(0,len(SIRV_abundance_dict[list(SIRV_abundance_dict)[0]]["norm_abund"])):
        k = [SIRV_abundance_dict[list(SIRV_abundance_dict)[0]]["tracking_ID"][row]]
        
        for sample in list(SIRV_abundance_dict):
            k.append(SIRV_abundance_dict[sample]["norm_abund"][row])
        
        writer.writerow(k)
        
