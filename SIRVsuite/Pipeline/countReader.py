import csv
import numpy as np
import re
import logging
import os
import gtfparse

log = logging.getLogger(__name__.split(".")[-1])

class countReader():
    # A class to read the count files for SIRVs and ERCC transcripts
    
    def __init__(self):
        self.spike_in_name_pattern = {"SIRV":{"gene":"","transcript":""}, "ERCC":{"gene":"","transcript":""}}
        
        # Provide regex gene or transcript IDs
        self.spike_in_name_pattern["SIRV"]["gene"] = "^SIRV\\d{1}$"
        self.spike_in_name_pattern["SIRV"]["transcript"] = "^SIRV\\d{3}$"
        self.spike_in_name_pattern["ERCC"]["gene"] = "^ERCC-\\d{5}$"
        self.spike_in_name_pattern["ERCC"]["transcript"] = "(DQ|EF)\\d{6}"

        self._annotation_path = "/".join(__file__.split("/")[:-3]) + "/SIRVsuite/Resources/SIRV_isoforms_multi-fasta-annotation_190418a.gtf"


    def read_counting_file(self, files=[], counting_method="mix2", spike_in_type=["ERCC","SIRV"], counting_type="transcript", q_unit="fpkm_chn"):
        
        counting_dict = {s:dict() for s in spike_in_type}
        # For now we support mix2 table format only. Additional format support can be added in later versions.
        if counting_method.lower() == "mix2":
            
            for path in files:
                try:
                    with open(os.path.abspath(path),'r') as MIX2_file:
                        MIX2_reader = csv.reader(MIX2_file,delimiter = "\t")
                        init = True
                        
                        for row in MIX2_reader:
                            row = np.array(row)
                            if init:
                                required_cols = np.array(["tracking_id", "gene_id", q_unit]) # trim to necessary info
                                row_lowercase = [r.lower() for r in row]
                                if (not set(required_cols) <= set(row_lowercase)):
                                    raise ValueError("Cannot load mix2 table format..")
                                
                                indices = [idx for col_r in required_cols for idx, col_i in enumerate(row_lowercase) if col_r==col_i] # find which indexes in row corresponds to desired columns
                                init = False
                                continue
                            
                            if re.match(self.spike_in_name_pattern["SIRV"][counting_type],row[indices[0]]) and ("SIRV" in spike_in_type):
                                identifier_idx = 0
                                spike_in = "SIRV"
                            
                            elif re.match(self.spike_in_name_pattern["ERCC"][counting_type],row[indices[0]]) and ("ERCC" in spike_in_type):
                                identifier_idx = 1
                                spike_in = "ERCC"
                            else:
                                continue

                            selected_columns = row[indices]

                            value = float(selected_columns[2])
                            '''
                            if q_unit.lower() == 'fpkm_chn':
                                value = float(selected_columns[2])
                            else:
                                log.warning("unknown quantity unit specified.. supported options are fpkm_chn or abundance")
                                break
                            '''
                            
                            if selected_columns[identifier_idx] not in counting_dict[spike_in].keys():
                                counting_dict[spike_in][selected_columns[identifier_idx]] = [value]
                            else:
                                counting_dict[spike_in][selected_columns[identifier_idx]].append(value)
                except:
                    raise ValueError("An error occured during loading file: %s"%(path))
                
        elif counting_method.lower() == "htseq":
            for path in files:
                try:
                    if "ERCC" in spike_in_type and counting_type=='transcript':
                        gtf_dict = gtfparse.read_gtf(self._annotation_path)
                        gtf_dict = gtf_dict[gtf_dict['seqname'].str.startswith("ERCC")][["gene_id", "transcript_id"]]

                    with open(path, "r") as f:
                            for line in f:
                                id, count = line.strip().split("\t")
                                
                                if "ERCC" in spike_in_type:
                                    if re.match(self.spike_in_name_pattern["ERCC"][counting_type], id):
                                        
                                        if (counting_type == "transcript"):
                                            id = gtf_dict[gtf_dict['transcript_id']==id].gene_id.item()

                                        if id not in counting_dict["ERCC"].keys(): 
                                            counting_dict["ERCC"][id] = [float(count)]
                                        else:
                                            counting_dict["ERCC"][id].append(float(count))

                                if "SIRV" in spike_in_type:
                                    if re.match(self.spike_in_name_pattern["SIRV"][counting_type], id):
                                        if id not in counting_dict["SIRV"].keys(): 
                                            counting_dict["SIRV"][id] = [float(count)]
                                        else:
                                            counting_dict["SIRV"][id].append(float(count))
                except:
                    raise ValueError("An error occured during loading file: %s"%(path))
                
        else:
            raise ValueError("Unsupported counting method format")

        return counting_dict