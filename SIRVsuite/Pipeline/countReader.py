import csv
import numpy as np
import re
import logging

log = logging.getLogger(__name__.split(".")[-1])

class countReader():
    def __init__(self, sample_sheet = None):
        self.spike_in_name_pattern = {"SIRV":{"gene":"","transcript":""}, "ERCC":{"gene":"","transcript":""}}
        #self.spike_in_ids["SIRV"]["gene"] = ["SIRV1","SIRV2","SIRV3","SIRV4","SIRV5","SIRV6","SIRV7"]
        #self.spike_in_ids["SIRV"]["transcript"] = ["SIRV101", "SIRV101", "SIRV101", "SIRV101", "SIRV101", "SIRV101", "SIRV102", "SIRV102", "SIRV102", "SIRV102", "SIRV103", "SIRV103", "SIRV103", "SIRV103", "SIRV103", "SIRV103", "SIRV105", "SIRV105", "SIRV105", "SIRV105", "SIRV105", "SIRV106", "SIRV106", "SIRV106", "SIRV107", "SIRV107", "SIRV107", "SIRV108", "SIRV108", "SIRV108", "SIRV109", "SIRV109", "SIRV109", "SIRV201", "SIRV201", "SIRV201", "SIRV201", "SIRV201", "SIRV201", "SIRV201", "SIRV201", "SIRV201", "SIRV201", "SIRV201", "SIRV202", "SIRV202", "SIRV202", "SIRV202", "SIRV202", "SIRV202", "SIRV202", "SIRV202", "SIRV202", "SIRV202", "SIRV202", "SIRV203", "SIRV203", "SIRV203", "SIRV203", "SIRV203", "SIRV204", "SIRV204", "SIRV204", "SIRV205", "SIRV206", "SIRV301", "SIRV301", "SIRV301", "SIRV301", "SIRV301", "SIRV302", "SIRV302", "SIRV303", "SIRV303", "SIRV303", "SIRV304", "SIRV304", "SIRV304", "SIRV304", "SIRV304", "SIRV304", "SIRV304", "SIRV304", "SIRV305", "SIRV305", "SIRV305", "SIRV306", "SIRV306", "SIRV306", "SIRV307", "SIRV307", "SIRV307", "SIRV307", "SIRV307", "SIRV308", "SIRV308", "SIRV308", "SIRV309", "SIRV309", "SIRV309", "SIRV310", "SIRV310", "SIRV310", "SIRV311", "SIRV403", "SIRV403", "SIRV403", "SIRV403", "SIRV404", "SIRV404", "SIRV404", "SIRV404", "SIRV405", "SIRV405", "SIRV406", "SIRV406", "SIRV408", "SIRV408", "SIRV408", "SIRV408", "SIRV408", "SIRV409", "SIRV409", "SIRV409", "SIRV410", "SIRV410", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV501", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV502", "SIRV503", "SIRV503", "SIRV503", "SIRV504", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV505", "SIRV506", "SIRV506", "SIRV507", "SIRV507", "SIRV507", "SIRV507", "SIRV507", "SIRV507", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV508", "SIRV509", "SIRV509", "SIRV509", "SIRV509", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV510", "SIRV511", "SIRV511", "SIRV512", "SIRV601", "SIRV601", "SIRV601", "SIRV601", "SIRV601", "SIRV601", "SIRV601", "SIRV601", "SIRV601", "SIRV602", "SIRV602", "SIRV602", "SIRV602", "SIRV602", "SIRV602", "SIRV602", "SIRV602", "SIRV603", "SIRV604", "SIRV604", "SIRV604", "SIRV604", "SIRV604", "SIRV604", "SIRV604", "SIRV604", "SIRV604", "SIRV604", "SIRV605", "SIRV605", "SIRV605", "SIRV605", "SIRV605", "SIRV605", "SIRV605", "SIRV605", "SIRV605", "SIRV606", "SIRV606", "SIRV606", "SIRV606", "SIRV607", "SIRV607", "SIRV607", "SIRV607", "SIRV608", "SIRV608", "SIRV608", "SIRV608", "SIRV609", "SIRV609", "SIRV609", "SIRV609", "SIRV610", "SIRV610", "SIRV610", "SIRV610", "SIRV610", "SIRV611", "SIRV611", "SIRV611", "SIRV612", "SIRV612", "SIRV612", "SIRV612", "SIRV612", "SIRV612", "SIRV612", "SIRV612", "SIRV612", "SIRV612", "SIRV613", "SIRV613", "SIRV613", "SIRV613", "SIRV613", "SIRV613", "SIRV614", "SIRV614", "SIRV614", "SIRV614", "SIRV614", "SIRV615", "SIRV615", "SIRV615", "SIRV616", "SIRV616", "SIRV616", "SIRV616", "SIRV617", "SIRV618", "SIRV701", "SIRV701", "SIRV701", "SIRV701", "SIRV701", "SIRV702", "SIRV702", "SIRV702", "SIRV702", "SIRV702", "SIRV702", "SIRV703", "SIRV703", "SIRV703", "SIRV703", "SIRV703", "SIRV704", "SIRV704", "SIRV704", "SIRV705", "SIRV705", "SIRV705", "SIRV705", "SIRV705", "SIRV706", "SIRV706", "SIRV706", "SIRV706", "SIRV706", "SIRV708", "SIRV708", "SIRV708", "SIRV708", "SIRV708", "SIRV708"]
        #self.spike_in_ids["ERCC"]["transcript"] = ["DQ459430", "DQ516784", "DQ516752", "DQ668364", "DQ883670", "EF011062", "DQ875385", "DQ883664", "DQ459420", "DQ883651", "DQ855004", "DQ854993", "DQ883689", "DQ459419", "DQ459431", "DQ516796", "DQ855001", "DQ459413", "DQ883656", "DQ883661", "EF011069", "DQ516783", "DQ516787", "DQ459424", "DQ516748", "DQ883671", "DQ516740", "DQ516785", "DQ516731", "DQ668366", "DQ459418", "DQ668356", "DQ516763", "DQ459426", "DQ516786", "DQ883653", "DQ459421", "DQ883654", "DQ668358", "DQ516754", "DQ516778", "DQ883650", "DQ516742", "DQ883673", "DQ883652", "DQ854991", "DQ516780", "DQ883682", "DQ883669", "DQ516791", "DQ459425", "DQ516759", "DQ459429", "DQ516758", "DQ459415", "DQ875387", "DQ516815", "DQ668365", "DQ854998", "DQ883685", "DQ459422", "DQ883663", "DQ668367", "DQ459412", "DQ854992", "DQ516782", "DQ459427", "EF011072", "DQ855003", "DQ516739", "EF011063", "DQ855000", "DQ516777", "DQ883646", "DQ668362", "DQ854995", "DQ875386", "DQ516790", "DQ883642", "DQ883659", "DQ854997", "DQ883643", "DQ839618", "DQ516795", "DQ883658", "DQ516750", "DQ668359", "DQ516779", "DQ668363", "DQ516776", "DQ516773", "DQ854994"]
        #self.spike_in_ids["ERCC"]["gene"] = ["ERCC-00002", "ERCC-00003", "ERCC-00004", "ERCC-00009", "ERCC-00012", "ERCC-00013", "ERCC-00014", "ERCC-00016", "ERCC-00017", "ERCC-00019", "ERCC-00022", "ERCC-00024", "ERCC-00025", "ERCC-00028", "ERCC-00031", "ERCC-00033", "ERCC-00034", "ERCC-00035", "ERCC-00039", "ERCC-00040", "ERCC-00041", "ERCC-00042", "ERCC-00043", "ERCC-00044", "ERCC-00046", "ERCC-00048", "ERCC-00051", "ERCC-00053", "ERCC-00054", "ERCC-00057", "ERCC-00058", "ERCC-00059", "ERCC-00060", "ERCC-00061", "ERCC-00062", "ERCC-00067", "ERCC-00069", "ERCC-00071", "ERCC-00073", "ERCC-00074", "ERCC-00075", "ERCC-00076", "ERCC-00077", "ERCC-00078", "ERCC-00079", "ERCC-00081", "ERCC-00083", "ERCC-00084", "ERCC-00085", "ERCC-00086", "ERCC-00092", "ERCC-00095", "ERCC-00096", "ERCC-00097", "ERCC-00098", "ERCC-00099", "ERCC-00104", "ERCC-00108", "ERCC-00109", "ERCC-00111", "ERCC-00112", "ERCC-00113", "ERCC-00116", "ERCC-00117", "ERCC-00120", "ERCC-00123", "ERCC-00126", "ERCC-00130", "ERCC-00131", "ERCC-00134", "ERCC-00136", "ERCC-00137", "ERCC-00138", "ERCC-00142", "ERCC-00143", "ERCC-00144", "ERCC-00145", "ERCC-00147", "ERCC-00148", "ERCC-00150", "ERCC-00154", "ERCC-00156", "ERCC-00157", "ERCC-00158", "ERCC-00160", "ERCC-00162", "ERCC-00163", "ERCC-00164", "ERCC-00165", "ERCC-00168", "ERCC-00170", "ERCC-00171"]
        
        self.spike_in_name_pattern["SIRV"]["gene"] = "^SIRV\d{1}$"
        self.spike_in_name_pattern["SIRV"]["transcript"] = "^SIRV\d{3}$"
        self.spike_in_name_pattern["ERCC"]["gene"] = "^ERCC-\d{5}$"
        self.spike_in_name_pattern["ERCC"]["transcript"] = "(DQ|EF)\d{6}"

    def read_counting_file(self, files = [], counting_method="mix2", counting_type="gene", spike_in_type=["ERCC","SIRV"]):
        
        counting_dict = {s:dict() for s in spike_in_type}
        
        if counting_method == "mix2":
            
            for path in files:
                
                try:
                    quantity_unit = "fpkm_chn"
                    with open(path,'r') as MIX2_file:
                        MIX2_reader = csv.reader(MIX2_file,delimiter = "\t")
                        init = True
                        
                        for row in MIX2_reader:
                            row = np.array(row)
                            if init:
                                required_cols = np.array(["tracking_ID", "gene_ID", "FPKM_CHN"]) # trim to necessary info
                                
                                if (not set(required_cols) <= set(row)):
                                    raise ValueError("Cannot load mix2 table format..")
                                
                                indices = (row[:,None] == required_cols).argmax(axis=0) # find which indexes in row corresponds to desired columns
                                init = False
                                continue
                            
                            if re.match(self.spike_in_name_pattern["SIRV"]["gene"],row[1]) and ("SIRV" in spike_in_type):
                                identifier_idx = 0
                                spike_in = "SIRV"
                            
                            elif re.match(self.spike_in_name_pattern["ERCC"]["gene"],row[1]) and ("ERCC" in spike_in_type):
                                identifier_idx = 1
                                spike_in = "ERCC"
                            else:
                                continue

                            selected_columns = row[indices]

                            if quantity_unit.lower() == 'fpkm_chn':
                                value = float(selected_columns[2])
                            else:
                                log.warning("unknown quantity unit specified.. supported options are fpkm_chn or abundance")
                                break
                            
                            if selected_columns[identifier_idx] not in counting_dict[spike_in].keys():
                                counting_dict[spike_in][selected_columns[identifier_idx]] = [value]
                            else:
                                counting_dict[spike_in][selected_columns[identifier_idx]].append(value)
                except:
                    raise ValueError("An error occured during loading file: %s"%(path))
                    
        else:
            raise ValueError("Unsupported counting method format")

        """
        elif counting_method == "htseq":
            for path in files:
                try:
                    with open(path, "r") as f:
                            for line in f:
                                id, count = line.strip().split("\t")
                                
                                if "ERCC" in spike_in_type:
                                    if ( (re.match(self.spike_in_name_pattern["ERCC"]["gene"],id)) or (re.match(self.spike_in_name_pattern["ERCC"]["transcript"],id))):
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
        """

        return counting_dict