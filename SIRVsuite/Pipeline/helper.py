import re
import csv
import sys
import logging

log = logging.getLogger(__name__.split(".")[-1])

def path_features(filePath):
    """
    A function for creating dictionary of parent directory name, file name, extension and path to the file directory using regex.
    """
    matching_pattern = "((?:(?:.*\\/)*(.*)\\/))*(?:(.*)\\.(.*))*"
    feature_match = re.match(matching_pattern, filePath)
    out_dict = {"parent_dir": feature_match.group(2),
                "file_name": feature_match.group(3),
                "extension": feature_match.group(4),
                "path": feature_match.group(1)}
    return out_dict

def __entry_check__(valid_dict, colname, value, module):
    if colname in valid_dict[module].keys():
        if value:
            if value.lower() not in valid_dict[module][colname]:
                raise ValueError("Invalid entry %s for %s column.. possible options for module %s are: %s"%(value, colname, module, ", ".join(valid_dict[module][colname])))
            else:
                return value.lower()
        else:
            raise ValueError("%s column not defined for all rows.. please define it and insert one of the possible values: %s"%(colname,", ".join(valid_dict[colname])))
    else:
        raise ValueError("colname not in valid_dict")

def read_sample_sheet(sheet_path, modules_to_execute = ["SIRV-concentration", "coverage", "ERCC-correlation"]):
    sample_sheet_info = """
    Sample sheet is required to have a following format:
        - used ";" separator
        - trailing whitespace or tab is allowed (trimmed during reading process)
        - UTF-8 encoding
        - every column must have a name which is predefined (see more info below)

    Allowed column names:
        General columns:
            sample_name - Supported values: any set of characters to identify samples (this will be printed in the final graphics).
            
            library_prep_type - Supported values: whole (whole transcriptome library prep) or qs (QuantSeq library preparation) 

        Concentration specific columns:
            counting_path - Supported values: valid path to count files.
            
            counting_method - Supported values: mix2, cufflinks or htseq. Defines, how the count file should be read.
            
            counting_feature - Supported values: gene or transcript. Defines, whether ERCC correlation plots and tables are quantified (gene counts) or 
                               ERCC correlation + SIRV heatmap and boxplot are quantified (transcript counts).
            
            replicate_group - Supported values: any set of characters to identify replicate groups (this will be printed in the final graphics). 
                              An optional column, if the same replicate group assigned to multiple samples, their mean value will be used for quantification 
    
        Coverage specific columns:
            alignment_path - Supported values: valid paths to a .bam file.
            
            read_orientation - Supported values: fwd, rev or none. Use "fwd" or "rev" for strand-specific libraries, "none" for non-strand specific libraries. 

    Any other columns will be ignored.
    """

    log.info("reading sample sheet")

    ## ALLOWED COLUMNS FOR DIFFERENT MODULES CAN BE DEFINED HERE ##

    # name required columns
    required_cols = dict()
    required_cols["SIRV-concentration"] = ["sample_name", "counting_path", "counting_feature"]
    required_cols["ERCC-correlation"] = ["sample_name", "counting_path", "counting_feature"]
    required_cols["coverage"] = ["sample_name", "alignment_path", "read_orientation"]

    # name optional columns
    optional_cols = dict()
    optional_cols["SIRV-concentration"] = ["replicate_group"]
    optional_cols["ERCC-correlation"] = []
    optional_cols["coverage"] = []

    ## -------------------------------------------------------- ##

    ## RESTRICTIONS FOR THE COLUMNS CAN BE DEFINED HERE ##
    
    value_restriction = dict()
    value_restriction["coverage"] = dict()
    value_restriction["coverage"]["read_orientation"] = ["fwd", "rev"]
    
    value_restriction["SIRV-concentration"] = dict()
    value_restriction["SIRV-concentration"]["counting_method"] = ["mix2", "htseq"]
    value_restriction["SIRV-concentration"]["counting_feature"] = ["transcript"]
    
    value_restriction["ERCC-correlation"] = dict()
    value_restriction["ERCC-correlation"]["counting_method"] = ["mix2", "htseq"]
    value_restriction["ERCC-correlation"]["counting_feature"] = ["gene","transcript"]

    ## ------------------------------------------------ ##

    sample_sheet_dict = dict()
  
    # select only columns for which a restriction has been defined
    restricted_cols = dict()
    for module in required_cols:
        restricted_cols[module] = [col for col in value_restriction[module].keys() if col in required_cols[module]]
    
    # select columns with defined restriction for all common modules
    common_restricted_cols = set()
    modules = list(required_cols.keys())
    for idx in range(len(modules)-1):
        common_restricted_cols = set(required_cols[modules[idx]]).intersection(set(required_cols[modules[idx+1]]))
    common_restricted_cols = [col for col in common_restricted_cols if col in value_restriction.keys()]

    try:
        with open(sheet_path, 'r') as sheet:

            sheet = filter(lambda row: row[0]!='#', sheet)

            csv.register_dialect('strip', skipinitialspace=True, delimiter = ';')

            header = [h.strip() for h in next(sheet).split(";")]

            # check if cols are complete for specific modes 
            available_modules = [module for module in required_cols.keys() if set(required_cols[module]) <= set(header)]

            sample_sheet_dict = {module:dict() for module in modules_to_execute if module in available_modules}
            
            reader = csv.DictReader(sheet, restkey=None, restval=None, dialect='strip', fieldnames=header)

            for row in reader:
                # delete first or last semicolon in the header
                if "" in row.keys():
                    del row['']

                if None in row.values():
                    raise ValueError("It seems there is a value missing in a column..")

                if None in row.keys():
                    raise ValueError("More columns than values in the row detected")

                # remove tabs
                row = {a:row[a].strip("\t ") for a in row.keys()}  
                
                # check for valid values for general columns (all modules)
                for common_col in common_restricted_cols:
                    __entry_check__(value_restriction, common_col, row.get(common_col))
                
                # check for valid values for a given column across modules
                for module in modules_to_execute:
                    if module in available_modules:
                        for col in restricted_cols[module]:
                            row[col] = __entry_check__(value_restriction, col, row.get(col), module = module)
                            sample_sheet_dict[module][row.get("sample_name")] = {val:row[val] for val in required_cols[module][1:] if val in row}

                        for col in optional_cols[module]:
                            if col in header:
                                sample_sheet_dict[module][row.get("sample_name")][col] = row.get(col)  
                    else:
                        raise ValueError("Cannot proceed with module: %s.. Please check required columns for the module in the sample sheet."%(module))
                    
                    if module == "concentration":
                        # c_methods = set([sample_sheet_dict[module][i]["counting_method"] for i in sample_sheet_dict[module].keys()])
                        c_features = set([sample_sheet_dict[module][i]["counting_feature"] for i in sample_sheet_dict[module].keys()])
                        # c_library_prep = set([sample_sheet_dict[module][i]["library_prep_type"] for i in sample_sheet_dict[module].keys()])
                        
                        # if len(c_methods) > 1:
                        #     raise ValueError("Cannot proceed with different counting_methods..")

                        if len(c_features) > 1:
                            raise ValueError("Cannot proceed with different counting_features..")

                        # if len(c_library_prep) > 1:
                        #     raise ValueError("Cannot proceed with different library_preps..")
        
    except Exception as e:
        log.error(str(e)+"\n"+sample_sheet_info)
        sys.exit(e)

    return sample_sheet_dict