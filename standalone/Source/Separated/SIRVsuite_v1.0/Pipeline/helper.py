import re
import csv
import sys

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
    
def __create_nested_dict__(input_list, init = 0):
        """
        Helper function using a recursive approach to create and initialize nested dictionary based on 2D-list of values.
        Length of input list defines the depth of a dictionary, the nested lists defines the keys at a corresponding level.
        Values are initialized based on init value. 
        
        Example: [[key1_level1,key2_level1,key3_level1],[[key1_level2,key2_level2]] as input_list creates a list corresponding to 
        
        {
         key1_level1: {key1_level2: 0, key2_level2: 0},
         key2_level1: {key1_level2: 0, key2_level2: 0},
         key3_level1: {key1_level2: 0, key2_level2: 0}
        }
        
        """

        input = input_list[0]
        input_list.pop(0)
        l = len(input_list)

        if (l > 0):
            temp = __create_nested_dict__(input_list)
            temp = [copy.deepcopy(temp) for i in range(len(input))]
            return dict(zip(input, temp))

        elif (l == 0):
            n = len(input)
            init_list = [init] * n
            return dict(zip(input, init_list))

def __entry_check__(valid_dict, colname, value):
    if colname in valid_dict.keys():
        if value:
            if value.lower() not in valid_dict[colname]:
                raise ValueError("Invalid entry %s for %s column.. possible options are: %s"%(value, colname, ", ".join(valid_dict[colname])))
            else:
                return value.lower()
        else:
            raise ValueError("%s column not defined.. please define it and insert one of the possible values: %s"%(colname,", ".join(valid_dict[colname])))
    else:
        raise ValueError("colname not in valid_dict")

def read_sample_sheet(sheet_path, modules_to_execute = ["concentration", "coverage"]):
    """
    A function for loading sample sheet in a following format:
        - used ";" separator
        - trailing whitespace or tab is allowed (trimmed during reading process)
        - UTF-8 encoding
        - every column must have a name which is predefined (see more info below)

    Input: path to the sample sheet file

    Allowed column names:
        General columns:
            sample_name - Supported value: any set of characters to identify samples (this will be printed in the final graphics).
            
            library_prep_type - supported values: whole (whole transcriptome library prep) or qs (QuantSeq library preparation) 

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

    ## ALLOWED COLUMNS FOR DIFFERENT MODULES CAN BE DEFINED HERE ##

    # name required columns
    required_cols = dict()
    required_cols["concentration"] = ["sample_name","counting_path","counting_method","counting_feature","library_prep_type"]
    required_cols["coverage"] = ["sample_name","alignment_path","read_orientation","library_prep_type"]
    
    # name optional columns
    optional_cols = dict()
    optional_cols["concentration"] = ["replicate_group"]
    optional_cols["coverage"] = []

    ## -------------------------------------------------------- ##

    ## RESTRICTIONS FOR THE COLUMNS CAN BE DEFINED HERE ##
    
    value_restriction = dict()
    value_restriction["library_prep_type"] = ["whole","qs"]
    value_restriction["read_orientation"] = ["fwd","rev","none"]
    value_restriction["counting_method"] = ["mix2","cufflinks","htseq"]
    value_restriction["counting_feature"] = ["gene","transcript"]

    ## ------------------------------------------------ ##

    sample_sheet_dict = dict()
  
    # select only columns for which a restriction has been defined
    restricted_cols = dict()
    for module in required_cols:
        restricted_cols[module] = [col for col in value_restriction.keys() if col in required_cols[module]]
    
    # select columns with defined restriction for all common modules
    common_restricted_cols = set()
    modules = list(required_cols.keys())
    for idx in range(len(modules)-1):
        common_restricted_cols = set(required_cols[modules[idx]]).intersection(set(required_cols[modules[idx+1]]))
    common_restricted_cols = [col for col in common_restricted_cols if col in value_restriction.keys()]

    try:
        with open(sheet_path, 'r') as sheet:

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
                # remove tabs
                row = {a:row[a].strip("\t ") for a in row.keys()}

                if None in row.values():
                    raise ValueError("Incorrect number of rows in a sample sheet..")
                
                # check for valid values for general columns (all modules)
                for common_col in common_restricted_cols:
                    __entry_check__(value_restriction, common_col, row.get(common_col))
                
                # check for valid values for a given column across modules
                for module in modules_to_execute:
                    if module in available_modules:
                        for col in restricted_cols[module]:
                            row[col] = __entry_check__(value_restriction, col, row.get(col))
                            sample_sheet_dict[module][row.get("sample_name")] = {val:row[val] for val in required_cols[module][1:] if val in row}

                        for col in optional_cols[module]:
                            if col in header:
                                sample_sheet_dict[module][row.get("sample_name")][col] = row.get(col)  
                    else:
                        raise ValueError("Cannot proceed with module: %s"%(module))

        c_methods = set([sample_sheet_dict["concentration"][i]["counting_method"] for i in sample_sheet_dict["concentration"].keys()])
        c_features = set([sample_sheet_dict["concentration"][i]["counting_feature"] for i in sample_sheet_dict["concentration"].keys()])
        c_library_prep = set([sample_sheet_dict["concentration"][i]["library_prep_type"] for i in sample_sheet_dict["concentration"].keys()])
        
        if len(c_methods) > 1:
            raise ValueError("Cannot proceed with different counting_methods..")

        if len(c_features) > 1:
            raise ValueError("Cannot proceed with different counting_features..")

        if len(c_library_prep) > 1:
            raise ValueError("Cannot proceed with different library_preps..")
        
    except ValueError as e:
        sys.exit(e)

    return sample_sheet_dict