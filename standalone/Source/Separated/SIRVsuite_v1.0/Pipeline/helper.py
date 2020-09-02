import re
import csv
import sys

def path_features(filePath):
    # function returning features of a given path
    matching_pattern = "((?:(?:.*\\/)*(.*)\\/))*(?:(.*)\\.(.*))*"
    feature_match = re.match(matching_pattern, filePath)
    out_dict = {"parent_dir": feature_match.group(2),
                "file_name": feature_match.group(3),
                "extension": feature_match.group(4),
                "path": feature_match.group(1)}
    return out_dict
    
def __create_nested_dict__(input_list, init = 0):
        """
        Helper method using a recursive approach to create and initialize nested dictionary based on 2D-list of values.
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
        if valid_dict[colname]:
            if value.lower() not in valid_dict[colname]:
                raise ValueError("Invalid entry %s for %s column.. possible options are: %s"%(value, colname, ", ".join(valid_dict[colname])))
        else:
            raise ValueError("%s column not defined.. please define it and insert one of the possible values: %s"%(colname,", ".join(valid_dict[colname])))
    else:
        raise ValueError("colname not in valid_dict")

def read_sample_sheet(sheet_path):

    sample_sheet_dict = {}
    
    # name required columns
    required_cols = dict()
    required_cols["concentration"] = ["sample_name","counting_path","counting_method","counting_feature","library_prep_type"]
    required_cols["coverage"] = ["sample_name","alignment_path","read_orientation","library_prep_type"]
    
    # name optional columns
    optional_cols = dict()
    optional_cols["concentration"] = ["replicate_group"]
    optional_cols["coverage"] = []

    # create a set of all possible columns in a sample sheet
    required_all = set()
    for mode in required_cols.keys():
        required_all = required_all.union(set(required_cols[mode]))

    # create dict from restricted values for a specified column
    value_restriction = {col_name:[] for col_name in required_all}
    value_restriction["library_prep_type"] = ["whole","qs"]
    value_restriction["read_orientation"] = ["fwd","rev","none"]
    value_restriction["counting_type"] = ["mix2","cufflinks","htseq"]
    value_restriction["counting_feature"] = ["gene","trasncript"]
    
    try:
        with open(sheet_path, 'r') as sheet:

            csv.register_dialect('strip', skipinitialspace=True, delimiter = ';')

            header = [h.strip() for h in next(sheet).split(";")]

            # check if cols are complete for specific modes 
            available_modes = [mode for mode in required_cols.keys() if set(required_cols[mode]) <= set(header)]
            
            reader = csv.DictReader(sheet, restkey=None, restval=None, dialect='strip', fieldnames=header)

            for row in reader:

                # remove tabs
                row = {a:row[a].strip("\t") for a in row.keys()}

                if None in row.values():
                    raise ValueError("Incorrect number of rows in a sample sheet..")

                __entry_check__(value_restriction, "library_prep_type", row.get("library_prep_type"))
                
                # check for valid values for a given column 
                if 'coverage' in available_modes:
                    __entry_check__(value_restriction, "library_prep_type", row.get("library_prep_type"))
                    sample_sheet_dict["coverage"] = {row.get("sample_name"):{val:row[val] for val in required_cols["coverage"][1:] if val in row}}

                if 'concentration' in available_modes:
                    __entry_check__(value_restriction, "counting_type", row.get("library_prep_type"))
                    __entry_check__(value_restriction, "counting_feature", row.get("counting_feature"))

                    sample_sheet_dict["concentration"] = {row.get("sample_name"):{val:row[val] for val in required_cols["concentration"][1:] if val in row}}
        
        
    except ValueError as e:
        sys.exit(e)

    return sample_sheet_dict