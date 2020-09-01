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


def read_sample_sheet(sheet_path):
   
    required_colnames = ["sample_name","counting_path","alignment_path","read_orientation"]
    
    try:
        with open(sheet_path) as sheet_path:
            csv.register_dialect('strip', skipinitialspace=True)
            reader = csv.DictReader(sheet_path, restkey=None, restval=None, delimiter = ';', dialect='strip')

            for row in reader:
                if not set(required_colnames) == row.keys():
                    raise ValueError("")
                if None in row.values():
                    raise ValueError("")
                
    except ValueError as e:
        sys.exit(e)

    sample_sheet_dict = None ## DELETE

    return sample_sheet_dict
