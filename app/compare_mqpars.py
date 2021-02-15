import pandas as pd
import collections
import xmltodict

#result_1 = parse_mqpar("example_par_files/UPS")
#result_2 = parse_mqpar("example_par_files/iPRG2015")
#result_3 = parse_mqpar("example_par_files/HER")

def parse_mqpar(folder):
    with open(os.path.join(folder, "mqpar.xml"), 'r') as file:
        data = file.read().replace('\n', '')
    result = xmltodict.parse(data)
    return flatten(result['MaxQuantParams'])

def compare_mqpars(result_1, result_2):
    agreement = (pd.Series(result_1) == pd.Series(result_2)) | (pd.Series(result_1).isna() & pd.Series(result_2).isna())
    disagree = agreement.index[~agreement]

    irrelevant = ['fastaFiles_FastaFileInfo_fastaFilePath', 
                    'filePaths_string', 'experiments_string', 
                    'fractions_short', 'ptms_boolean', 
                    'paramGroupIndices_int', 
                    'referenceChannel_string']
    
    meaningful_disagreement = [i for i in disagree if i not in irrelevant]

    for i in meaningful_disagreement:
        print(i)
        print("Result 1: {}".format(result_1[i]))
        print("Result 2: {}".format(result_2[i]))
        print("")

    return #meaningful_disagreement

# https://stackoverflow.com/questions/6027558/flatten-nested-dictionaries-compressing-keys

def flatten(d, parent_key='', sep='_'):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)