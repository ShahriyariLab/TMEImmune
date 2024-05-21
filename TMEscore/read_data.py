import pandas as pd
import os
from cmapPy.pandasGEXpress.parse_gct import parse

def read_data(path):

    # gct, txt, csv, first column as genes
    file_name, file_extension = os.path.splitext(path)
    if file_extension == '.txt':
        df = pd.read_table(path, sep = "\t") 
        genes = df.iloc[:, 0]
    elif file_extension == '.csv':
        df = pd.read_csv(path)
        genes = df.iloc[:, 0]
    elif file_extension == '.gct':
        df = parse(path)
        df = df.data_df
        genes = df.index
        df.insert(loc=0, column='gene', value=genes)
    else:
        raise TypeError(file_extension + " file is not supported")

    if genes.duplicated().any():
        d_ind = genes[genes.duplicated()].index
        if len(d_ind) == 1:
            d_ind = d_ind[0]+1
        else:
            d_ind = tuple([d_ind[i] +1 for i in range(len(d_ind))])
        raise ValueError("Duplicate genes in row " + str(d_ind))
    
    if genes.isnull().any():
        raise ValueError("Exist NA's in the input genes")

    # normalization? log2, min-max, z-score

    return (df)

