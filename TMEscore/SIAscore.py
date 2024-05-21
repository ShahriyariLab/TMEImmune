import pandas as pd
import numpy as np

def sia_score(df, marker = "C1QA"): # C1QA, C1QB, C1QC
    genes = df.iloc[:, 0]
    df.index = genes
    df = df.iloc[:,1:]
    sig = ["CD8A", marker]
    if not all(i in list(genes) for i in sig):
        raise ValueError("Genes needed to get the SIA score do not exist in the input data")

    sia_score = df.loc['CD8A', :]/df.loc[marker, :]
    return sia_score
