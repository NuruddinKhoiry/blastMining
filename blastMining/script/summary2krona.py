#!/usr/bin/env python3

import pandas as pd

def summary2krona(table_in, table_out):
    df = pd.read_csv(table_in, delimiter='\t')
    
    col = df.columns.to_list()[::-1]
    df2 = df[col]
    df3 = df2.Taxa.apply(lambda x: pd.Series(str(x).split(";")))
    
    df4 = pd.concat([df2, df3], axis=1)
    col += ['K','P','C','O','F','G','S']
    df4.columns = col
    
    df4.drop('Taxa',  axis=1, inplace=True)
    
    df4.to_csv(table_out, header=True, index=None, sep='\t')