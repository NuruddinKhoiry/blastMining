#!/usr/bin/env python3
"""
Copyright 2022 Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

https://github.com/NuruddinKhoiry/blastMining
This file is a part of blastMining. blastMining is a free software: you can redistribute it and/or modify
it under the terms of GNU General Public License v3.0. blastMining is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""

import pandas as pd

def summary2krona(table_in, table_out):
    df = pd.read_csv(table_in, delimiter='\t')
    
    col = df.columns.to_list()[::-1]
    df2 = df[col]
    df3 = df2.Taxa.apply(lambda x: pd.Series(str(x).split(";")))
    
    df4 = pd.concat([df2, df3], axis=1)
    col += ['Kingdom','Phylum','Class','Order','Family','Genus','Species']
    df4.columns = col
    
    df4.drop('Taxa',  axis=1, inplace=True)
    
    df4.to_csv(table_out, header=True, index=None, sep='\t')