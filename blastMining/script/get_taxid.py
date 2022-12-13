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

def get_taxid(df, out):

    df['Kingdom'] = df['Kingdom'].str.replace('k__','')
    df['Phylum'] = df['Phylum'].str.replace('p__','')
    df['Class'] = df['Class'].str.replace('c__','')
    df['Order'] = df['Order'].str.replace('o__','')
    df['Family'] = df['Family'].str.replace('f__','')
    df['Genus'] = df['Genus'].str.replace('g__','')
    df['Species'] = df['Species'].str.replace('s__','')

    for col in range(1, len(df.columns)):
        for row in range(0, df.shape[0]):
            if df.iat[row,1] == '':
                df.iat[row,1] = 'Root' 
            if df.iat[row,col] == '':
                if df.iat[row,col-1] != '':
                    df.iat[row,col] = df.iat[row,col-1]
                else:
                    df.iat[row,col] = df.iat[row,col-1]

    df.to_csv(out, header=False, index=None, sep='\t') 

