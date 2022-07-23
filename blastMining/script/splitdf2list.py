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
import numpy as np

def splitdf2list(df, column, chunks):
    df[column] = df[column].astype('category')
    
    list_seq = df[column].unique().tolist()
    
    k, m = divmod(len(list_seq), chunks)
    x = list(list_seq[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(chunks))
    
    if chunks > len(df[column].unique()):
        chunksize = len(df[column].unique())
    else:
        chunksize = chunks
    
    df_split = []
    for i in range(chunksize):
        df_split.append(df[df[column].isin(x[i])])
    
    return df_split