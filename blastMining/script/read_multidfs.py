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
import glob

def read_multidfs(OUTPUT, EXTENSION):
    OUTPUT=OUTPUT
    EXTENSION=EXTENSION
    def df_generator(OUTPUT, EXTENSION):
        for file in glob.glob(OUTPUT+'/' + EXTENSION):
            aux=pd.read_csv(file, sep='\t')
            yield aux

    generator=df_generator(OUTPUT, EXTENSION)
    df_list = []
    for table in generator:
        df_list.append(table)
        
    DF = pd.concat(df_list)
    return DF