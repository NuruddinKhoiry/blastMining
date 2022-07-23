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

def summary_df(DF, sample_name):
    DF['Taxa'] = DF[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']].apply(';'.join, axis=1)
    DT = DF.drop(columns=['Kingdom','Phylum','Class','Order','Family','Genus', 'Species'])
    summary = pd.DataFrame(DT.groupby(['Taxa'])['Taxa'].count())
    summary.columns = [str(sample_name)]
    summary_DF = summary.reset_index(level=0)
    return(summary_DF)