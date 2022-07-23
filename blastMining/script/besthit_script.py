#!/usr/bin/env python3
"""
Copyright 2022 Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

https://github.com/NuruddinKhoiry/blastMining
This file is a part of blastMining. blastMining is a free software: you can redistribute it and/or modify
it under the terms of GNU General Public License v3.0. blastMining is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""

import numpy as np
import pandas as pd
from fastnumbers import fast_forceint

def besthit(blast, pident, evalue,topN):
    
    blast[['pident', 'bitscore', 'evalue', 'mismatch']] = blast[['pident', 'bitscore', 'evalue', 'mismatch']].apply(pd.to_numeric)
    
    final_df = []
    for i, sub_df in blast.groupby(blast.qseqid.ne(blast.qseqid.shift()).cumsum()):
        df_eval = sub_df.loc[(sub_df.evalue <= evalue) & (sub_df.pident >= pident)]
        if df_eval.shape[0] < 1:
            sub_df.reset_index(inplace=True, drop=True)
            df_qseqid = sub_df.iloc[:1,:1]
            df_qseqid[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = pd.DataFrame([['k__',
                'p__','c__','o__','f__','g__', 's__']], index=df_qseqid.index)
            final_df.append(df_qseqid)
        else:
            if df_eval.shape[0] > topN:
                df_eval = df_eval.iloc[:topN,:]
            else:
                df_eval = df_eval
                
            taxa = df_eval.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
            taxa = taxa.reset_index(drop=True)
            final_df.append(taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']])
            
    DF = pd.concat(final_df)
    DF = DF.reset_index(inplace=False, drop=True)
    
    return(DF)