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

def voteSpecies(blast, pident, evalue, topN):
    blast[['pident', 'bitscore', 'evalue', 'mismatch']] = blast[['pident', 'bitscore', 'evalue', 'mismatch']].apply(pd.to_numeric)
    final_df = []
    for i, sub_df in blast.groupby(blast.qseqid.ne(blast.qseqid.shift()).cumsum()):
        df_sort = sub_df.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
        df_eval = df_sort.loc[(df_sort.evalue <= evalue) & (df_sort.pident >= pident)]
        if df_eval.shape[0] < 1:
            df_sort.reset_index(inplace=True, drop=True)
            df_qseqid = df_sort.iloc[:1,:1]
            df_qseqid[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = pd.DataFrame([['k__',
                'p__','c__','o__','f__','g__', 's__']], index=df_qseqid.index)
            final_df.append(df_qseqid)
        else:
            if df_eval.shape[0] > topN:
                df_eval = df_eval.iloc[:topN,:]
            else:
                df_eval = df_eval
                
            pidentity = [fast_forceint(x) for x in df_eval.pident.tolist()]
            tax_df = df_eval.loc[df_eval.pident.astype(float) >= sum(pidentity)/len(pidentity),: ]
            tx, nm = np.unique(tax_df[['Species']], return_counts=True)
            if nm.tolist().count(max(nm)) == 1:
                taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                final_df.append(tax_df.loc[tax_df['Species'] == taxa][['qseqid',
                            'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']].drop_duplicates())
            else:
                indices = [index for index, value in enumerate(nm.tolist()) if value == max(nm.tolist())]
                selected = [tx.tolist()[i] for i in indices]
                df_selected = tax_df[tax_df.Species.isin(selected)]
                taxa = df_selected.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
                taxa = taxa.reset_index(drop=True)
                final_df.append(taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']])
                
    DF = pd.concat(final_df)
    DF = DF.reset_index(inplace=False, drop=True)
    return(DF)