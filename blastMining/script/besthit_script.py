#!/usr/bin/env python3

import numpy as np
import pandas as pd
from fastnumbers import fast_forceint


def besthit(blast, tax, pident, evalue,topN):
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
                
            tax_df = pd.merge(df_eval, tax, left_index=True, right_index=True)
            taxa = tax_df.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
            taxa = taxa.reset_index(drop=True)
            final_df.append(taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']])
    DF = pd.concat(final_df)
    DF = DF.reset_index(inplace=False, drop=True)
    return(DF)