#!/usr/bin/env python3

import numpy as np
import pandas as pd
from fastnumbers import fast_forceint

def vote(blast, tax, evalue, tax_level, topN):
    tax_level = [fast_forceint(x) for x in tax_level]
    blast[['pident', 'bitscore', 'evalue', 'mismatch']] = blast[['pident', 'bitscore', 'evalue', 'mismatch']].apply(pd.to_numeric)
    final_df = []
    for i, sub_df in blast.groupby(blast.qseqid.ne(blast.qseqid.shift()).cumsum()):
        df_sort = sub_df.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
        df_eval = df_sort.loc[df_sort.evalue.astype(float) <= evalue, : ]
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
                
            pident = [fast_forceint(x) for x in df_eval.pident.tolist()]
            df_filt = df_eval.loc[df_eval.pident.astype(float) >= sum(pident)/len(pident),: ]
            tax_df = pd.merge(df_filt, tax, left_index=True, right_index=True)
            if max(pident) >= tax_level[0]:
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

            elif max(pident) >= tax_level[1] and max(pident) < tax_level[0]:
                tx, nm = np.unique(tax_df[['Genus']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    genus = tax_df.loc[tax_df['Genus'] == taxa][['qseqid','Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    genus[['Species']] = 's__'
                    genus = genus.drop_duplicates()
                    final_df.append(genus)
                else:
                    indices = [index for index, value in enumerate(nm.tolist()) if value == max(nm.tolist())]
                    selected = [tx.tolist()[i] for i in indices]
                    df_selected = tax_df[tax_df.Species.isin(selected)]
                    taxa = df_selected.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
                    taxa = taxa.reset_index(drop=True)
                    genus = taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    genus[['Species']] = 's__'
                    final_df.append(genus)

            elif max(pident) >= tax_level[2] and max(pident) < tax_level[1]:
                tx, nm = np.unique(tax_df[['Family']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    family = tax_df.loc[tax_df['Family'] == taxa][['qseqid','Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    family[['Genus','Species']] = ['g__', 's__']
                    family = family.drop_duplicates()
                    final_df.append(family)
                else:
                    indices = [index for index, value in enumerate(nm.tolist()) if value == max(nm.tolist())]
                    selected = [tx.tolist()[i] for i in indices]
                    df_selected = tax_df[tax_df.Species.isin(selected)]
                    taxa = df_selected.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
                    taxa = taxa.reset_index(drop=True)
                    family = taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    family[['Genus','Species']] = ['g__', 's__']
                    final_df.append(family)

            elif max(pident) >= tax_level[3] and max(pident) < tax_level[2]:
                tx, nm = np.unique(tax_df[['Order']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    order = tax_df.loc[tax_df['Order'] == taxa][['qseqid','Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    order[['Family','Genus','Species']] = ['f__','g__', 's__']
                    order = order.drop_duplicates()
                    final_df.append(order)
                else:
                    indices = [index for index, value in enumerate(nm.tolist()) if value == max(nm.tolist())]
                    selected = [tx.tolist()[i] for i in indices]
                    df_selected = tax_df[tax_df.Species.isin(selected)]
                    taxa = df_selected.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
                    taxa = taxa.reset_index(drop=True)
                    order = taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    order[['Family','Genus','Species']] = ['f__','g__', 's__']
                    final_df.append(order)

            elif max(pident) >= tax_level[4] and max(pident) < tax_level[3]:
                tx, nm = np.unique(tax_df[['Class']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    clas = tax_df.loc[tax_df['Class'] == taxa][['qseqid','Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    clas[['Order','Family','Genus','Species']] = ['o__','f__','g__', 's__']
                    clas = clas.drop_duplicates()
                    final_df.append(clas)
                else:
                    indices = [index for index, value in enumerate(nm.tolist()) if value == max(nm.tolist())]
                    selected = [tx.tolist()[i] for i in indices]
                    df_selected = tax_df[tax_df.Species.isin(selected)]
                    taxa = df_selected.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
                    taxa = taxa.reset_index(drop=True)
                    clas = taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    clas[['Order','Family','Genus','Species']] = ['o__','f__','g__', 's__']
                    final_df.append(clas)

            elif max(pident) >= tax_level[5] and max(pident) < tax_level[4]:
                tx, nm = np.unique(tax_df[['Phylum']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    phylum = tax_df.loc[tax_df['Phylum'] == taxa][['qseqid','Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    phylum[['Class','Order','Family','Genus','Species']] = ['c__','o__','f__','g__', 's__']
                    phylum = phylum.drop_duplicates()
                    final_df.append(phylum)
                else:
                    indices = [index for index, value in enumerate(nm.tolist()) if value == max(nm.tolist())]
                    selected = [tx.tolist()[i] for i in indices]
                    df_selected = tax_df[tax_df.Species.isin(selected)]
                    taxa = df_selected.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
                    taxa = taxa.reset_index(drop=True)
                    phylum = taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    phylum[['Class','Order','Family','Genus','Species']] = ['c__','o__','f__','g__', 's__']
                    final_df.append(phylum)

            elif max(pident) >= tax_level[6] and max(pident) < tax_level[5]:
                tx, nm = np.unique(tax_df[['Kingdom']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    kingdom = tax_df.loc[tax_df['Kingdom'] == taxa][['qseqid','Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    kingdom[['Phylum','Class','Order','Family','Genus','Species']] = ['p__','c__','o__','f__','g__', 's__']
                    kingdom = kingdom.drop_duplicates()
                    final_df.append(kingdom)
                else:
                    indices = [index for index, value in enumerate(nm.tolist()) if value == max(nm.tolist())]
                    selected = [tx.tolist()[i] for i in indices]
                    df_selected = tax_df[tax_df.Species.isin(selected)]
                    taxa = df_selected.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
                    taxa = taxa.reset_index(drop=True)
                    kingdom = taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                    kingdom[['Phylum','Class','Order','Family','Genus','Species']] = ['p__','c__','o__','f__','g__', 's__']
                    final_df.append(kingdom)

            else:
                taxa = tax_df.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
                taxa = taxa.reset_index(drop=True)
                unasg = taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                unasg[['Kingdom','Phylum','Class','Order','Family','Genus','Species']] = ['k__','p__','c__','o__','f__','g__', 's__']
                final_df.append(unasg)
    DF = pd.concat(final_df)
    DF = DF.reset_index(inplace=False, drop=True)
    return(DF)