import os
import sys
import numpy as np
import pandas as pd
from fastnumbers import fast_forceint
import argparse


def add_arguments(parser):
    parser.add_argument("-i", "--input", type=str, required=True, 
        help='''blast.out file. Please use this blast outfmt 6 ONLY:
               ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")''')
    parser.add_argument("-e","--evalue",dest="evalue",action="store",default=1e-3,type=float,
        help='''Threshold of evalue 
               (Ignore hits if their evalues are above this threshold)
               [default=1-e3]''')
    parser.add_argument("-n","--topN", dest="topN",action="store",default=10,type=int,
        help="Top N hits used for voting [default=10]")
    parser.add_argument("-txl","--taxa_level", dest="taxa_level",action="store",default=[99,97,95,90,85,80,75],type=list,
        help='''P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species
                [default=99,97,95,90,85,80,75]''')
    parser.add_argument("-o", "--output", type=str, required=True, 
        help="output")  
    return parser

def average(lst):
    return sum(lst) / len(lst)

def vote(blast, tax, topN, evalue, tax_level):
    final_df = []
    for i, sub_df in blast.groupby(blast.qseqid.ne(blast.qseqid.shift()).cumsum()):
        df_sort = sub_df.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
        df_toN = df_sort.iloc[:topN,:]
        if df_toN.shape[0] >= topN:
            df_eval = df_toN.loc[df_toN.evalue.astype(float) <= evalue, :]
            pident = [fast_forceint(x) for x in df_eval.pident.tolist()]
            df_filt = df_eval.loc[df_eval.pident.astype(float) >= average(pident)]
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
                    
            elif max(pident) >= tax_level[1] & max(pident) < tax_level[0]:
                tx, nm = np.unique(tax_df[['Genus']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    genus = tax_df.loc[tax_df['Genus'] == taxa][['qseqid',
                            'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']].drop_duplicates()
                    genus[['Species']] = 's__'
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

            elif max(pident) >= tax_level[2]  & max(pident) < tax_level[1]:
                tx, nm = np.unique(tax_df[['Family']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    family = tax_df.loc[tax_df['Family'] == taxa][['qseqid',
                            'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']].drop_duplicates()
                    family[['Genus','Species']] = ['g__', 's__']
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

            elif max(pident) >= tax_level[3]  & max(pident) < tax_level[2]:
                tx, nm = np.unique(tax_df[['Order']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    order = tax_df.loc[tax_df['Order'] == taxa][['qseqid',
                            'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']].drop_duplicates()
                    order[['Family','Genus','Species']] = ['f__','g__', 's__']
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
                    
            elif max(pident) >= tax_level[4] & max(pident) < tax_level[3]:
                tx, nm = np.unique(tax_df[['Class']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    clas = tax_df.loc[tax_df['Class'] == taxa][['qseqid',
                            'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']].drop_duplicates()
                    clas[['Order','Family','Genus','Species']] = ['o__','f__','g__', 's__']
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
                    
            elif max(pident) >= tax_level[5] & max(pident) < tax_level[4]:
                tx, nm = np.unique(tax_df[['Phylum']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    phylum = tax_df.loc[tax_df['Phylum'] == taxa][['qseqid',
                            'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']].drop_duplicates()
                    phylum[['Class','Order','Family','Genus','Species']] = ['c__','o__','f__','g__', 's__']
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

            elif max(pident) >= tax_level[6] & max(pident) < tax_level[5]:
                tx, nm = np.unique(tax_df[['Kingdom']], return_counts=True)
                if nm.tolist().count(max(nm)) == 1:
                    taxa = tx.tolist()[nm.tolist().index(max(nm.tolist()))]
                    kingdom = tax_df.loc[tax_df['Kingdom'] == taxa][['qseqid',
                            'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']].drop_duplicates()
                    kingdom[['Phylum','Class','Order','Family','Genus','Species']] = ['p__','c__','o__','f__','g__', 's__']
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
                taxa = tab_df.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
                unasg = taxa.loc[:0,:][['qseqid', 'Kingdom','Phylum','Class','Order','Family','Genus', 'Species']]
                unasg[['Kingdom','Phylum','Class','Order','Family','Genus','Species']] = ['k__','p__','c__','o__','f__','g__', 's__']
                final_df.append(unasg)   

        else:
            df_toN2 = df_toN.assign(Kingdom=['k__'],Phylum=['p__'],
                          Class=['c__'],Order=['o__'],Family=['f__'],Genus=['g__'],Species=['s__'])
            final_df.append(df_toN2[['qseqid','Kingdom','Phylum','Class','Order','Family','Genus','Species']])
    return(pd.concat(final_df))

def main(args):
    os.system("cat "+str(args.input)+" | cut -f 9 | taxonkit lineage | taxonkit reformat -P | csvtk -H -t cut -f 1,3 > output.tmp")
    
    blast = pd.read_csv(args.input, sep='\t', header=None, dtype=str)
    blast.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid"]
    
    tax = pd.read_csv('output.tmp', sep='\t', header=None, dtype=str)
    tax.columns = ['staxid', 'lineage']
    tax[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = tax['lineage'].str.split(';', expand=True)
    tax.drop('lineage', axis=1, inplace=True)
    
    taxa_level = [fast_forceint(x) for x in args.taxa_level]
    
    DF = vote(blast=blast, tax=tax, evalue=args.evalue, topN=args.topN, tax_level=taxa_level)
    DF.to_csv(args.output, header=True, index=None, sep='\t')
    
    os.system("rm output.tmp")
    

if __name__ == '__main__':
    main()