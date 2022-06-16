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
    parser.add_argument("-o", "--output", type=str, required=True, 
        help="output")  
    return parser

def average(lst):
    return sum(lst) / len(lst)

def average(lst):
    return sum(lst) / len(lst)

def voteSpecies(blast, tax, topN, evalue):
    final_df = []
    for i, sub_df in blast.groupby(blast.qseqid.ne(blast.qseqid.shift()).cumsum()):
        df_sort = sub_df.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
        df_toN = df_sort.iloc[:topN,:]
        if df_toN.shape[0] >= topN:
            df_eval = df_toN.loc[df_toN.evalue.astype(float) <= evalue, :]
            pident = [fast_forceint(x) for x in df_eval.pident.tolist()]
            df_filt = df_eval.loc[df_eval.pident.astype(float) >= average(pident)]
            tax_df = pd.merge(df_filt, tax, left_index=True, right_index=True)
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
    
    DF = voteSpecies(blast=blast, tax=tax, evalue=args.evalue, topN=args.topN)
    DF.to_csv(args.output, header=True, index=None, sep='\t')
    
    os.system("rm output.tmp")

if __name__ == '__main__':
    main()