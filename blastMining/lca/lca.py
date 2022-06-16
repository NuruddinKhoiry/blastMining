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

def lca(blast,topN,evalue):
    final_df = []
    for i, sub_df in blast.groupby(blast.qseqid.ne(blast.qseqid.shift()).cumsum()):
        df_sort = sub_df.sort_values(['pident', 'bitscore', 'evalue', 'mismatch'], ascending=[False, False, True, True])
        df_toN = df_sort.iloc[:topN,:]
        if df_toN.shape[0] >= topN:
            df_eval = df_toN.loc[df_toN.evalue.astype(float) <= evalue, :]
            pident = [fast_forceint(x) for x in df_eval.pident.tolist()]
            df_filt = df_eval.loc[df_eval.pident.astype(float) >= average(pident)]
            df_gb = df_filt.groupby('qseqid').agg(lambda x: ','.join(set(x)))
            df_gb.reset_index(inplace=True)
            final_df.append(df_gb[['qseqid','staxid']])
    return(pd.concat(final_df))


def main(args):
    
    blast = pd.read_csv(args.input, sep='\t', header=None, dtype=str)
    blast.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid"]
    
    DF = lca(blast=blast, evalue=args.evalue, topN=args.topN)
    DF.to_csv('LCA', header=None, index=None, sep='\t')
    
    os.system("bash blastMining_lca.sh LCA")
    os.system("bash blastMining_lca2.sh tmp_lca")
    
    lca1 = pd.read_csv('tmp_lca', sep='\t', header=None, dtype=str)
    lca1.columns = ['staxid', 'lca']
    lca2 = pd.read_csv('tmp_lca2', sep='\t', header=None, dtype=str)
    lca2.columns = ['lca', 'lineage']
    
    DT2 = DF.merge(lca1).merge(lca2)
    DT2[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = DT2['lineage'].str.split(';', expand=True)
    DT2.drop(['staxid','lineage'], axis=1, inplace=True)
    
    DT2.to_csv(args.output, header=True, index=None, sep='\t')
    
    os.system("rm tmp_lca")
    os.system("rm tmp_lca2")
    os.system("rm LCA")  

if __name__ == '__main__':
    main()