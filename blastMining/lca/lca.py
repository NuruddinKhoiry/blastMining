#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
from fastnumbers import fast_forceint
from blastMining.script import summary_df
from blastMining.script import lca_script
from blastMining.script import summary2krona
import argparse

def add_arguments(parser):

    parser.add_argument("-i", "--input", type=str, required=True, 
        help='''blast.out file. Please use this blast outfmt 6 ONLY:
("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")''')
    parser.add_argument("-e","--evalue",dest="evalue",action="store",default=1e-3,type=float,
        help='''Threshold of evalue
(Ignore hits if their evalues are above this threshold)
[default=1-e3]''')
    parser.add_argument("-pi","--pident", dest="pident",action="store",default=97,type=int,
        help='''Threshold of p. identity
(Ignore hits if their p. identities are below this threshold)
[default=97]''')
    parser.add_argument("-n","--topN", dest="topN",action="store",default=10,type=int,
        help="Top N hits used for LCA calculation [default=10]")
    parser.add_argument("-sm","--sample_name",dest="sample_name",action="store",default='sample',type=str,
        help='Sample name in the print out table [default="sample"]')
    parser.add_argument("-kp","--krona_plot",dest="krona_plot", default=argparse.SUPPRESS,action='store_true',
        help='Draw krona plot')
    parser.add_argument("-o", "--output", type=str, required=True, 
        help="output") 
    parser.add_argument('-v','--version', action='version', version='blastMining v.1.0.0')
    
    return parser

def main(args):

    if len(sys.argv) == 2:
        os.system('blastMining lca -h')
        sys.exit()
        
    blast = pd.read_csv(args.input, sep='\t', header=None, dtype=str)
    blast.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid"]
    
    DF = lca_script.lca(blast=blast, evalue=args.evalue, pident=args.pident, topN=args.topN)
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
    
    DT2.to_csv(str(args.output+'.tsv'), header=True, index=None, sep='\t')
    
    SD = summary_df.summary_df(DT2, args.sample_name)
    SD.to_csv(str(args.output+'.summary'), header=True, index=None, sep='\t')
    
    os.system("rm tmp_lca")
    os.system("rm tmp_lca2")
    os.system("rm LCA")

    if hasattr(args, 'krona_plot'):
        summary2krona.summary2krona(str(args.output+'.summary'), str(args.output+'.krona'))
        os.system("ktImportText "+str(args.output+'.krona')+' -o '+str(args.output+'.html'))    

if __name__ == '__main__':
    main()