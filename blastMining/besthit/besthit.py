#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
from fastnumbers import fast_forceint
from blastMining.script import summary_df
from blastMining.script import besthit_script
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
        help="Top N hits used for sorting [default=10]")
    parser.add_argument("-sm","--sample_name",dest="sample_name",action="store",default='sample',type=str,
        help='Sample name in the print out table [default="sample"]')
    parser.add_argument("-kp","--krona_plot",dest="krona_plot", default=argparse.SUPPRESS, action='store_true',
        help='Draw krona plot')
    parser.add_argument("-o", "--output", type=str, required=True, 
        help="output")
    parser.add_argument('-v','--version', action='version', version='blastMining v.1.0.0')
    
    return parser

def main(args):
 
    os.system("cat "+str(args.input)+" | cut -f 9 | taxonkit lineage | taxonkit reformat -P | csvtk -H -t cut -f 1,3 > output.tmp")
    
    blast = pd.read_csv(args.input, sep='\t', header=None, dtype=str)
    blast.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid"]
    
    tax = pd.read_csv('output.tmp', sep='\t', header=None, dtype=str)
    tax.columns = ['staxid', 'lineage']
    tax[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = tax['lineage'].str.split(';', expand=True)
    tax.drop('lineage', axis=1, inplace=True)
    
    DF = besthit_script.besthit(blast=blast, tax=tax, pident=args.pident, evalue=args.evalue, topN=args.topN)
    DF.to_csv(str(args.output+'.tsv'), header=True, index=None, sep='\t')
    
    SD = summary_df.summary_df(DF, args.sample_name)
    SD.to_csv(str(args.output+'.summary'), header=True, index=None, sep='\t')
    
    if hasattr(args, 'krona_plot'):
        summary2krona.summary2krona(str(args.output+'.summary'), str(args.output+'.krona'))
        os.system("ktImportText "+str(args.output+'.krona')+' -o '+str(args.output+'.html'))
    
if __name__ == '__main__':
    main()