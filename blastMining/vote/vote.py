#!/usr/bin/env python3
"""
Copyright 2022 Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

https://github.com/NuruddinKhoiry/blastMining
This file is a part of blastMining. blastMining is a free software: you can redistribute it and/or modify
it under the terms of GNU General Public License v3.0. blastMining is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""

import os
import sys
import numpy as np
import pandas as pd
import time
import argparse
from fastnumbers import fast_forceint
from blastMining.script import summary_df
from blastMining.script import vote_script
from blastMining.script import summary2krona
from blastMining.script import splitdf2list
from blastMining.script import create_dir
from blastMining.script import read_multidfs
from blastMining.script import get_taxid

list_tax ="99,97,95,90,85,80,75"
def add_arguments(parser):
    
    parser.add_argument('-v','--version', action='version', version='blastMining v.1.2.0')

    parser.add_argument("-i", "--input", type=str, required=True, 
        help='''blast.out file. Please use this blast outfmt 6 ONLY:
("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
[required]''')

    parser.add_argument("-o", "--outdir", dest="outdir", type=str, required=True,
        help="Output directory\n[required]")
    
    parser.add_argument("-e","--evalue",dest="evalue", action="store", default=1e-3, type=float,
        help='''Threshold of evalue
(Ignore hits if their evalues are above this threshold)
[default=1-e3]''')
    
    parser.add_argument("-txl","--taxa_level", dest="taxa_level", default=list_tax, type=str,
        help='''P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species
A comma separated list of integers as an argument
[default=99,97,95,90,85,80,75]''')
    
    parser.add_argument("-n","--topN", dest="topN",action="store",default=10,type=int,
        help="Top N hits used for voting\n[default=10]")
    
    parser.add_argument("-sm","--sample_name",dest="sample_name",action="store",default='sample',type=str,
        help='Sample name in the print out table\n[default="sample"]')

    parser.add_argument("-j", "--jobs", dest="jobs", type=int, default=1,
        help='Number of jobs to run parallelly\n[default=1]')
        
    parser.add_argument("-p", "--prefix", dest="prefix", type=str, default='vote_method',
        help="Output prefix\n[default='vote_method']")
    
    parser.add_argument("-kp","--krona_plot",dest="krona_plot", default=argparse.SUPPRESS,action='store_true',
        help='Draw krona plot\n[default=False]')

    parser.add_argument("-rm", "--rm_tmpdir", dest="rm_tmpdir", default=argparse.SUPPRESS, action='store_true',
        help='Remove temporary directory (TMPDIR)\n[default=False]')
    
    return parser

def main(args):

    start = time.time()
    
    Green='\033[0;32m'
    CEND = '\033[0m'
    
    create_dir.create_dir(args.outdir)
    
    blast = pd.read_csv(args.input, sep='\t', header=None, dtype=str)
    blast.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid"]
    
    run_taxonkit = "cat "+str(args.input)+" | cut -f 9 | taxonkit lineage -j "+str(args.jobs)+" | taxonkit reformat -P -j "+str(args.jobs)+" | csvtk -H -t cut -f 1,3 -j "+str(args.jobs)+" > "+str('./'+os.path.join(args.outdir, "TMPDIR")+"/TAXONKIT.out")
    os.system(run_taxonkit)
    
    tax = pd.read_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+"/TAXONKIT.out"), sep='\t', header=None, dtype=str)
    tax.columns = ['staxid', 'lineage']
    tax[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = tax['lineage'].str.split(';', expand=True)
    tax.drop('lineage', axis=1, inplace=True)
    
    df_merge = pd.merge(blast, tax, left_index=True, right_index=True)
    
    sdf = splitdf2list.splitdf2list(df_merge, 'qseqid', args.jobs)
    
    list_taxa = str(args.taxa_level)
    list_taxa = list_taxa.replace(' ','').replace('[','').replace(']','')
    
    for i in range(len(sdf)):
        pd.DataFrame(sdf[i]).to_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+'/_'+str(i)+'.BLASTDF'), index=None, header=True, sep = '\t')            
        with open(str('./')+os.path.join(args.outdir, "TMPDIR")+'/'+str('run_vote')+'_'+str(i)+'.sh', 'w') as output:
            y = "run_vote.py -i "+str('./'+os.path.join(args.outdir, "TMPDIR")+'/_'+str(i)+'.BLASTDF')+" -txl "+str(list_taxa)+" -e "+str(args.evalue)+" -n "+str(args.topN)+" -o "+str('./')+os.path.join(args.outdir, "TMPDIR")+"/_"+str(i)+".VOTE"
            output.write('#!/bin/bash'+'\n'+y)
    
    if int(args.jobs) == 1:
        run_vote = "bash ./"+os.path.join(args.outdir, "TMPDIR")+"/run_vote_0.sh"
        os.system(run_vote)
    else:
        run_vote = "parallel --will-cite -j "+str(round(args.jobs))+" bash ::: ./"+os.path.join(args.outdir, "TMPDIR")+"/run_vote*.sh"
        os.system(run_vote)
    
    multi_dfs = read_multidfs.read_multidfs(str('./'+os.path.join(args.outdir, "TMPDIR")), '*.VOTE')
    get_taxid.get_taxid(multi_dfs, str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.GETTAXID"))

    run_get_taxid = "cat " + str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.GETTAXID") + " | cut -f 8 | taxonkit name2taxid -j "+str(args.jobs)+" | awk '!a[$1]++' > "+str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.TAXID")
    os.system(run_get_taxid)

    df2 = pd.read_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.TAXID"), sep='\t', header= None, dtype=str)
    df2.columns = ['Species', 'staxid']
    df2['Species']=df2['Species'].astype(str)
    
    multi_dfs['Species']=multi_dfs['Species'].astype(str)
    FD_df = multi_dfs.merge(df2, on="Species")

    FD_df = FD_df[['qseqid', 'staxid', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]

    df4 = read_multidfs.read_multidfs(str('./'+os.path.join(args.outdir, "TMPDIR")), '*.VOTE')
    df5 = FD_df[["qseqid","staxid"]].merge(df4, on="qseqid", how='left')

    df5.to_csv(str('./'+args.outdir+'/'+args.prefix+'.tsv'), header=True, index=None, sep='\t')
    
    SD = summary_df.summary_df(df5, args.sample_name)
    SD.to_csv(str('./'+args.outdir+'/'+args.prefix+'.summary'), header=True, index=None, sep='\t')
    
    if hasattr(args, 'krona_plot'):
        print(Green+'\n')
        summary2krona.summary2krona(str('./'+args.outdir+'/'+args.prefix+'.summary'), str('./'+args.outdir+'/'+args.prefix+'.krona'))
        os.system("ktImportText "+str('./'+args.outdir+'/'+args.prefix+'.krona')+' -o '+str('./'+args.outdir+'/'+args.prefix+'.html'))
        print(''+CEND)
        
    remove_sh = 'rm '+os.path.join(args.outdir, "TMPDIR")+"/*.sh"
    os.system(remove_sh)
    
    if hasattr(args, 'rm_tmpdir'):
        rm_tmpdir = 'rm -r '+str('./'+os.path.join(args.outdir, "TMPDIR"))
        os.system(rm_tmpdir)

    end = time.time()
    
    print(Green+'\nFinish in '+str(round(end-start,3))+' s\n'+CEND)

if __name__ == '__main__':
    main()