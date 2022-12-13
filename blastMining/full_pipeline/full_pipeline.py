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
import argparse
import time
import numpy as np
import pandas as pd
from argparse import RawTextHelpFormatter
from fastnumbers import fast_forceint
from blastMining.script import summary_df
from blastMining.script import vote_script
from blastMining.script import voteSpecies_script
from blastMining.script import lca_script
from blastMining.script import besthit_script
from blastMining.script import summary2krona
from blastMining.script import splitdf2list
from blastMining.script import create_dir
from blastMining.script import read_multidfs
from blastMining.script import get_taxid

list_tax ="99,97,95,90,85,80,75"
def add_arguments(parser):

    parser.add_argument('-v','--version', action='version', version='blastMining v.1.2.0')

    parser.add_argument("-i", "--input", dest="input", type=str, required=True, 
        help='input FASTA\n[required]')

    parser.add_argument("-o", "--outdir", dest="outdir", type=str, required=True,
        help="Output directory\n[required]")
		
    parser.add_argument("-bp", "--blast_param", type=str, required=True,
        default='-db nt -num_threads 1 -max_target_seqs 10', 
        help='''BLAST parameters:
Note: "-outfmt" has been defined by the package, you don't need to add it
[default="-db nt -num_threads 1 -max_target_seqs 10"]''')
				
    parser.add_argument("-m","--mining",dest="mining",action="store",default='vote',type=str,
        help='''blastMining method 
Available methods={'vote','voteSpecies','lca','besthit'}
[default='vote']''')
				
    parser.add_argument("-e","--evalue",dest="evalue",action="store",default=1e-3,type=float,
        help='''Threshold of evalue 
(Ignore hits if their evalues are above this threshold)
[default=1-e3]''')
			   
    parser.add_argument("-pi","--pident", dest="pident",action="store",default=97,type=int,
        help='''Threshold of p. identity 
(Ignore hits if their p. identities are below this threshold)
[default=97]
**Required** for "voteSpecies, lca, and besthit methods"
**Not compatible** with "vote method"''')
				
    parser.add_argument("-txl","--taxa_level", dest="taxa_level",action="store",default=list_tax,type=str,
        help='''P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species
[default=99,97,95,90,85,80,75]
**Required** for "vote method"
**Not compatible** with "voteSpecies, lca, and besthit methods"''')		   
				
    parser.add_argument("-n","--topN", dest="topN",action="store",default=10,type=int,
        help="Top N hits used for voting\n[default=10]")
    
    parser.add_argument("-sm","--sample_name",dest="sample_name",action="store",default='sample',type=str,
        help='Sample name in the print out table\n[default="sample"]')

    parser.add_argument("-j", "--jobs", dest="jobs", type=int, default=1,
        help='Number of jobs to run parallelly\n[default=1]')
        
    parser.add_argument("-p", "--prefix", dest="prefix", type=str, default='blastMining',
        help="Output prefix\n[default='blastMining']")

    parser.add_argument("-kp", "--krona_plot", dest="krona_plot", default=argparse.SUPPRESS, action='store_true',
        help='Draw krona plot\n[default=False]')
        
    parser.add_argument("-rm", "--rm_tmpdir", dest="rm_tmpdir", default=argparse.SUPPRESS, action='store_true',
        help='Remove temporary directory (TMPDIR)\n[default=False]')
		
    return parser

def main(args):

    start = time.time()
    
    CRED = '\033[91m'
    Green='\033[0;32m'
    CEND = '\033[0m'
    
    av_methods = ('vote','voteSpecies','lca','besthit')
    if args.mining not in av_methods:
        print(CRED+'\n'+str(args.mining+" method is NOT VALID")+CEND)
        print(CRED+'\n'+"Available methods={'vote','voteSpecies','lca','besthit'}\n"+CEND)
        sys.exit()
    
    create_dir.create_dir(args.outdir)
        
    print(Green+"\nRunning BLAST ...\n"+CEND)
    print("blastn -query "+str(args.input)+" "+str(args.blast_param)+" -out "+'./'+os.path.join(args.outdir, "TMPDIR")+"/"+str(args.prefix)+"_BLASTN.out"+' -outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid"')
    os.system("blastn -query "+str(args.input)+" "+str(args.blast_param)+" -out "+'./'+os.path.join(args.outdir, "TMPDIR")+"/"+str(args.prefix)+"_BLASTN.out"+' -outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid"')
    if not os.path.exists(str('./'+os.path.join(args.outdir, "TMPDIR")+"/"+str(args.prefix)+"_BLASTN.out")):
        sys.exit()
    print(Green+"\nRunning BLAST: **FINISH**\n"+CEND)
    
    print(Green+"\nRunning blastMining ...\n"+CEND)
    blast = pd.read_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+"/"+str(args.prefix)+"_BLASTN.out"), sep='\t', header=None, dtype=str)
    blast.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid"]
    
    if args.mining == 'vote':
        print(Green+"\nYou're using **vote method**\n"+CEND)

        run_taxonkit = "cat "+str('./'+os.path.join(args.outdir, "TMPDIR")+"/"+str(args.prefix)+"_BLASTN.out")+" | cut -f 9 | taxonkit lineage -j "+str(args.jobs)+" | taxonkit reformat -P -j "+str(args.jobs)+" | csvtk -H -t cut -f 1,3 -j "+str(args.jobs)+" > "+str('./'+os.path.join(args.outdir, "TMPDIR")+"/TAXONKIT.out")
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
        
    elif args.mining == 'voteSpecies':
    
        print(Green+"\nYou're using **voteSpecies method**\n"+CEND)

        run_taxonkit = "cat "+str('./'+os.path.join(args.outdir, "TMPDIR")+"/"+str(args.prefix)+"_BLASTN.out")+" | cut -f 9 | taxonkit lineage -j "+str(args.jobs)+" | taxonkit reformat -P -j "+str(args.jobs)+" | csvtk -H -t cut -f 1,3 -j "+str(args.jobs)+" > "+str('./'+os.path.join(args.outdir, "TMPDIR")+"/TAXONKIT.out")
        os.system(run_taxonkit)
        
        tax = pd.read_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+"/TAXONKIT.out"), sep='\t', header=None, dtype=str)
        tax.columns = ['staxid', 'lineage']
        tax[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = tax['lineage'].str.split(';', expand=True)
        tax.drop('lineage', axis=1, inplace=True)
        
        df_merge = pd.merge(blast, tax, left_index=True, right_index=True)
        
        sdf = splitdf2list.splitdf2list(df_merge, 'qseqid', args.jobs)
        
        for i in range(len(sdf)):
            pd.DataFrame(sdf[i]).to_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+'/_'+str(i)+'.BLASTDF'), index=None, header=True, sep = '\t')            
            with open(str('./')+os.path.join(args.outdir, "TMPDIR")+'/'+str('run_voteSpecies')+'_'+str(i)+'.sh', 'w') as output:
                y = "run_voteSpecies.py "+str('./'+os.path.join(args.outdir, "TMPDIR")+'/_'+str(i)+'.BLASTDF')+' '+str(args.pident)+" "+str(args.evalue)+" "+str(args.topN)+" "+str('./')+os.path.join(args.outdir, "TMPDIR")+"/_"+str(i)+".VOTESPECIES"
                output.write('#!/bin/bash'+'\n'+y)
        
        if int(args.jobs) == 1:
            run_voteSpecies = "bash ./"+os.path.join(args.outdir, "TMPDIR")+"/run_voteSpecies_0.sh"
            os.system(run_voteSpecies)
        else:
            run_voteSpecies = "parallel --will-cite -j "+str(round(args.jobs))+" bash ::: ./"+os.path.join(args.outdir, "TMPDIR")+"/run_voteSpecies*.sh"
            os.system(run_voteSpecies)
        
        DF = read_multidfs.read_multidfs(str('./'+os.path.join(args.outdir, "TMPDIR")), '*.VOTESPECIES')
        DF.to_csv(str('./'+args.outdir+'/'+args.prefix+'.tsv'), header=True, index=None, sep='\t')
        
        SD = summary_df.summary_df(DF, args.sample_name)
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
		
    elif args.mining == 'lca':
        print(Green+"\nYou're using **lca method**\n"+CEND)
        
        sdf = splitdf2list.splitdf2list(blast, 'qseqid', args.jobs)
    
        for i in range(len(sdf)):
            pd.DataFrame(sdf[i]).to_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+'/_'+str(i)+'.BLASTDF'), index=None, header=True, sep = '\t')            
            with open(str('./')+os.path.join(args.outdir, "TMPDIR")+'/'+str('run_lca')+'_'+str(i)+'.sh', 'w') as output:
                y = "run_lca.py "+str('./'+os.path.join(args.outdir, "TMPDIR")+'/_'+str(i)+'.BLASTDF')+' '+str(args.pident)+" "+str(args.evalue)+" "+str(args.topN)+" "+str('./')+os.path.join(args.outdir, "TMPDIR")+"/_"+str(i)+".LCA"
                output.write('#!/bin/bash'+'\n'+y)
                
        if int(args.jobs) == 1:
            run_lca = "bash ./"+os.path.join(args.outdir, "TMPDIR")+"/run_lca_0.sh"
            os.system(run_lca)
        else:
            run_lca = "parallel --will-cite -j "+str(round(args.jobs))+" bash ::: ./"+os.path.join(args.outdir, "TMPDIR")+"/run_lca*.sh"
            os.system(run_lca)
        
        DF = read_multidfs.read_multidfs(str('./'+os.path.join(args.outdir, "TMPDIR")), '*.LCA')
        DF.to_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+'/'+'LCA'), header=None, index=None, sep='\t')
        
        os.system("bash blastMining_lca.sh "+str('./'+os.path.join(args.outdir, "TMPDIR")+'/'+'LCA ')+str(args.jobs)+" "+str(args.outdir))
        os.system("bash blastMining_lca2.sh "+str('./'+os.path.join(args.outdir, "TMPDIR")+'/'+'tmp_lca ')+str(args.jobs)+" "+str(args.outdir))
        
        lca1 = pd.read_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+'/'+'tmp_lca'), sep='\t', header=None, dtype=str)
        lca1.columns = ['staxid', 'lca']
        lca2 = pd.read_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+'/'+'tmp_lca2'), sep='\t', header=None, dtype=str)
        lca2.columns = ['lca', 'lineage']
        
        DF2 = pd.concat([lca1,lca2], axis = 1)
        DF2.columns = ['staxids','lca_x','lca_y','lineage']
        DF3 = pd.concat([DF.reset_index(drop=True), DF2], axis = 1)
        DF3[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = DF3['lineage'].str.split(';', expand=True)
        DF3.drop(['staxid','staxids','lca_x','lca_y','lineage'], axis=1, inplace=True)
        DF3.to_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.df3"), header=True, index=None, sep='\t')
        get_taxid.get_taxid(DF3, str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.GETTAXID"))

        run_get_taxid = "cat " + str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.GETTAXID") + " | cut -f 8 | taxonkit name2taxid -j "+str(args.jobs)+" | awk '!a[$1]++' > "+str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.TAXID")
        os.system(run_get_taxid)

        df2 = pd.read_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.TAXID"), sep='\t', header= None, dtype=str)
        df2.columns = ['Species', 'staxid']
        df2['Species']=df2['Species'].astype(str)
        
        DF3['Species']=DF3['Species'].astype(str)
        FD_df = DF3.merge(df2, on="Species")

        FD_df = FD_df[['qseqid', 'staxid', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]

        df4 = pd.read_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+"/temp.df3"), sep='\t', header=0)
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
	
    elif args.mining == 'besthit':
        print(Green+"\nYou're using **besthit method**\n"+CEND)
        
        run_taxonkit = "cat "+str('./'+os.path.join(args.outdir, "TMPDIR")+"/"+str(args.prefix)+"_BLASTN.out")+" | cut -f 9 | taxonkit lineage -j "+str(args.jobs)+" | taxonkit reformat -P -j "+str(args.jobs)+" | csvtk -H -t cut -f 1,3 -j "+str(args.jobs)+" > "+str('./'+os.path.join(args.outdir, "TMPDIR")+"/TAXONKIT.out")
        os.system(run_taxonkit)
        
        tax = pd.read_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+"/TAXONKIT.out"), sep='\t', header=None, dtype=str)
        tax.columns = ['staxid', 'lineage']
        tax[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = tax['lineage'].str.split(';', expand=True)
        tax.drop('lineage', axis=1, inplace=True)
        
        df_merge = pd.merge(blast, tax, left_index=True, right_index=True)
        
        sdf = splitdf2list.splitdf2list(df_merge, 'qseqid', args.jobs)
        
        for i in range(len(sdf)):
            pd.DataFrame(sdf[i]).to_csv(str('./'+os.path.join(args.outdir, "TMPDIR")+'/_'+str(i)+'.BLASTDF'), index=None, header=True, sep = '\t')            
            with open(str('./')+os.path.join(args.outdir, "TMPDIR")+'/'+str('run_besthit')+'_'+str(i)+'.sh', 'w') as output:
                y = "run_besthit.py "+str('./'+os.path.join(args.outdir, "TMPDIR")+'/_'+str(i)+'.BLASTDF')+' '+str(args.pident)+" "+str(args.evalue)+" "+str(args.topN)+" "+str('./')+os.path.join(args.outdir, "TMPDIR")+"/_"+str(i)+".BESTHIT"
                output.write('#!/bin/bash'+'\n'+y)
        
        if int(args.jobs) == 1:
            run_besthit = "bash ./"+os.path.join(args.outdir, "TMPDIR")+"/run_besthit_0.sh"
            os.system(run_besthit)
        else:
            run_besthit = "parallel --will-cite -j "+str(round(args.jobs))+" bash ::: ./"+os.path.join(args.outdir, "TMPDIR")+"/run_besthit*.sh"
            os.system(run_besthit)
        
        DF = read_multidfs.read_multidfs(str('./'+os.path.join(args.outdir, "TMPDIR")), '*.BESTHIT')
        DF.to_csv(str('./'+args.outdir+'/'+args.prefix+'.tsv'), header=True, index=None, sep='\t')
        
        SD = summary_df.summary_df(DF, args.sample_name)
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
    print(Green+'\nFinish in '+str(round(end-start, 3))+' s\n'+CEND)
    
if __name__ == '__main__':
    main()