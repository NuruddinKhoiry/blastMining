import os
import sys
import argparse
import numpy as np
import pandas as pd
from fastnumbers import fast_forceint
from blastMining.script import summary_df
from blastMining.script import vote_script
from blastMining.script import voteSpecies_script
from blastMining.script import lca_script
from blastMining.script import besthit_script
from argparse import RawTextHelpFormatter

def add_arguments(parser):
    
    parser.add_argument("-i", "--input", type=str, required=True, 
        help='input FASTA')
		
    parser.add_argument("-bp", "--blast_param", type=str, required=True,
        default='-db nt -num_threads 1 -max_target_seqs 10', 
        help='''BLAST parameters:\n
				Note: "-outfmt" has been defined by the package, you don't need to add it\n 
                (-outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid")\n
			    [default="-db nt -num_threads 1 -max_target_seqs 10"]''')
				
    parser.add_argument("-m","--mining",dest="mining",action="store",default='vote',type=str,
        help='''blast mining method\n 
               Available methods={'vote','voteSpecies','lca','besthit'}\n
               [default='vote']''')
				
    parser.add_argument("-e","--evalue",dest="evalue",action="store",default=1e-3,type=float,
        help='''Threshold of evalue\n 
               (Ignore hits if their evalues are above this threshold\n)
               [default=1-e3]''')
			   
    parser.add_argument("-pi","--pident", dest="pident",action="store",default=97,type=int,
        help='''Threshold of p. identity\n 
                (Ignore hits if their p. identities are below this threshold)\n
				[default=97]\n
                **Not compatible** with "vote method"''')
				
    parser.add_argument("-txl","--taxa_level", dest="taxa_level",action="store",default=[99,97,95,90,85,80,75],type=list,
        help='''P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species\n
                [default=99,97,95,90,85,80,75]\n
				**Required** for "vote method"''')		   
				
    parser.add_argument("-n","--topN", dest="topN",action="store",default=10,type=int,
        help="Top N hits used for voting [default=10]")
    
    parser.add_argument("-sm","--sample_name",dest="sample_name",action="store",default='sample',type=str,
        help='Sample name in the print out table [default="sample"]')
        
    parser.add_argument("-o", "--output", type=str, required=True,
        help="output")
		
    return parser

def main(args):
    
    os.system('echo "Running...\n"')
    os.system("echo blastn -query "+str(args.input)+" "+str(args.blast_param)+" -out "+str(args.output)+".out"+' -outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid"')
    os.system("blastn -query "+str(args.input)+" "+str(args.blast_param)+" -out "+str(args.output)+".out"+' -outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid"')
    blast = pd.read_csv(str(args.output+".out"), sep='\t', header=None, dtype=str)
    blast.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid"]
    
    if args.mining == 'vote':
        print("You're using vote method ...\n")
        os.system("cat "+str(args.output+".out")+" | cut -f 9 | taxonkit lineage | taxonkit reformat -P | csvtk -H -t cut -f 1,3 > output.tmp")
        
        tax = pd.read_csv('output.tmp', sep='\t', header=None, dtype=str)
        tax.columns = ['staxid', 'lineage']
        tax[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = tax['lineage'].str.split(';', expand=True)
        tax.drop('lineage', axis=1, inplace=True)
        DF = vote_script.vote(blast=blast, tax=tax, evalue=args.evalue, tax_level=args.taxa_level, topN=args.topN)
        DF.to_csv(str(args.output+'.tsv'), header=True, index=None, sep='\t')
        
        SD = summary_df.summary_df(DF, args.sample_name)
        SD.to_csv(str(args.output+'.summary'), header=True, index=None, sep='\t')
		
    elif args.mining == 'voteSpecies':
        print("You're using voteSpecies method ...\n")
        os.system("cat "+str(args.output+".out")+" | cut -f 9 | taxonkit lineage | taxonkit reformat -P | csvtk -H -t cut -f 1,3 > output.tmp")
        
        tax = pd.read_csv('output.tmp', sep='\t', header=None, dtype=str)
        tax.columns = ['staxid', 'lineage']
        tax[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = tax['lineage'].str.split(';', expand=True)
        tax.drop('lineage', axis=1, inplace=True)
        
        DF = voteSpecies_script.voteSpecies(blast=blast, tax=tax, pident=args.pident, evalue=args.evalue, topN=args.topN)
        DF.to_csv(str(args.output+'.tsv'), header=True, index=None, sep='\t')
        
        SD = summary_df.summary_df(DF, args.sample_name)
        SD.to_csv(str(args.output+'.summary'), header=True, index=None, sep='\t')
		
    elif args.mining == 'lca':
        print("You're using lca method ...\n")
        
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
	
    elif args.mining == 'besthit':
        print("You're using besthit method ...\n")
        
        os.system("cat "+str(args.output+".out")+" | cut -f 9 | taxonkit lineage | taxonkit reformat -P | csvtk -H -t cut -f 1,3 > output.tmp")
        
        tax = pd.read_csv('output.tmp', sep='\t', header=None, dtype=str)
        tax.columns = ['staxid', 'lineage']
        tax[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']] = tax['lineage'].str.split(';', expand=True)
        tax.drop('lineage', axis=1, inplace=True)
        
        DF = besthit_script.besthit(blast=blast, tax=tax, pident=args.pident, evalue=args.evalue, topN=args.topN)
        DF.to_csv(str(args.output+'.tsv'), header=True, index=None, sep='\t')
        
        SD = summary_df.summary_df(DF, args.sample_name)
        SD.to_csv(str(args.output+'.summary'), header=True, index=None, sep='\t')
    
    else:
        print(str(args.mining+" method is NOT VALID"))

if __name__ == '__main__':
    main()