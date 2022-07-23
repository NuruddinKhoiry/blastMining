#!/usr/bin/env python3
"""
Copyright 2022 Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

https://github.com/NuruddinKhoiry/blastMining
This file is a part of blastMining. blastMining is a free software: you can redistribute it and/or modify
it under the terms of GNU General Public License v3.0. blastMining is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""

import sys
import pandas as pd
import argparse
from fastnumbers import fast_forceint
from blastMining.script import vote_script

def main():
    
    parser = argparse.ArgumentParser(description='run_vote', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-i", "--input", type=str, required=True, 
        help="Input")

    parser.add_argument("-o", "--output", dest="output", type=str, required=True,
        help="Output")
    
    parser.add_argument("-e","--evalue", dest="evalue", action="store", default=1e-3, type=float,
        help="E-value")
    
    parser.add_argument("-txl","--taxa_level", dest="taxa_level", default=[99,97,95,90,85,80,75], type=list,
        help="Taxa level")
    
    parser.add_argument("-n","--topN", dest="topN",action="store",default=10,type=int,
        help="Top N hits")
        
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    df_merge = pd.read_csv(args.input, sep='\t', header=0)
    
    tax_lvl = args.taxa_level
    tax_level = [fast_forceint(x) for x in tax_lvl]

    DF = vote_script.vote(blast=df_merge, evalue=args.evalue, tax_level=tax_level, topN=args.topN)
    
    DF.to_csv(args.output, header=True, index=None, sep='\t')
    
if __name__ == '__main__':
    main()