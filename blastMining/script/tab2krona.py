#!/usr/bin/env python3
'''
convert TABLE.summary to TABLE.krona

***
This script is a part of blastMining program
***

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)
'''

import pandas as pd
import argparse
from argparse import RawTextHelpFormatter
import sys

def tab2krona(table_in, table_out):
    df = pd.read_csv(table_in, delimiter='\t')

    df.drop('staxid',axis=1, inplace=True)
    
    col = df.columns.to_list()[::-1]
    df2 = df[col]
    df3 = df2.Taxa.apply(lambda x: pd.Series(str(x).split(";")))
    
    df4 = pd.concat([df2, df3], axis=1)
    col += ['Kingdom','Phylum','Class','Order','Family','Genus','Species']
    df4.columns = col
    
    df4.drop('Taxa',  axis=1, inplace=True)
    
    df4.to_csv(str(table_out+'.krona'), header=True, index=None, sep='\t')

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-v", "--version",
                        help="print version and exit",action="version",version='summary2krona v.1.2.0')
    parser.add_argument("-i", "--input", required=True, dest='input',
                        help="input table")
    parser.add_argument("-o", "--output", default='OUTPUT', type=str, dest='output',
                        help="output name [default = 'OUTPUT']")
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
        
    tab2krona(args.input, args.output)
    
if __name__ == '__main__':
    main()