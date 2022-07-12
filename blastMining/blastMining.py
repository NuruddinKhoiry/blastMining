#!/usr/bin/env python3
"""
blastMining v.1.0.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

BLAST outfmt 6 only:
("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")  

"""

__version__ = 'v.1.0.0'

import argparse
import sys
import os
from argparse import RawTextHelpFormatter

from blastMining.utility import utility
import importlib

COMMANDS = ('vote','voteSpecies','lca','besthit','full_pipeline')

def main():
    parser = argparse.ArgumentParser(description=__doc__, prog='blastMining', formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--debug', dest='debug', action='store_true', default=False, help='debug mode output [%(default)s]')
    
    subparsers = parser.add_subparsers()
    for command_name in COMMANDS:
        module = importlib.import_module('blastMining.' + command_name)
        subparser = subparsers.add_parser(command_name, 
                                        help=module.__doc__.split('\n')[1], 
                                        description=module.__doc__, 
                                        formatter_class=RawTextHelpFormatter)
        subparser.set_defaults(func=module.main)
        module.add_arguments(subparser)        
      
    '''
    ArgumentParser.parse_args(args=None, namespace=None)
    args - List of strings to parse. The default is taken from sys.argv.
    '''
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    if len(sys.argv) == 2 and str(sys.argv[1]) == 'vote':
        os.system('blastMining vote -h')
        sys.exit()
    if len(sys.argv) == 2 and str(sys.argv[1]) == 'voteSpecies':
        os.system('blastMining voteSpecies -h')
        sys.exit() 
    if len(sys.argv) == 2 and str(sys.argv[1]) == 'lca':
        os.system('blastMining lca -h')
        sys.exit()
    if len(sys.argv) == 2 and str(sys.argv[1]) == 'besthit':
        os.system('blastMining besthit -h')
        sys.exit()    
    if len(sys.argv) == 2 and str(sys.argv[1]) == 'full_pipeline':
        os.system('blastMining full_pipeline -h')
        sys.exit()

    args = parser.parse_args()
    args.logger = utility.get_logger(debug=args.debug)
    
    if not hasattr(args, 'func'):
        parser.error('Please provide the name of a subcommand to run')
    else:        
        args.func(args)

if __name__ == '__main__':
    main()
