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
from blastMining.script import voteSpecies_script

def main():
    
    df_merge = pd.read_csv(sys.argv[1], sep='\t', header=0, dtype=str)
    
    DF = voteSpecies_script.voteSpecies(blast=df_merge, pident=float(sys.argv[2]), evalue=float(sys.argv[3]), topN=int(sys.argv[4]))
    DF.to_csv(sys.argv[5], header=True, index=None, sep='\t')
    
if __name__ == '__main__':
    main()