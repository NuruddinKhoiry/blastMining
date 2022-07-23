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

def create_dir(OUTDIR):
    
    CRED = '\033[91m'
    CEND = '\033[0m'
    
    if os.path.exists(str(OUTDIR)):
        print(CRED+'\n'+str(OUTDIR)+' diretory is already EXIST. Please use another name\n'+CEND)
        sys.exit()
    else:
        os.makedirs(str(OUTDIR))
        os.makedirs(os.path.join(str(OUTDIR), "TMPDIR"))