#!/usr/bin/env python3

import os
import re
import sys
import logging

def get_logger(debug=False):
    # CRITICAL > ERROR > WARNING > INFO > DEBUG
    formatter = logging.Formatter(
        "\n%(asctime)s - %(name)s - %(levelname)s - \n%(message)s")

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)  # must be DEBUG, then 'ch' below works.
    # logger.setFormatter(formatter)

    fh = logging.FileHandler(os.path.basename(sys.argv[0]) + '.log')
    if debug:
        fh.setLevel(logging.DEBUG)
    else:
        fh.setLevel(logging.INFO)  # INFO level goes to the log file
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    if debug:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)  # only INFO level will output on screen
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger