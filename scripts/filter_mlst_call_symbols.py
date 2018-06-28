#!/usr/bin/env python
import logging
logger = logging.getLogger()

import argparse
import sys
import os
import string
import re

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filters uncertainty symbols ?,* and + from MLST call files.")
    parser.add_argument("-o", "--output", type=str, help="Output folder with filtered files.")
    parser.add_argument("files", nargs="+", help="MLST call files to filter.")
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    # make out folder if needed:
    if not os.path.isdir(param.output):
        os.makedirs(param.output)
    # remove all but digits:
    all_chars = string.maketrans('','')
    nodigs = all_chars.translate(all_chars, string.digits)
    for file in param.files:
        with open(file) as f_in, open(os.path.join(param.output,os.path.basename(file)), "wb") as f_out:
            header = f_in.readline()
            f_out.write(header)
            genotype = f_in.readline().strip().split()
            genotype = "\t".join([genotype[0]] + [x.translate(all_chars, nodigs) for x in genotype[1:]])
            # genotype = "\t".join([genotype[0]] + [x for x in genotype[1:]])
            while re.search(pattern='\t\t', string = genotype):
                genotype = re.sub(pattern='\t\t', repl='\t-1\t', string=genotype)
            f_out.write("%s\n" % genotype)
