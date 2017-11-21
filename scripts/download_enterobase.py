#!/usr/bin/env python
import logging
logger = logging.getLogger()

import argparse
import sys
import os
import urllib
import gzip
import shutil
from multiprocessing import Pool

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Downloads and unzips all FASTA files from a specific enterobase MLST specific. ")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output folder for FASTA files.")
    parser.add_argument("-s", "--species", type=str, required=True, choices=["S","Y","E"], help="Choose the target scheme: (S)almonella, (Y)ersinia, or (E)scherichia/Shigella.")
    parser.add_argument("-y", "--type", type=str, required=True, choices=["cg","wg"], help="Choose the MLST scheme type, cgMLST (cg) or wgMLST (wg).")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads, to download in parallel. Use responsibly.")
    # parser.add_argument("files", nargs="+", help="Fasta files")
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    basefolder = os.path.dirname(sys.argv[0])
    # Read the file for the corresponding scheme
    sp = {"E":"ESC", "S":"SAL", "Y":"YER"}
    verbose = {"S":"Salmonella", "Y":"Yersinia", "E":"Escherichia/Shigella"}
    sp_code = "%swgMLST" % sp[param.species]
    tp_code = "%sMLSTv1" % param.type
    filename = os.path.join(basefolder, "%s.txt" % sp_code)

    def download_and_gunzip_locus(locus):
        fasta_file =  os.path.join(param.output, "%s.fa" % locus)
        gzip_file =  os.path.join(param.output, "%s.fa.gz" % locus)
        if os.path.isfile(fasta_file):
            logger.info("File for locus %s already exists, skipping ..." % locus)
            return
        logger.info("Downloading locus %s" % locus)
        url = "http://enterobase.warwick.ac.uk/download_data?species=%s&scheme=%s&allele=%s" % (sp_code, tp_code, locus)
        if not os.path.isfile(gzip_file):
            urllib.urlretrieve(url, gzip_file)
        with gzip.open(gzip_file, 'rb') as f_in, open(fasta_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(gzip_file)

    all_loci = []
    with open(filename) as f:
        all_loci = [obj[1].split(" ")[0] for obj in [l.strip().split("\t") for l in f if l] if obj[0].startswith(param.type)]

    logging.info("Downloading %d FASTA files for %s %sMLST scheme, it might take a while ..." % (len(all_loci), verbose[param.species], param.type))
    if not os.path.isdir(param.output):
        os.makedirs(param.output)

    # download all:
    if param.threads > 1:
      # process in parallel:
        p = Pool(param.threads)
        p.map(download_and_gunzip_locus, all_loci)
    else:
      # serial, useful for debugging
        for locus in all_loci:
            download_and_gunzip_locus(locus)



    # locus = l.strip()

            # sys.exit()
