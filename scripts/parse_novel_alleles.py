#!/usr/bin/env python2
import logging
logger = logging.getLogger()
from collections import defaultdict
import argparse
import sys
import os
from Bio import SeqIO

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Given a list of FASTA files with novel alleles found with MentaLiST, output a FASTA with a unique list of novel alleles.")
    # parser.add_argument("-s", nargs="+", help="New scheme fasta files, to compare if the novel allees are present.")
    parser.add_argument("-f", nargs="+", help="Fasta files with novel alleles.")
    parser.add_argument("-o", type=str, help="Output Fasta file with alleles above the threshold requirement(s).")
    parser.add_argument("-t", "--threshold", type=int, default=5, help="Minimum number of different samples to appear, to include a novel allele in the output fasta.")
    parser.add_argument("-m", "--mutation", type=int, default=0, help="Also include if novel allel has equal or less than this number of mutations, regardless of times seen. Disabled by default.")
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    logger.info("Reading the new alleles  ...")
    # novel = defaultdict(lambda : defaultdict(int))
    novel = defaultdict(lambda : defaultdict(list))
    loci = set()
    for f in param.f:
        # get mutations:
        with open(f[:-2] + "txt") as mutfile:
            mutations = {locus:nmut for locus, ab, nmut, desc in [l.strip().split("\t") for l in mutfile]}
        logger.debug("Opening file %s ..." % f)
        for seq_record in SeqIO.parse(f, "fasta"):
            locus = seq_record.id
            dna = str(seq_record.seq)
            loci.add(locus)
            novel[locus][dna].append(int(mutations[locus]))

    logger.info("Writing output ...")
    output_fasta = []
    print("Locus\tAlleles found\tSamples x (mutations)")
    for locus in sorted(loci):
        output_fasta.extend([(locus,seq) for seq, l in novel[locus].items() if (len(l) >= param.threshold or min(l) <= param.mutation)])
        print("%s\t%d\t%s" % (locus, len(novel[locus]), ", ".join(["%dx (%d)" % (len(l),min(l)) for seq, l in sorted(novel[locus].items(), key=lambda x:len(x[1]), reverse=True)])))
    if param.o:
        with open(param.o, "wb") as f:
            for locus, seq in output_fasta:
                print >> f, (">%s\n%s" % (locus, seq))

    logger.info("Done.")
