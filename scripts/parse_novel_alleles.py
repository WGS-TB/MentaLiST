#!/usr/bin/env python
import logging
logger = logging.getLogger()
from collections import defaultdict, Counter
import argparse
import sys
import os
from Bio import SeqIO

def mut_list_to_str(l, locus):
    c = Counter([x[1] for x in l])
    if len(c) == 1:
        return str(list(c)[0])
    return ", ".join(["%dx (%d)" % (times, mut) for (mut, times) in c.items()])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Given a list of FASTA files with novel alleles found with MentaLiST, output a FASTA with a unique list of novel alleles.")
    # parser.add_argument("-s", nargs="+", help="New scheme fasta files, to compare if the novel allees are present.")
    parser.add_argument("-f", nargs="+", help="Fasta files with novel alleles.")
    parser.add_argument("-o", type=str, help="Output Fasta file with alleles above the threshold requirement(s).")
    parser.add_argument("-t", "--threshold", type=int, default=1, help="Minimum number of different samples to appear, to include a novel allele in the output fasta.")
    parser.add_argument("-m", "--mutation", type=int, default=0, help="Also include if novel allel has equal or less than this number of mutations, regardless of times seen. Disabled by default.")
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    logger.info("Reading the new alleles  ...")
    novel = defaultdict(lambda : defaultdict(list))
    loci = set()
    for f in param.f:
        # get mutations:
        mutations = defaultdict(lambda : defaultdict(int))
        novel_ids = defaultdict(lambda : defaultdict(str))
        with open(f[:-2] + "txt") as mutfile:
            mutfile.readline()
            for l in mutfile:
                sample, locus, novel_id, ab, nmut, desc = l.strip().split("\t")
                mutations[sample][locus] = int(nmut)
                novel_ids[sample][locus] = novel_id
        logger.debug("Opening file %s ..." % f)
        
        for seq_record in SeqIO.parse(f, "fasta"):
            locus = "_".join(seq_record.id.split("_")[:-1])
            novel_id = seq_record.id.split("_")[-1]
            # if locus == "Rv0551c":
            #     import pdb; pdb.set_trace()
            dna = str(seq_record.seq)
            loci.add(locus)
            for sample in mutations.keys():
                # if mutations[sample][locus] > 0:
                if novel_ids[sample][locus] == novel_id:
                    novel[locus][dna].append((sample,mutations[sample][locus]))

    logger.info("Writing output ...")
    output_fasta = []
    output_report = []
    with open(param.o + ".txt", "w") as f:
        f.write("Locus\tAlleles found\tSamples x (mutations)\n")
        for locus in sorted(loci):
            output_fasta.extend([(locus,seq) for seq, lst in novel[locus].items() if (len(lst) >= param.threshold or min([x[1] for x in lst]) <= param.mutation)])
            f.write("%s\t%d\t%s\n" % (locus, len(novel[locus]), ", ".join(["%dx (%s)" % (len(lst), mut_list_to_str(lst, locus)) for seq, lst in sorted(novel[locus].items(), key=lambda x:len(x[1]), reverse=True)])))
            for seq, lst in sorted(novel[locus].items(), key=lambda x:len(x[1]), reverse=True):
                output_report.append((locus, len(lst), ",".join(["%s" % x[0] for x in lst])))
    with open(param.o + ".fa", "w") as f:
        for locus, seq in output_fasta:
            f.write(">%s\n%s\n" % (locus, seq))
    with open(param.o + ".samples.txt", "w") as f:
        f.write("Locus\tCount\tSamples\n")
        for tuple in output_report:
            f.write("%s\t%dx\t%s\n" % tuple)

    logger.info("Done.")
