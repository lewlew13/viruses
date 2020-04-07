#!/usr/bin/env python

__author__ = "Ben Dickins"
__version__ = "2.0.0"

import re, sys

if __name__ == "__main__":

    # ask the user for input
    homo_or_dinuc = str(input("Examine homopolymer (h) or dinucleotide (d) repeats? ")) or "h"
    min_rep_num = str(input("Minimum number of repeats [default=3]: ") or 3)
    output_name = str(input("Please write name of output file (no spaces and .txt will be appended): ")) or "output"

    # based on input choose motif list or
    if homo_or_dinuc == "h":
        motifs = ["A","C","G","T"]
    elif homo_or_dinuc == "d":
        motifs = ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"]
    else:
        print("Please answer the 1st question with 'h' or 'd'.")
        sys.exit()

    # read the FASTA file into a variable (seq)
    with open(sys.argv[1], 'r') as fasta:
        all_lines = fasta.readlines()
        chrom, seq = 'null', ''
        for line in all_lines:
            line = line.rstrip()
            if line.startswith('>'):
                chrom = line[1:]
            else:
                seq += line

    # analyse for repeats
    with open(output_name + '.txt', 'w') as outfile:
        # header lines (X2)
        # print("#settings: homo/di: {}, min_rep: {}\ngenome\tstart\tend\tmotif\tstrand".format(homo_or_dinuc, min_rep_num), file=outfile)
        # loop through each motif
        for unit in motifs:
            regex = re.compile("({}){{{},}}".format(unit, min_rep_num), re.IGNORECASE)
            discovery = [(m.start()+1, m.end(), m.group()) for m in re.finditer(regex, seq)]
            for tup in discovery:
                print(chrom, "\t", sep = "", end = "", file = outfile)
                print("\t".join([str(t) for t in tup]), end = "", file = outfile)
                print("\t+", file = outfile)

# dead darling:
# pat = re.compile(r'(A){3,}|(C){3,}|(G){3,}|(T){3,}', re.IGNORECASE)
