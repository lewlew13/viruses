#!/usr/bin/env python

from sys import argv
from Bio import SeqIO

# stores all the CDS entries
all_entries = []

with open(argv[1], 'r') as GBFile:
    GBcds = SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)
    for cds in GBcds:
        if cds.seq is not None:
            cds.id = cds.name
            cds.description = ''
            all_entries.append(cds)

# write file
SeqIO.write(all_entries, "{}.fasta".format(argv[1].split(".g")[0]), "fasta")

# acknowledgement: https://bioinformatics.stackexchange.com/a/4366
