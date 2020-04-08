#!/usr/bin/env python3

__author__ = "Ben Dickins"
__email__ = "ben.dickins@cantab.net"
__status__ = "Prototype"
__version__ = "1.0.0"

from sys import argv
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# define regex to clear list format
rx = re.compile("(\'|\[|\])")

# functions
def overlap(left, right):
	pquity = min(left[1], right[1]) - max(left[0], right[0]) # no +1 required because python has -1'd from start
	if pquity > 0:
		return pquity
	else:
		return 0

def getlen(sometuple):
	return sometuple[1] - sometuple[0] # pythonic numbering

# main program loop
if __name__ == "__main__":

	with open(argv[1], 'r') as repeats:
		replocs = []
		for line in repeats:
			line = line.split(' ')
			try:
				replocs.append((int(line[0]), int(line[1])))
			except:
				pass

	outtab = open(argv[2].split('_')[0]+'.txt', 'w')
	print("GOA\tRepOver\tRepCov\tGene\tLocus\tGeneLen\tProduct", file=outtab)

	with open(argv[2], 'r') as genome:

		for seq_record in SeqIO.parse(genome, 'genbank'):
			for feature in seq_record.features:
				if feature.type == 'CDS':
					goaccn = [item for item in feature.qualifiers.get('db_xref',[]) if "GOA" in item]
					if goaccn == []:
						goaccn = 'NA'
					cdstuple = feature.location.start, feature.location.end
					hits = sum([overlap(reptuple,cdstuple) > 0 for reptuple in replocs])
					coverage = sum([overlap(reptuple,cdstuple) for reptuple in replocs])
					genedata = [goaccn, hits, coverage, feature.qualifiers.get('gene','NA'), \
					feature.qualifiers['locus_tag'], getlen(cdstuple), feature.qualifiers['product']]
					print('\t'.join(rx.sub("", str(x)) for x in genedata), file=outtab)
	outtab.close()
