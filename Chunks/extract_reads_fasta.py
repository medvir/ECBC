#!/opt/python2.7/bin/python2.7
import sys
import gzip
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator








ids_file = sys.argv[1]
ids = set([l.strip() for l in open(ids_file)])

fasta_file = sys.argv[2]
fa_handle = open(fasta_file)

out_file = sys.argv[2].replace('.fasta', '_matching.fasta')
output_handle = open(out_file, 'w+')

c = 0
for title, seq in FastqGeneralIterator(fa_handle):
	print(seq)
	break
	if title.split()[0] in ids:
		c += 1
		output_handle.write("@%s\n%s\n" % (title, seq))
		if c % 10000 == 0:
			print >> sys.stderr, 'written %d reads' % c

output_handle.close()
