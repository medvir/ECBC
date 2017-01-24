#! /opt/python2.7/bin/python2.7

### Version 2
'''
Assigns the subtypes of a demultiplexed sample based on I1 and R2 reads.
Extract Error Correcting BarCodes and add in front of R1.
Stitches R1 and R2 togoether with PandaSeq.
'''

### import
import sys
import gzip
import os
import re
import subprocess
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2


### Arguments and variables
ECBC_length = 12
random_nucs = 4 # at start of R1, between ECBC and constant for klMA
minimal_overlap = 10 # for pandaseq
readfile1  = sys.argv[1] # format NAME_S#_L001_R1_001.fastq[.gz]
readfile2 = readfile1.replace("L001_R1_001", "L001_R2_001")
indexfile1 = readfile1.replace("L001_R1_001", "L001_I1_001")

### Ig subtypes and read files
igs = ['IgG1', 'IgG2', 'IgG3', 'IgG4', 'undetIgG', 'IgM', 'IgA', 'k', 'l', 'undetklMA']
reads = ['R1', 'R2', 'I1']

### Sample_name
sample_name = os.path.split(readfile1)[1]
sample_name = sample_name.replace("_L001_R1_001.fastq.gz", "")
sample_name = sample_name.replace("_L001_R1_001.fastq", "")
print(sample_name)

### Read primers and convert to str
fwd_H = [str(s.seq) for s in SeqIO.parse(open('%s/%s' % (sys.path[0], 'primer_fwd_H.fasta')),'fasta')]
fwd_k = [str(s.seq) for s in SeqIO.parse(open('%s/%s' % (sys.path[0], 'primer_fwd_k.fasta')),'fasta')]
fwd_l = [str(s.seq) for s in SeqIO.parse(open('%s/%s' % (sys.path[0], 'primer_fwd_l.fasta')),'fasta')]

### Functions
def run_child(cmd, exe='/bin/bash'):
	'''
	use subrocess.check_output to run an external program with arguments
	run this function with run_child('string that you want to run')
	'''
	try:
		output = subprocess.check_output(cmd, universal_newlines=True,
		shell=True,
		#executable=exe,
		stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError as ee:
		output = None
	return output

def hamming_dist(str1, str2):
	'''
	returns hamming distance of two strings of identical length
	'''
	#assert len(str1) == len(str2), 'hamming_dist: length differs' # not checking for equal lenght at the moment!
	return sum([1 for a, b in zip(str1, str2) if a != b])

def trim_primer(seq, set):
	'''
	trim primer
	'''
	for f in set:
		l = len(f)
		#w = count_wobbles(f)
		if len(seq) >= l:
			print f
			print seq[random_nucs:l+random_nucs]
			if hamming_dist(seq[random_nucs:l+random_nucs], f) <= 5:
			#if alignment_penalty(seq[:l], f) <= m + w:
				return l + random_nucs
				break

### Read sequecing files
if readfile1.endswith('.gz'):
	itR1 = FastqGeneralIterator(gzip.open(readfile1))
else:
	itR1 = FastqGeneralIterator(open(readfile1))

if readfile2.endswith('.gz'):
	itR2 = FastqGeneralIterator(gzip.open(readfile2))
else:
	itR2 = FastqGeneralIterator(open(readfile2))

if indexfile1.endswith('.gz'):
	itI1 = FastqGeneralIterator(gzip.open(indexfile1))
else:
	itI1 = FastqGeneralIterator(open(indexfile1))

### stores all handles that will be used to write to files
handle_dict = {}
for i in igs:
	for r in reads:
		handle_dict['%s_%s' % (i, r)] = open('%s_%s_%s.fastq' % (sample_name, i, r), 'w+')

### Iterator
for it_obj in izip(itR1, itR2, itI1):
	assert it_obj[0][0].split()[0] == it_obj[1][0].split()[0] and it_obj[1][0].split()[0] == it_obj[2][0].split()[0], 'Read ids do not match'
	R1s = it_obj[0][1].upper()
	R2s = it_obj[1][1].upper()
	I1s = it_obj[2][1].upper()
	
	# # # # # # # # # # # # # # # # # # # # # #
	# it_obj[0][x] = Read 1                   #
	# it_obj[1][x] = Read 2                   #
	# it_obj[2][x] = Index 1                  #
	#                                         #
	# it_obj[x][0] = Name                     #
	# it_obj[x][1] = Sequence --> R1s/R2s/I1s #
	# it_obj[x][2] = Quality                  #
	# # # # # # # # # # # # # # # # # # # # # #
	
### Assinging Ig subtypes
	### klMA index is 'TAAGGCGAGAGC', allowing 1 mismatch
	if hamming_dist(I1s[0:12], 'TAAGGCGAGAGC') <= 1:
	    
		### Igl primer is 'ECBC_NNNN_CYAGTGTGGCCTTGTTGGCTTGR'
		if hamming_dist(R2s[ECBC_length+random_nucs:39], 'CYAGTGTGGCCTTGTTGGCTTGR') <= 4:
			ig = 'l'
			R2 = '@%s' % it_obj[0][0], it_obj[1][1][39:], '+', it_obj[1][2][39:] # R2 sequence and quality trim
			n = trim_primer(it_obj[0][1], fwd_l)
			R1 = '@%s' % it_obj[0][0], it_obj[0][1][n:], '+', it_obj[0][2][n:] # R1 sequence and quality trim
						
		### Igk primer is 'ECBC_NNNN_xxxxxxxxxxxxxxxxxxxxxxx'
		elif hamming_dist(R2s[ECBC_length+random_nucs:39], 'xxxxxxxxxxxxxxxxxxxxxxx') <= 3:
			ig = 'k'
			R2 = '@%s' % it_obj[0][0], it_obj[1][1][39:], '+', it_obj[1][2][39:] # R2 sequence and quality trim
			n = trim_primer(it_obj[0][1], fwd_k)
			R1 = '@%s' % it_obj[0][0], it_obj[0][1][n:], '+', it_obj[0][2][n:] # R1 sequence and quality trim

		### IgM primer is 'ECBC_NNNN_NNNNGGTTGGGGCGGATGCACTCC',
		elif hamming_dist(R2s[ECBC_length+random_nucs:39], 'GGTTGGGGCGGATGCACTCC') <= 2:
			ig = 'IgM'
			R2 = '@%s' % it_obj[0][0], it_obj[1][1][39:], '+', it_obj[1][2][39:] # R2 sequence and quality trim
			n = trim_primer(it_obj[0][1], fwd_H)
			R1 = '@%s' % it_obj[0][0], it_obj[0][1][n:], '+', it_obj[0][2][n:] # R1 sequence and quality trim

        ### IgA primer is 'ECBC_NNNN_GCATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCKRCAGCACCCMSCMAGATGGGAACGTGGTCRTC',
		elif hamming_dist(R2s[ECBC_length+random_nucs:88], 'GAYGACCACGTTCCCATCTKGSKGGGTGCTGYMGAGGCTCAGCGGGAAGACCTTGGGGCTGGTCGGGGATGC') <= 10:
			ig = 'IgA'
			R2 = '@%s' % it_obj[0][0], it_obj[1][1][88:], '+', it_obj[1][2][88:] # R2 sequence and quality trim
			n = trim_primer(it_obj[0][1], fwd_H)
			R1 = '@%s' % it_obj[0][0], it_obj[0][1][n:], '+', it_obj[0][2][n:] # R1 sequence and quality trim

		### undetermined for klMA
		else:
			ig = 'undetklMA'
			R2 = '@%s' % it_obj[0][0], it_obj[1][1], '+', it_obj[1][2] # R2 no trim
			R1 = '@%s' % it_obj[0][0], it_obj[0][1], '+', it_obj[0][2] # R1 no trim
	
	 	### ECBC for klMA is [:ECBC_length] of R2, add ECBC sequence in front of R1
		R1 = '@%s' % it_obj[0][0], it_obj[1][1][:ECBC_length] + R1[1], '+', it_obj[1][2][:ECBC_length] + R1[3]
		
	else: 		
		### IgG1 constant is 'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGG'
		if hamming_dist(R2s[0:61], 'CCCCAGAGGTGCTCTTGGAGGAGGGTGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAGGC') <= 2:
			ig = 'IgG1'
			R2 = '@%s' % it_obj[0][0], it_obj[1][1][61:], '+', it_obj[1][2][61:] # R2 sequence and quality trim
			
		### IgG2 constant is 'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGA'
		elif hamming_dist(R2s[0:61], 'TCTCGGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAGGC') <= 1:
			ig = 'IgG2'
			R2 = '@%s' % it_obj[0][0], it_obj[1][1][61:], '+', it_obj[1][2][61:] # R2 sequence and quality trim
			
		### IgG3 constant is 'GCTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGG' or 'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGG'
		elif hamming_dist(R2s[0:61], 'CCCCAGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAAGC') <= 1 or hamming_dist(R2s[0:61], 'CCCCAGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAGGC') <= 1:
			ig = 'IgG3'	
			R2 = '@%s' % it_obj[0][0], it_obj[1][1][61:], '+', it_obj[1][2][61:] # R2 sequence and quality trim
			
		### IgG4 constant is 'GCTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGA' or 'GCTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGA'
		elif hamming_dist(R2s[0:61], 'TCTCGGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACGGATGGGCCCTTGGTGGAAGC') <= 1 or hamming_dist(R2s[0:61], 'TCTCGGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAAGC') <= 1:
			ig = 'IgG4'
			R2 = '@%s' % it_obj[0][0], it_obj[1][1][61:], '+', it_obj[1][2][61:] # R2 sequence and quality trim
			
		### undet for IgG		
		else:
			ig = 'undetIgG'
			R2 = '@%s' % it_obj[0][0], it_obj[1][1], '+', it_obj[1][2] # R2 no trim
	
		n = trim_primer(it_obj[0][1], fwd_H)
		print n
		R1 = '@%s' % it_obj[0][0], it_obj[0][1][n:], '+', it_obj[0][2][n:] # R1 sequence and quality trim
			
		### ECBC for IgGs is equal to I1, add ECBC sequence in front of R1
		R1 = '@%s' % it_obj[0][0], it_obj[2][1][:ECBC_length] + R1[1], '+', it_obj[2][2][:ECBC_length] + R1[3]

	### writing R1
	handle_here = handle_dict['%s_R1' % ig]
	handle_here.write('\n'.join(R1) + '\n')
	
	### writing R2
	handle_here = handle_dict['%s_R2' % ig]
	handle_here.write('\n'.join(R2) + '\n')

### close files
for i in igs:
	for r in reads:
		handle_dict['%s_%s' % (i, r)].close()

### remove temp files
for i in igs:
	os.remove('%s_%s_I1.fastq' % (sample_name, i))

### count sequences in fasta and fastq files
print(run_child('%s/%s' % (sys.path[0], 'seq_count.sh')))


### run PandaSeq 
for ig in igs:
	fwd = '%s_%s_R1.fastq' % (sample_name, ig)
	rev = '%s_%s_R2.fastq' % (sample_name, ig)
	out = '%s_%s_ECBC-panda.fasta' % (sample_name, ig)
	cmd = 'pandaseq -o %s -f %s -r %s -w %s' % (minimal_overlap, fwd, rev, out)
	print(cmd)
	run_child(cmd)
	
	### delete R1 and R2
	#os.remove(fwd)
	#os.remove(rev)

### count sequences in fasta and fastq files
print(run_child('%s/%s' % (sys.path[0], 'seq_count.sh')))
