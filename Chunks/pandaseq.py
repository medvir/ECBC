#!/opt/python2.7/bin/python2.7

### Version xx
'''
Assigns the subtypes of a demultiplexed sample based on the I1 and the R2 read. Writes _R1, and I1 files.
- IgG1_R1/R2/I1
- IgG2_R1/R2/I1
- IgG3_R1/R2/I1
- IgG4_R1/R2/I1
- IgM_R1/R2/I1
- IgA1_R1/R2/I1
- IgA2_R1/R2/I1
- IgA_R1/R2/I1
- k_R1/R2/I1
- l_R1/R2/I1
- undetI1_R1/R2/I1
- undet24_R1/R2/I1
- undetklMA_R1/R2/I1
'''

import sys
import gzip
import os
import re
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import subprocess


#############
# Functions #
#############

def run_child(cmd, exe='/bin/bash'):
	'''
	use subrocess.check_output to run an external program with arguments
	run this function with run_child (string that you want to run)^
	'''
	try:
		output = subprocess.check_output(cmd, universal_newlines=True,
        shell=True,
	    #executable=exe,
        stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError as ee:
		output = None
	return output
	

#############
# Arguments #
#############

readfile1 = sys.argv[1]
minimal_overlap = 10


########
# Code #
########

### files
readfile2 = readfile1.replace("_R1_001.fastq", "_R2_001.fastq")
indexfile1 = readfile1.replace("_R1_001.fastq", "_I1_001.fastq")
indexfile2 = readfile1.replace("_R1_001.fastq", "_I2_001.fastq")

### unzip if necessary
if readfile1.endswith('.gz'):
	cmd = 'gunzip -c %s > R1_fasta' % readfile1
	run_child(cmd)
else:
	R1_fasta = readfile1

if readfile2.endswith('.gz'):
	cmd = 'gunzip -c %s > R2_fasta' % readfile1
	run_child(cmd)
else:
	R2_fasta = readfile2


name = readfile1.replace("_R1_001.fastq", "")	
cmd = 'pandaseq -o %s -f %s -r %s -w %s_panda.fasta' % (minimal_overlap, R1_fasta, R2_fasta, name)
print cmd
run_child(cmd)


