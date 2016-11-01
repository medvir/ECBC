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
	run this function with run_child(‘string that you want to run’)
	'''
   try:
       output = subprocess.check_output(cmd, universal_newlines=True,
       shell=True,
#        executable=exe,
       stderr=subprocess.STDOUT)
   except subprocess.CalledProcessError as ee:
       output = None
   return output

def hamming_dist(str1, str2):
    '''
	returns hamming distance of two strings of identical length
	'''
    assert len(str1) == len(str2), 'hamming_dist: length differs'
    return sum([1 for a, b in zip(str1, str2) if a != b])


#############
# Arguments #
#############

readfile1, readfile2, indexfile1 = sys.argv[1:4]


########
# Code #
########

### open and unzip if necessary
if sys.argv[1].endswith('.gz'):
	itR1 = FastqGeneralIterator(gzip.open(readfile1))
else:
	itR1 = FastqGeneralIterator(open(readfile1))

if sys.argv[2].endswith('.gz'):
    itR2 = FastqGeneralIterator(gzip.open(readfile2))
else:
    itR2 = FastqGeneralIterator(open(readfile2))

if sys.argv[3].endswith('.gz'):
    itI1 = FastqGeneralIterator(gzip.open(indexfile1))
else:
	itI1 = FastqGeneralIterator(open(indexfile1))


		
		

###
fn = os.path.split(readfile1)[1]
sample_name = re.search('(.*)_L001_.*', fn).group(1)

igs = ['IgG1', 'IgG2', 'IgG3', 'IgG4', 'IgM', 'IgA', 'IgA1', 'IgA2', 'k', 'l', 'undetI1', 'undet24', 'undetklMA']
file_type = ['R1', 'R2', 'I1']

### stores all handles that will be used to write to files
handle_dict = {}
for i in igs:
    for ft in file_type:
        handle_dict['%s_%s' % (i, ft)] = open('%s_%s_%s.fastq' % (sample_name, i, ft), 'w+')

for it_obj in izip(itR1, itR2, itI1):
    assert it_obj[0][0].split()[0] == it_obj[1][0].split()[0] \
        and it_obj[1][0].split()[0] == it_obj[2][0].split()[0], \
        'Read ids do not match'

    R1s = it_obj[0][1].upper()
    R2s = it_obj[1][1].upper()
    I1s = it_obj[2][1].upper()

    # IgG1 index is 'AaGAGCACCTCt'
    if I1s.startswith('AAGAGCACCTCT'):
        ig = 'IgG1'

    # IgG3 index is 'AgGAGCACCTCt'
    elif I1s.startswith('AGGAGCACCTCT'):
        ig = 'IgG3'

    # IgG24 index is 'AgGAGCACCTCc'
    elif I1s.startswith('AGGAGCACCTCC'):

        # IgG2 constant starts with 'GCcTCC',
        # therefore 4th nucleotide reverse complement on R2 is 'G'
        if R2s[3] == 'G':
            ig = 'IgG2'

        # IgG4 constant starts with 'GCtTCC',
        # therefore 4th nucleotide reverse complement on R2 is 'A'
        elif R2s[3] == 'A':
            ig = 'IgG4'

        # remaining reads are undetermined
        else:
            ig = 'undet24'

    # klMA index is 'TTCTCGTGGAGA'
    elif hamming_dist(I1s[0:12], 'TTCTCGTGGAGA') <= 2:

        # IgM primer is 'NNNNGGTTGGGGCGGATGCACTCC',
        if hamming_dist(R2s[4:24], 'GGTTGGGGCGGATGCACTCC') <= 2:
            ig = 'IgM'

        # IgA1 primer (MM) is 'NNNNgGCGAtGACCACGTTCCCATC',
        elif hamming_dist(R2s[4:11], 'GGCGATG') <= 0:
            ig = 'IgA1'

        # IgA2 primer (MM) is 'NNNNtGCGAcGACCACGTTCCCATC',
        elif hamming_dist(R2s[4:11], 'TGCGACG') <= 0:
            ig = 'IgA2'

        # IgA only (MS) primer is 'NNNNTTGGGGCTGGTCGGGGATGC',
        elif hamming_dist(R2s[4:24], 'TTGGGGCTGGTCGGGGATGC') <= 2:
            ig = 'IgA'

        # k primer is 'NNNNAGATGGTGCAGCCACAGTTC',
        elif hamming_dist(R2s[4:24], 'AGATGGTGCAGCCACAGTTC') <= 2:
            ig = 'k'

        # l primer is 'NNNNTAGGACGGTSASCTTGGTCC' or 'NNNNGAGGACGGTCAGCTGGGTGC',
        elif hamming_dist(R2s[4:24], 'TAGGACGGTSASCTTGGTCC') <= 4 or hamming_dist(R2s[4:24], 'GAGGACGGTCAGCTGGGTGC') <= 2:
            ig = 'l'

        # remaining reads are undetermined
        else:
            ig = 'undetklMA'

     # remaining reads are undetermined
    else:
        ig = 'undetI1'

    for read_here, ft in zip(it_obj, file_type):
        handle_here = handle_dict['%s_%s' % (ig, ft)]
        top = list(read_here)
        top[0] = '@%s' % top[0]
        top.insert(2, '+')
        top.append('')
        handle_here.write('\n'.join(top))


### close files
for i in igs:
    for ft in file_type:
        handle_dict['%s_%s' % (i, ft)].close()
		
		
		
		
		
		
run_child('pandaseq -o %s -f %s -r %s -w %s_panda.fasta' % (minimal, R1, R2, Name))

