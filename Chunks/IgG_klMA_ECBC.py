#! /opt/python2.7/bin/python2.7

### Version 0.1
'''
Assigns the subtypes of a demultiplexed sample based on I1 and R2 reads.
Extract Error Correcting BarCodes and add in front of R1.
Stitches R1 and R2 togoether with PandaSeq.
Runs ab-analysis (Brian Briney)
Runs clonify (Brian Briney)
'''

### import
import sys
import gzip
import os
import re
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import subprocess



# # # # # # # # # # # # #
# Arguments & Variables #
# # # # # # # # # # # # #

ECBC_length = 12
readfile1  = sys.argv[1] # format NAME_S#_L001_R1_001.fastq[.gz]

readfile2 = readfile1.replace("L001_R1_001", "L001_R2_001")
indexfile1 = readfile1.replace("L001_R1_001", "L001_I1_001")

### Ig subtypes and read files
igs = ['IgG1', 'IgG2', 'IgG3', 'IgG4', 'undetIgG', 'IgM', 'IgA', 'IgA1', 'IgA2', 'k', 'l',  'undetklMA'] # do not change order!
reads = ['R1', 'R2', 'I1']

### sample_name
sample_name = os.path.split(readfile1)[1]
sample_name = sample_name.replace("_L001_R1_001.fastq.gz", "")
sample_name = sample_name.replace("_L001_R1_001.fastq", "")
print(sample_name)



# # # # # # #
# Functions #
# # # # # # #

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



# # # # # # # # # # # # #
# Assinging Ig subtypes #
# # # # # # # # # # # # #

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

for it_obj in izip(itR1, itR2, itI1):
	assert it_obj[0][0].split()[0] == it_obj[1][0].split()[0] and it_obj[1][0].split()[0] == it_obj[2][0].split()[0], 'Read ids do not match'
	
	R1s = it_obj[0][1].upper()
	R2s = it_obj[1][1].upper()
	I1s = it_obj[2][1].upper()


### Assinging to ECBC, IgG or klMA
	### klMA index is 'TAAGGCGAGAGC', allowing 1 mismatch'
	if hamming_dist(I1s[0:12], 'TAAGGCGAGAGC') <= 1:
	    
		### Igl primer is '12*ECBC_4*N_CYAGTGTGGCCTTGTTGGCTTGR'
		if hamming_dist(R2s[16:39], 'CYAGTGTGGCCTTGTTGGCTTGR') <= 4:
			ig = 'l'
			
		### Igk primer is '12*ECBC_4*N_xxxxxxxxxxxxxxxxxxxxxxx'
		elif hamming_dist(R2s[16:39], 'xxxxxxxxxxxxxxxxxxxxxxx') <= 3:
			ig = 'k'		

		### IgM primer is 'NNNNGGTTGGGGCGGATGCACTCC',
		elif hamming_dist(R2s[16:39], 'GGTTGGGGCGGATGCACTCC') <= 2:
			ig = 'IgM'

        ### IgA1 primer (MM) is 'NNNNgGCGAtGACCACGTTCCCATC',
		elif hamming_dist(R2s[16:39], 'GGCGATG') <= 0:
			ig = 'IgA1'

        ### IgA2 primer (MM) is 'NNNNtGCGAcGACCACGTTCCCATC',
		elif hamming_dist(R2s[16:39], 'TGCGACG') <= 0:
			ig = 'IgA2'

        ### IgA only (MS) primer is 'GCATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCKRCAGCACCCMSCMAGATGGGAACGTGGTCRTC',
		elif hamming_dist(R2s[16:88], 'GAYGACCACGTTCCCATCTKGSKGGGTGCTGYMGAGGCTCAGCGGGAAGACCTTGGGGCTGGTCGGGGATGC') <= 10:
			ig = 'IgA'

		### undet for klMA
		else:
			ig = 'undetklMA'		
	
	else: 		
		### IgG1 constant is 'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGG'
		if hamming_dist(R2s[0:61], 'CCCCAGAGGTGCTCTTGGAGGAGGGTGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAGGC') <= 2:
			ig = 'IgG1'
			
		### IgG2 constant is 'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGA'
		elif hamming_dist(R2s[0:61], 'TCTCGGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAGGC') <= 1:
			ig = 'IgG2'
			
		### IgG3 constant is 'GCTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGG' or 'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGG'
		elif hamming_dist(R2s[0:61], 'CCCCAGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAAGC') <= 1 or hamming_dist(R2s[0:61], 'CCCCAGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAGGC') <= 1:
			ig = 'IgG3'	
			
		### IgG4 constant is 'GCTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGA' or 'GCTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGA'
		elif hamming_dist(R2s[0:61], 'TCTCGGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACGGATGGGCCCTTGGTGGAAGC') <= 1 or hamming_dist(R2s[0:61], 'TCTCGGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAAGC') <= 1:
			ig = 'IgG4'
		
		### undet for IgG		
		else:
			ig = 'undetIgG'
	
	### writing R1, R2, I1
	for read_here, r in zip(it_obj, reads):
		handle_here = handle_dict['%s_%s' % (ig, r)]
		top = list(read_here)
		top[0] = '@%s' % top[0]
		top.insert(2, '+')
		top.append('')
		handle_here.write('\n'.join(top))
		
### close files
for i in igs:
	for r in reads:
		handle_dict['%s_%s' % (i, r)].close()

### count sequencs in fasta and fastq files
print('\n*** count after subtype assignment ***\n')
print(run_child('seq_count.sh'))



# # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # #  # # # # #
# Adding ECBC in front of R1, removing 4N, writing new files _ECBC-R1.fastq #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

### for IgGs ECBC is equal to I1 (or [:ECBC_length] of I1), write new file _ECBC-R1.fastq
for i in igs[:5]:
	itR1 = FastqGeneralIterator(open('%s_%s_R1.fastq' % (sample_name, i)))
	itI1 = FastqGeneralIterator(open('%s_%s_I1.fastq' % (sample_name, i)))
	file = open('%s_%s_ECBC-R1.fastq' % (sample_name, i), 'w+')
	for it_obj in izip(itR1, itI1):
		assert it_obj[0][0].split()[0] == it_obj[1][0].split()[0], 'Read ids do not match'
		#top = '@%s' % it_obj[0][0], it_obj[1][1][:ECBC_length] + it_obj[0][1][4:], '+', it_obj[1][2][:ECBC_length] + it_obj[0][2][4:]
		top = '@%s' % it_obj[0][0], it_obj[0][1][4:], '+', it_obj[0][2][4:] ################### NOT ADDING ECBC !!!!!!!
		file.write('\n'.join(top) + '\n')
	file.close()

### for klMA ECBC is [:ECBC_length] of R2
for i in igs[5:]:
	itR1 = FastqGeneralIterator(open('%s_%s_R1.fastq' % (sample_name, i)))
	itR2 = FastqGeneralIterator(open('%s_%s_R2.fastq' % (sample_name, i)))
	file = open('%s_%s_ECBC-R1.fastq' % (sample_name, i), 'w+')
	for it_obj in izip(itR1, itR2):
		assert it_obj[0][0].split()[0] == it_obj[1][0].split()[0], 'Read ids do not match'
		#top = '@%s' % it_obj[0][0], it_obj[1][1][:ECBC_length] + it_obj[0][1][4:], '+', it_obj[1][2][:ECBC_length] + it_obj[0][2][4:]
		top = '@%s' % it_obj[0][0], it_obj[0][1][4:], '+', it_obj[0][2][4:] ################### NOT ADDING ECBC !!!!!!!
		file.write('\n'.join(top) + '\n')
	file.close()

### delete I1 and R1
for i in igs:
	os.remove('%s_%s_I1.fastq' % (sample_name, i))
	os.remove('%s_%s_R1.fastq' % (sample_name, i))



# # # # # # # # # # # # # # # # # # # # # # # # # # #
# run PandaSeq, writing new files _ECBC-panda.fasta #
# # # # # # # # # # # # # # # # # # # # # # # # # # #

minimal_overlap = 10
for ig in igs:
	fwd = '%s_%s_ECBC-R1.fastq' % (sample_name, ig)
	rev = '%s_%s_R2.fastq' % (sample_name, ig)
	name = '%s_%s_ECBC-panda.fasta' % (sample_name, ig)
	cmd = 'pandaseq -o %s -f %s -r %s -w %s' % (minimal_overlap, fwd, rev, name)
	print(cmd)
	run_child(cmd)
	
	### delete ECBC-R1 and R2
	os.remove(fwd)
	os.remove(rev)

### count sequencs in fasta and fastq files
print('\n*** count after pandaseq ***\n')
print(run_child('seq_count.sh'))	



# # # # # # # # # # # # #
# collapse unique ECBC  #
# # # # # # # # # # # # #


for i in igs[0:1]: ################### ONLY TESTING IgG1 SAMPLE !!!!!!!
    fn = '%s_%s_ECBC-panda.fasta' % (sample_name, ig)
    itECBC_panda = FastqGeneralIterator(open(fn))

    ### stores all handles that will be used to write to files
    handle_dict = {}
    for g in groups:
	    handle_dict['%s_%s' % (i, g)] = open('%s_%s_%s.fastg' % (sample_name, i, g), 'w+')
        
    for it_obj in itECBC_panda:
        group = it_obj[0][1][:12].upper()
        print(group)

	### writing
	for read_here, r in zip(it_obj, reads):
		handle_here = handle_dict['%s_%s' % (i, g)]
		top = list(read_here)
		top[0] = '>%s' % top[0]
        top.append('')
        print(top)
        handle_here.write('\n'.join(top))

    ### close files
    for g in groups:
	    handle_dict['%s_%s' % (i, g)].close()

### count sequencs in fasta and fastq files
print('\n*** count after subtype assignment ***\n')
print(run_child('seq_count.sh'))























exit()
# # # # # # # # # # # # #
# collapse unique reads #
# # # # # # # # # # # # #

for ig in igs:
	#cmd = '''grep -v M0 %s_%s_ECBC-panda.fasta | sort | uniq -c | sort -nr | awk '$1>=10 {print ">"$1"\n"$2 }' > %s_%s_uniq.fasta''' % (sample_name, ig, sample_name, ig)
	cmd = '''grep -v M0 %s_%s_ECBC-panda.fasta | sort | uniq -c | sort -nr | awk '$1>=1 {print ">"$1 }' > %s_%s_uniq.fasta''' % (sample_name, ig, sample_name, ig)
	print(cmd)
	run_child(cmd)

### count sequencs after collapse unique reads
print('\n*** count after pandaseq ***\n')
print(run_child('seq_count.sh'))	







exit()

# # # # # # # # # # #
# run ab-analysis   #
# # # # # # # # # # #

run_child('ln -sv /data/AbX/bin/abanalysis/igblastn_linux')
run_child('ln -sv /data/AbX/bin/abanalysis/internal_data/')
run_child('ln -sv /data/AbX/bin/abanalysis/database/')

json_dir = './%s_ab-analysis' % sample_name
#for i in igs:
for i in igs[0:1]: ################### ONLY TESTING IgG1 SAMPLE !!!!!!!
	f = '%s_%s_ECBC-panda.fasta' % (sample_name, i)
	cmd = 'python2.7 /data/AbX/bin/abanalysis/ab_analysis.py -i %s -o %s -u %s' % (f, json_dir, ECBC_length)
	print(cmd)
	run_child(cmd)
	
	### delete .fastaprocessed directory
	cmd = 'rm -rf %s_%s_ECBC-panda.fastaprocessed' % (sample_name, i)
	run_child(cmd)











exit()
# # # # # # # # # # # # #
# import into Mongo DB  #
# # # # # # # # # # # # #

#su -c 'service mongod start'

MongoDB = '%s_MongoDB' % sample_name
MongoDB_log = '%s_MongoDB.log' % sample_name
cmd = 'python2.7 /data/AbX/bin/abanalysis/mongoimport.py -f %s -d %s -l %s' % (json_dir, MongoDB, log_file)
print(cmd)
run_child(cmd)



# # # # # # # # #
# run clonify   #
# # # # # # # # #


