#!/usr/bin/env python3

import os
import sys
import subprocess
import pandas as pd
from time import localtime, strftime

def extract_genes(df):
    '''Split V/J genes and alleles at the first *
    Homsap IGHV3-7*01 F, or Homsap IGHV3-7*02 F -> Homsap IGHV3-7
    '''
    vgene = df['V-GENE and allele'].apply(lambda x: str(x).split('*')[0])
    jgene = df['J-GENE and allele'].apply(lambda x: str(x).split('*')[0])
    df = df.copy()
    df['VGENE',] = vgene
    df['JGENE',] = jgene
    return df

def extract_barcode(df):
    '''First 12 nt of sequence are barcode'''
    bc = df['Sequence'].apply(lambda x: str(x)[:12])
    df = df.copy()
    df['barcode'] = bc
    return df

def run_child(cmd, exe='/bin/bash'):
    '''use subrocess.check_output to run an external program with arguments'''
    try:
        output = subprocess.check_output(cmd, universal_newlines=True,
        shell=True,
#        executable=exe,
        stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as ee:
        sys.exit("Execution of %s failed with returncode %d: %s" % (cmd, ee.returncode, ee.output))
        sys.exit(cmd)
        output = None
    return output

folder_name = sys.argv[1]
#folder_name = '/data/AbX/experiments/161006/IMGT_download/AK170_1_S1_IgG1_ECBC_panda_a'
summary_file = os.path.join(folder_name, '1_Summary.txt')
nt_file = os.path.join(folder_name, '3_Nt-sequences.txt')
sel_cols = ['Sequence ID', 'V-GENE and allele', 'J-GENE and allele', 'Sequence',
            'CDR3-IMGT length']

now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('%s reading summary file' % now, file=sys.stderr)
sys.stderr.flush()
imgt_ann = pd.read_csv(filepath_or_buffer=summary_file, delimiter='\t',
            header=0, low_memory=True)
mask = imgt_ann.Functionality.str.startswith('productive')

imgt_ann = (imgt_ann.loc[mask, sel_cols]
            .pipe(extract_genes)
            .pipe(extract_barcode)
)

now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('%s reading nt file' % now, file=sys.stderr)
sys.stderr.flush()
imgt_seqs = (pd.read_csv(nt_file, delimiter='\t', header=0, low_memory=True)
            .loc[:, ['Sequence ID', 'V-D-J-REGION']]
)

now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('%s grouping' % now, file=sys.stderr)
sys.stderr.flush()
grouped_all = (imgt_ann.pipe(pd.groupby, ('VGENE', 'JGENE', 'CDR3-IMGT length', 'barcode'))
                      .size().sort_values(ascending=False)
)

stem = folder_name.rstrip('/').split('/')[-1]
all_abs = '%s_all_antibodies.fasta' % stem
if os.path.exists(all_abs):
    os.remove(all_abs)

for idx, seqs in grouped_all.iteritems():
    if seqs < 3:
        break
    sys.stderr.flush()
    vgene, jgene, cdr_length, barcode = idx
    if vgene != 'Homsap IGHV3-23' or jgene != 'Homsap IGHJ6':
        continue
    seq_name = '_'.join(idx)
    now = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print('%s: doing %s with %d sequences' % (now, seq_name, seqs), file=sys.stderr)
    first_group = imgt_ann[(imgt_ann['VGENE'] == vgene) & \
                           (imgt_ann['JGENE'] == jgene) & \
                           (imgt_ann['CDR3-IMGT length'] == cdr_length) & \
                           (imgt_ann['barcode'] == barcode)]

    to_write = pd.merge(first_group, imgt_seqs, on='Sequence ID')

    h = open('g.fasta', 'w')
    for i, s in to_write.iterrows():
        print('>%s' % s['Sequence ID'], file=h)
        print('%s' % s['V-D-J-REGION'], file=h)
    h.close()
    run_child('muscle -in g.fasta -out msa.fasta')
    os.remove('g.fasta')
    run_child('cons -sequence msa.fasta -name %s-%d_reads -outseq tmp.fasta'
               % (seq_name.replace(' ', '-'), seqs))
    os.remove('msa.fasta')
    run_child('cat tmp.fasta >> %s' % all_abs)
    os.remove('tmp.fasta')
