#!/usr/bin/env python

import os
import re
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
from time import localtime, strftime

def extract_genes(df):
    ''' Split V/J genes and alleles at the first *
    Homsap IGHV3-7*01 F, or Homsap IGHV3-7*02 F -> Homsap IGHV3-7
    '''
    vgene = df['V-GENE and allele'].apply(lambda x: "_or_".join(sorted(set(re.findall('IG.V.-[0-9]+', x)))))
    jgene = df['J-GENE and allele'].apply(lambda x: "_or_".join(sorted(set(re.findall('IG.J[0-9]+', x)))))
    df = df.copy()
    df['VGENE',] = vgene
    df['JGENE',] = jgene
    return df

def extract_barcode(df):
    '''First 21 nt of sequence are barcode'''
    bc = df['Sequence'].apply(lambda x: str(x)[:21])
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
chain = sys.argv[2]
if chain not in ['HC', 'LC']:
    sys.exit('Usage: %s folder_path [HC|LC]' % sys.argv[0])

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
if chain == 'HC':
    imgt_seqs = (pd.read_csv(nt_file, delimiter='\t', header=0, low_memory=True)
                .loc[:, ['Sequence ID', 'V-D-J-REGION']]
                )
    imgt_seqs.rename(index=str, columns={'V-D-J-REGION': 'Ab_sequence'}, inplace=True)
else:
    imgt_seqs = (pd.read_csv(nt_file, delimiter='\t', header=0, low_memory=True)
                .loc[:, ['Sequence ID', 'V-J-REGION']]
                )
    imgt_seqs.rename(columns={'V-J-REGION': 'Ab_sequence'}, inplace=True)

imgt_seqs = imgt_seqs.assign(nt_length = imgt_seqs['Ab_sequence'].str.len())
imgt_seqs.fillna(0, inplace=True)
imgt_seqs['nt_length'] = imgt_seqs['nt_length'].astype(int)

imgt_ann = pd.merge(imgt_seqs, imgt_ann, on='Sequence ID')

now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('%s grouping' % now, file=sys.stderr)
sys.stderr.flush()
grouped_all = (imgt_ann.pipe(pd.groupby, ('VGENE', 'JGENE', 'CDR3-IMGT length', 'barcode', 'nt_length')))

stem = folder_name.rstrip('/').split('/')[-1]
all_abs = '%s_all_antibodies.fasta' % stem
if os.path.exists(all_abs):
    os.remove(all_abs)
#n_groups = grouped_all.size().shape[0]

#n = 0
for idx, group in grouped_all:
#    n += 1
    seqs = group.shape[0]
    if seqs < 3:
        continue
    sys.stderr.flush()
    vgene, jgene = idx[:2]

    # heavy chain
    if chain == 'HC':
        if 'IGHV3-23' not in vgene or 'IGHJ6' not in jgene: ### heavy chain gene selection
            continue
    # light chain
    else:
        if 'IGLV2-23' not in vgene or ('IGLJ3' not in jgene and 'IGLJ5' not in jgene):  ### ligth chain gene selection
            continue

    seq_name = '_'.join([str(x) for x in idx])
    now = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print('%s: doing %s with %d sequences' % (now, seq_name, seqs), file=sys.stderr)
#    to_write = imgt_ann[(imgt_ann['VGENE'] == vgene) & \
#                        (imgt_ann['JGENE'] == jgene) & \
#                        (imgt_ann['CDR3-IMGT length'] == cdr_length) & \
#                        (imgt_ann['barcode'] == barcode) & \
#						(imgt_ann['nt_length'] == nt_len)]

    
    h = open('%s_g.fasta' % stem, 'w')
    for i, s in group.iterrows():
        print('>%s' % s['Sequence ID'], file=h)
        print('%s' % s['Ab_sequence'], file=h)
    h.close()
    run_child('muscle -in %s_g.fasta -out %s_msa.fasta' % (stem, stem))
    os.remove('%s_g.fasta' % stem)
    # zfill needs a higher value to sort correctly if > 999 reads per group
    run_child('cons -sequence %s_msa.fasta -name %s_reads-%s -outseq %s_tmp.fasta'
               % (stem, str(seqs).zfill(3), seq_name.replace(' ', '-'), stem))
    os.remove('%s_msa.fasta' % stem)
    run_child('cat %s_tmp.fasta >> %s' % (stem, all_abs))
    os.remove('%s_tmp.fasta' % stem)

abs_seqs = SeqIO.to_dict(SeqIO.parse(all_abs, 'fasta'))
new_abs_seqs = [abs_seqs[k] for k in sorted(abs_seqs.keys(), reverse=True)]
SeqIO.write(new_abs_seqs, 'sorted-%s' % all_abs, 'fasta')
