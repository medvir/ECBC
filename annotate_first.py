#!/usr/bin/env python3

import sys
import pandas as pd

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
    '''First 8 nt of sequence are barcode'''
    bc = df['Sequence'].apply(lambda x: str(x)[:8])
    df = df.copy()
    df['barcode'] = bc
    return df


summary_file = '/data/AbX/experiments/161006/IMGT_download/AK170_1_S1_IgG1_ECBC_panda_a/1_Summary.txt'
nt_file = '/data/AbX/experiments/161006/IMGT_download/AK170_1_S1_IgG1_ECBC_panda_a/3_Nt-sequences.txt'
sel_cols = ['Sequence ID', 'V-GENE and allele', 'J-GENE and allele', 'Sequence',
            'CDR3-IMGT length']

imgt_ann = (pd.read_csv(filepath_or_buffer=summary_file, delimiter='\t',
            header=0, low_memory=False)
            .loc[:, sel_cols]
            .pipe(extract_genes)
            .pipe(extract_barcode)
)

imgt_seqs = (pd.read_csv(nt_file, delimiter='\t', header=0, low_memory=False)
            .loc[:, ['Sequence ID', 'V-D-J-REGION']]
)


grouped_all = (imgt_ann.pipe(pd.groupby, ('VGENE', 'JGENE', 'CDR3-IMGT length', 'barcode'))
                      .size().sort_values(ascending=False)
)

first_group = imgt_ann[(imgt_ann['VGENE'] == 'Homsap IGHV3-11') & \
                       (imgt_ann['JGENE'] == 'Homsap IGHJ5') & \
                       (imgt_ann['CDR3-IMGT length'] == '12') & \
                       (imgt_ann['barcode'] == 'agacagct')]

to_write = pd.merge(first_group, imgt_seqs, on='Sequence ID')

h = open('g.fasta', 'w')
for idx, s in to_write.iterrows():
    print('>%s' % s['Sequence ID'], file=h)
    print('%s' % s['V-D-J-REGION'], file=h)
h.close()
