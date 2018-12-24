#!/usr/bin/env python3
"""Simple test facility"""

import os
import pytest
#from Bio.Seq import Seq
#from Bio import SeqIO

from ECBC.annotate_first import check_clustering_and_trim

test_msa_fasta_01 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'test_sample_104_caccaaaaccaactgaaattc_msa.fasta')

test_msa_fasta_02 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'AK170_3_06_B_IgG2_IGHV3-23_IGHJ6_22_aaaaacacaagctaaaacaac_387_12501_msa.fasta')

test_msa_fasta_03 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 'AK170_3_06_B_IgG2_IGHV3-23_IGHJ6_22_aaacactagaacagccttagc_387_12541_msa.fasta')

assert os.path.exists(test_msa_fasta_01)

def test_check_clustering_and_trim():
    list_file_name = check_clustering_and_trim(test_msa_fasta_01)
    #print("I am in the test!")
    assert not list_file_name
    
    
    list_file_name = check_clustering_and_trim(test_msa_fasta_02)
    assert len(list_file_name) == 3
    
    list_file_name = check_clustering_and_trim(test_msa_fasta_03)
    assert len(list_file_name) == 4