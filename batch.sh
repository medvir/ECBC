#!/bin/bash


### IgG_klMA_ECBC_V2.py
samples=("AK170-1_S1_L001_R1_001.fastq.gz
AK170-2_S2_L001_R1_001.fastq.gz
AK170-3_S3_L001_R1_001.fastq.gz
AK170-4_S4_L001_R1_001.fastq.gz
AK170-5_S5_L001_R1_001.fastq.gz")

for s in $samples; do
	echo $s
	#/data/AbX/ECBC/IgG_klMA_ECBC_V2.py $s
done


### annotate_first.py
samples=("/data/AbX/experiments/161006/IMGT_download/AK170_1_S1_l_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_2_S2_l_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_3_S3_l_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_4_S4_l_ECBC_panda
/data/AbX/experiments/161006/IMGT_download/AK170_5_S5_l_ECBC_panda")

for s in $samples; do
	echo $s
	#/data/AbX/ECBC/annotate_first.py $s
done


### annotate_first_new.py
samples=("/data/AbX/experiments/161006/IMGT_download/AK170_1_S1_IgA_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_1_S1_IgG1_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_1_S1_IgG2_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_1_S1_IgG3_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_1_S1_IgG4_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_2_S2_IgA_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_2_S2_IgG1_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_2_S2_IgG2_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_2_S2_IgG3_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_2_S2_IgG4_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_3_S3_IgA_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_3_S3_IgG1_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_3_S3_IgG2_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_3_S3_IgG3_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_3_S3_IgG4_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_4_S4_IgA_ECBC_panda
/data/AbX/experiments/161006/IMGT_download/AK170_4_S4_IgG1_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_4_S4_IgG2_ECBC_panda
/data/AbX/experiments/161006/IMGT_download/AK170_4_S4_IgG3_ECBC_panda
/data/AbX/experiments/161006/IMGT_download/AK170_4_S4_IgG4_ECBC_panda
/data/AbX/experiments/161006/IMGT_download/AK170_5_S5_IgA_ECBC_panda
/data/AbX/experiments/161006/IMGT_download/AK170_5_S5_IgG1_ECBC_panda_fasta
/data/AbX/experiments/161006/IMGT_download/AK170_5_S5_IgG2_ECBC_panda
/data/AbX/experiments/161006/IMGT_download/AK170_5_S5_IgG3_ECBC_panda
/data/AbX/experiments/161006/IMGT_download/AK170_5_S5_IgG4_ECBC_panda")

for s in $samples; do
	echo $s
	/data/AbX/ECBC/annotate_first_new.py $s
done