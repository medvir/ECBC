#!/bin/bash

samples=("AK170-1_S1_L001_R1_001.fastq.gz
AK170-2_S2_L001_R1_001.fastq.gz
AK170-3_S3_L001_R1_001.fastq.gz
AK170-4_S4_L001_R1_001.fastq.gz
AK170-5_S5_L001_R1_001.fastq.gz")

for s in $samples; do
	/data/AbX/ECBC/IgG_klMA_ECBC_V2.py $s
done