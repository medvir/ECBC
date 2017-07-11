#!/bin/bash

smalt index -k 13 -s 2 11G5_HC_CDR3 /data/AbX/ECBC/Chunks/11G5_HC_CDR3.fasta
smalt map -n 28 -x -y 0.85 -f samsoft -o 5_IgA_HC_CDR3.sam 11G5_HC_CDR3 /data/AbX/ECBC/IMGT_grouped/S5_IgA/S5_IgA_ntCDR3.fasta
samtools view -Su 5_IgA_HC_CDR3.sam | samtools sort - 5_IgA_HC_CDR3
samtools index 5_IgA_HC_CDR3.bam

samtools flagstat 5_IgA_HC_CDR3.bam | grep "mapped ("


#samtools tview5_IgA_HCgl_CDR3.bam /data/AbX/ECBC/Chunks/11G5_HC_CDR3.fasta

#samtools view -f 45_IgA_HCgl_CDR3.bam > unmapped.sam

samtools view -b -F 4 5_IgA_HC_CDR3.bam > 5_IgA_HC_CDR3_mapped85.bam
samtools bam2fq -s 5_IgA_HC_CDR3_mapped85.fasta 5_IgA_HC_CDR3_mapped85.bam

grep '>' 5_IgA_HC_CDR3_mapped85.fasta | sed 's/>/ /' > 5_IgA_CDR3_mapped85.txt

#grep -e YYYYYYY -e YYFYYYY -e YYYYFYY /data/AbX/ECBC/IMGT_grouped/S1_IgA/1_Summary.txt > /data/AbX/ECBC/IMGT_grouped/1_IgA_Ymotivs.txt

#searching for barcods in original panda file
#grep -i -A1 '^cgctcacgtaaa' AK170-1_S1_IgG1_ECBC-panda.fasta | sed '/^-/ d' | wc -l

#searching for barcodes in IMGT output file
#grep 'agggaaacctacactcatcgc' 1_Summary.txt | grep -v unproductive | grep 'IGHV3-23' | grep 'IGHJ6' | cut -f2,29 | sed -e 's/^/>/' | sed -e 's/\t/\n/' > /data/AbX/ECBC/ECBC21/results/1_IgG3_agggaaacctacactcatcgc_6.fasta

#grep 'IGLV2-23' AK170_1_S1_l_ECBC21_comb/1_Summary.txt | grep -e 'IGLJ3' -e 'IGLJ5' | grep -v unproductive | awk -F'\t' '{if ($18 == 16) print $0 }' |wc -l

#grep 'IGLV2-23' AK170_1_S1_l_ECBC21_comb/1_Summary.txt | grep -e 'IGLJ3' -e 'IGLJ5' | grep -v unproductive | awk -F'\t' '{if ($18 >= 15 && $18 <= 17) print $0 }' |wc -l

#grep IGHJ6 AK170_1_S1_IgG1_ECBC21_comb/1_Summary.txt | grep IGHV3-23 | awk -F'\t' '{if ($18 == 21) print $0 }' | awk -F'\t' '{if ($3 ~ /^productive/) print $0 }' | wc -l