#!/bin/bash

samples=("AK170-1_S1_IgG1_ECBC-panda.fasta
AK170-1_S1_IgG4_ECBC-panda.fasta
AK170-1_S1_l_ECBC-panda.fasta")

for s in $samples; do
	echo $s
	grep -v M0127 $s | awk '{print substr ($0, 0, 12)}' | sort | uniq -c | sort -nr > ${s}.txt
	head -10 ${s}.txt | awk '{print $2}' > ${s}.indices
		
	cat ${s}.indices | while read i
	do
	   echo $i
	   grep -B 1 $i $s | grep -v -e '--' > ${s}_${i}.fasta
	done
	
	rm ${s}.indices
done

