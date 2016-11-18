#! /bin/bash


list=$(ls | grep ECBC.fastq)

for i in $list; do
	echo $i
	seqtk seq -A $i | grep -v ">" | sort | uniq -c | sort -nr | awk '$1>=10 {print ">"$1"\n"$2 }'
done


exit


grep -B 1 AGACAGCTCATG AK170-1_S1.gz_IgG1_ECBC.fastq | grep "@M0"  | sed 's/@M0/M0/' > AK170_IgG1_reads.txt

