#! /bin/bash
# merge reads with pandaseq

### make list of files to analyse
if [ -d $1 ]; then list=$(ls $1 | grep .fastq); else list=$1; fi

### set minimal overlap = 1 if third argument not given
if [ -z $2 ]; then minimal=1; else minimal=$2; fi

for i in $list; do
	
	R1=$i
	R2=$(echo $R1 | sed 's/_R1_001.fastq/_R2_001.fastq/')

	### unzip fastq file if necessary
	if [[ $i =~ \.gz$ ]]
		then
			gunzip -c $R1 > R1.fastq
			gunzip -c $R2 > R2.fastq
			R1=$(echo $R1 | sed 's/.gz//')
			R2=$(echo $R2 | sed 's/.gz//')
		fi 
	
	### merge reads with pandaseq
	name=$(basename $R1 | sed 's/_L001_R1_001.fastq//')
	echo merging $R1 $R2 $minimal $name
	pandaseq -o $minimal -f $R1 -r $R2 -w ${name}_panda.fasta	

done