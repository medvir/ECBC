#! /bin/bash
# merge reads with pandaseq

### make list of files to analyse
if [ -d $1 ]; then files=$(find $1 | grep _R1_001.fastq); else files=$1; fi

### set minimal overlap = 10 if third argument not given
if [ -z $2 ]; then min_overlap=10; else min_overlap=$2; fi

for i in $files; do
	
	R1=$i
	R2=$(echo $R1 | sed 's/_R1_001.fastq/_R2_001.fastq/')
	name=$(basename $R1 | sed 's/_L001_R1_001.fastq.gz//' | sed 's/_L001_R1_001.fastq//')
	echo merging $R1 $R2 $min_overlap $name
	
	### unzip fastq file if necessary
	if [[ $i =~ \.gz$ ]]
		then
			gunzip -c $R1 > ${name}_R1.fastq
			gunzip -c $R2 > ${name}_R2.fastq
			
			### merge reads with pandaseq
			pandaseq -o $min_overlap -f ${name}_R1.fastq -r ${name}_R2.fastq -w ${name}_panda.fasta
			#rm ${name}_R1.fastq
			#rm ${name}_R1.fastq
		else
			### merge reads with pandaseq
			pandaseq -o $min_overlap -f $R1 -r $R2 -w ${name}_panda.fasta	
		fi
done