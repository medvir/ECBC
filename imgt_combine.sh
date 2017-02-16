#!/bin/bash

# Combine IMGT two analysis files from the same sample that had to be analysed in increments of 500'000 sequences
# Usage `imgt_combine.sh directory_a directory_b joined_directory`
# Does not work for 11_Parameters.txt, as the headers are different.

dir_a=$1
dir_b=$2
dir_out=$3

files=$(ls $dir_a)
mkdir -p $dir_out

for file in $files; do
	head_a=$(head -1 $dir_a/$file)
	head_b=$(head -1 $dir_b/$file)
	if [ "$head_a" = "$head_b" ] ; then
		cat $dir_a/$file > $dir_out/$file
		tail -n +2 -q $dir_b/$file >> $dir_out/$file
		echo OK $dir_a/$file and $dir_b/$file into $dir_out/$file
	else
		echo Failed to join $dir_a/$file and $dir_b/$file
	fi
done

### remove directories
#rm -rf $dir_a
#rm -rf $dir_b