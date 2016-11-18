#! /bin/bash
list=$(ls | grep "\.fastq$")
sum=0
for i in $list; do
	n=$(wc -l $i | cut -f 1 -d " ")
	n=$(($n / 4))
	sum=$(($sum + $n))
	echo -e $n "\t" $i 
done 
if (( $sum > 0 )); then
	echo -e $sum "\t Total"
	echo
fi

list=$(ls | grep "\.fasta$")
sum=0
for i in $list; do
	n=$(wc -l $i | cut -f 1 -d " ")
	n=$(($n / 2))
	sum=$(($sum + $n))
	echo -e $n "\t" $i 
done 
if (( $sum > 0 )); then
	echo -e $sum "\t Total"
	echo
fi