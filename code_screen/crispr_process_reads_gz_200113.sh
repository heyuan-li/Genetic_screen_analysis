#!/bin/bash

# unzip and concatenate fastq files, then trim adapters using cutadapt
# Orgininal Author Nicholas J. Neill
# Adapted by Heyuan Li for CRISPR screen 200113

#$1 is path to sample key

cd /storage/westbrook/grail/westbrook/hl_crispr_minipool_screen/fastq/processed_fastq/

##concatenate fastq files

echo Merging fastq files.
for i in $(cut -f 1 $1|tail -n +2);
do
	cat ../raw_fastq/${i}* > ${i}_R1.fastq.gz
	echo ${i} merged.
done
echo Finished merging fastq files!


##trim adapters from reads
echo Started trimming 5 prime adapters.
for i in $(cut -f 1 $1|tail -n +2);
do
	echo ${i}
	for j in $(ls ${i}_R1.fastq.gz);
	do
		cutadapt -g TGGAAAGGACGAAAC -e .1 -O 8 -o 5trim_R1_${i}.fastq.gz ${j};
	done > 5trim_R1_${i}_report.txt;
done
echo Finished!

echo Started trimming 3 prime adapters.
for i in $(cut -f 1 $1|tail -n +2);
do
	echo ${i}
	for j in $(ls 5trim_R1_${i}.fastq.gz);
	do
		cutadapt -a TAGAGCTAGAAATAG -e .1 -O 8 -o 3trim_R1_${i}.fastq.gz ${j};
	done > 3trim_R1_"${i}"_report.txt;
done
echo Finished!

##organize files
echo Organizing results.
#cleanup fastq files
rm 5trim*.fastq.gz
for i in $(cut -f 1 $1|tail -n +2);
do
	mv 3trim_R1_${i}.fastq.gz R1_${i}_trimmed.fastq.gz;
done

for i in $(cut -f 1 $1|tail -n +2);
do
       echo ${i}
       awk '/Total reads processed/ {print $(NF)}
	    /Reads with adapters/ {print $(NF-1)}' 5trim_R1_${i}_report.txt
       awk '/Reads with adapters/ {print $(NF-1)}' 3trim_R1_${i}_report.txt;
done > global_trim_report.txt

#concatenate report files
for i in $(cut -f 1 $1|tail -n +2);
do
       for j in $(ls *_R1_${i}_report.txt);
       do
		cat ${j};
       done > R1_${i}_trim_report.txt;
done
#remove original reports
rm 5trim*report.txt
rm 3trim*report.txt
echo Finished!
