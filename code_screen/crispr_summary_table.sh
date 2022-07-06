#!/bin/bash

# Orgininal Author Nicholas J. Neill
# Adapted by Heyuan Li for CRISPR screen 200113

#$1 is path to sample key
#$2 is path to library (trimmed fasta)
#$3 is path to global trim_report

#extract read species counts and alignment score for each sample
echo Generating mapped species lists.
for i in $(cut -f 1 $1|tail -n +2);
do
	echo ${i}
	for j in $(ls aligned_sam/${i}_*sam);
	do
		awk '/YT:Z:/ && /AS:i:/' ${j}|cut -f 3,10,12
	done|sort|uniq -c|sort -k 2,2 -k 1,1nr -k 4,4n > summary_table/${i}_mapped_species.txt;
done
echo Finished!

echo Generating unmapped species lists.
for i in $(cut -f 1 $1|tail -n +2);
do
	echo ${i}
	for j in $(ls aligned_sam/${i}_*sam);
	do
		awk '/YT:Z:/ && !/AS:i:/' ${j}l|cut -f 10
	done|sort|uniq -c|sort -k 1nr > summary_table/${i}_unmapped_species.txt;
done
echo Finished!

#use R to build summary tables
echo Building summary tables.
Rscript script/summary_table.r $1 $2 $3
echo Finished!
