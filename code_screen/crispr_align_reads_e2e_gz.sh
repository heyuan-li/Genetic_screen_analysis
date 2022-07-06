#!/bin/sh
# Orgininal Author Nicholas J. Neill
# Adapted by Heyuan Li for CRISPR screen 200113

#$1 is library filename (no extension)
#$2 is path to sample key

#convert library to fasta format
Rscript scripts/make_library.r library/$1

#trim library
cutadapt -g TGGAAAGGACGAAAC -e .1 -O 8 -o library/$1_5trim.fasta library/$1.fasta > library/$1_trim_report_5.txt
cutadapt -a TAGAGCTAGAAATAG -e .1 -O 8 -o library/$1_trimmed.fasta library/$1_5trim.fasta > library/$1_trim_report_3.txt
cat library/$1_trim_report_5.txt library/$1_trim_report_3.txt > library/$1_trim_report.txt
rm library/$1_5trim.fasta library/$1_trim_report_5.txt library/$1_trim_report_3.txt

#call bowtie2 to generate index files
echo Building indexed library.
bowtie2-build -q library/$1_trimmed.fasta library/$1
bowtie2-build -q library/$1.fasta library/$1
echo Finished!

#align trimmed reads to index library
echo Aligning reads.
for i in $(cut -f 1 $2|tail -n +2);
do
	for j in $(ls fastq/processed_fastq/R1_${i}_trimmed.fastq.gz);
	do
		echo ${j}
		bowtie2 -p 16 -N 1 -L 9 --gbar 1 --min-score C,-20,0 --mp 2,2 --rdg 0,2 --rfg 0,2 -x library/$1 -U ${j} -S aligned_sam/${i}_aligned.sam
	done 2> aligned_sam/alignment_error_report.txt
done
echo Finished!
