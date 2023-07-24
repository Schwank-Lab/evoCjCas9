#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_cutadaptlog.txt"
for  filename in *R2_001.fastq.gz; do
	shortname="${filename::-16}"
	cutadapt -j 0 -e 2 -g GAAACACC -o "Output/"$shortname"_5trim.fastq.gz" $filename >> $datum"_cutadaptlog.txt"
	#cutadapt -j 0 -e 2 -a GTTCTAG -m 22 -M 23 --discard-untrimmed -o "Output/"$shortname"_Protospacer.fastq.gz" "Output/"$shortname"_5trim.fastq.gz" >> $datum"_cutadaptlog.txt"
done
for  filename in *R1_001.fastq.gz; do
	shortname="${filename::-16}"
	revcompname=$shortname"_REVCOMP_001.fastq.gz"
	seqkit seq -r -p $filename | gzip -c > $revcompname
	cutadapt -j 0 -e 2 -a AGCTTGGC --discard-untrimmed -o "Output/"$shortname"_3trim.fastq.gz" $revcompname >> $datum"_cutadaptlog.txt"
	#cutadapt -j 0 -e 2 -g CTAGATCT --discard-untrimmed -o "Output/"$shortname"_randombarcode.fastq.gz" $revcompname >> $datum"_cutadaptlog.txt"
done
