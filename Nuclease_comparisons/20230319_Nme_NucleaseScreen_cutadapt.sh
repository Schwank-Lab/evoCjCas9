#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_Nme_cutadaptlog.txt"
for  filename in *_R2_001.fastq.gz; do
	shortname="${filename::-16}"
	cutadapt -j 0 -g AAGGACGAAACACCG -o "Output/"$shortname"_5trim.fastq.gz" $filename >> $datum"_Nme_cutadaptlog.txt"
	cutadapt -j 0 -a GTTGTAGCTCCCTT -m 22 -M 22 --discard-untrimmed -o "Output/"$shortname"_Spacer.fastq.gz" "Output/"$shortname"_5trim.fastq.gz" >> $datum"_Nme_cutadaptlog.txt"
done
for  filename in *_R1_001.fastq.gz; do
	shortname="${filename::-16}"
	cutadapt -j 0 -g TCGTCGTAGCT --discard-untrimmed -o "Output/"$shortname"_3trim.fastq.gz" $filename >> $datum"_Nme_cutadaptlog.txt"
done
