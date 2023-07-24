#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_cutadaptlog.txt"
for  filename in *R2_001.fastq.gz; do
	shortname="${filename::-25}"
	cutadapt -j 0 -g AAGGACGAAACACCG -o "Output/"$shortname"_5trim.fastq.gz" $filename >> $datum"_cutadaptlog.txt"
	cutadapt -j 0 -a GTTCTAGTCCCTGA --discard-untrimmed -o "Output/"$shortname"_Spacer.fastq.gz" "Output/"$shortname"_5trim.fastq.gz" >> $datum"_cutadaptlog.txt"
done
for  filename in *R1_001.fastq.gz; do
	shortname="${filename::-25}"
	cutadapt -j 0 -g tcgtcgtagct --discard-untrimmed -o "Output/"$shortname"_3trim.fastq.gz" $filename >> $datum"_cutadaptlog.txt"
done
