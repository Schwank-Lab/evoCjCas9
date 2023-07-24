#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_cutadaptlog.txt"
for  filename in *R1.fastq.gz; do
	shortname="${filename::-12}"
	cutadapt -j 0 -g AAGGACGAAACACCG -o "Output/"$shortname"_5trim.fastq.gz" $filename >> $datum"_1421_cutadaptlog.txt"
	cutadapt -j 0 -a GTTCTAGTCCCTGA -m 22 -M 22 --discard-untrimmed -o "Output/"$shortname"_Spacer.fastq.gz" "Output/"$shortname"_5trim.fastq.gz" >> $datum"_1421_cutadaptlog.txt"
done
for  filename in *R2.fastq.gz; do
	shortname="${filename::-12}"
	cutadapt -j 0 -g tcgtcgtagct --discard-untrimmed -o "Output/"$shortname"_3trim.fastq.gz" $filename >> $datum"_1421_cutadaptlog.txt"
done
