#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_cutadaptlog.txt"
for  filename in *R1.fastq.gz; do
	shortname="${filename::-12}"
	cutadapt -j 0 -g AAGGACGAAACACCG -o "Output/"$shortname"_5trim.fastq.gz" $filename >> $datum"_1422_cutadaptlog.txt"
	cutadapt -j 0 -a GTCTTAGTACTCTGG -m 22 -M 22 --discard-untrimmed -o "Output/"$shortname"_Spacer.fastq.gz" "Output/"$shortname"_5trim.fastq.gz" >> $datum"_1422_cutadaptlog.txt"
done
for  filename in *R2.fastq.gz; do
	shortname="${filename::-12}"
	cutadapt -j 0 -g ACCGCGTAGAC --discard-untrimmed -o "Output/"$shortname"_3trim.fastq.gz" $filename >> $datum"_1422_cutadaptlog.txt"
done
