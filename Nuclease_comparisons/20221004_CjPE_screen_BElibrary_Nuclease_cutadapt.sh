#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_Nuclease_BELib_cutadaptlog.txt"
for  filename in *R2_001.fastq.gz; do
	shortname="${filename::-16}"
	cutadapt -j 0 -e 2 -g GAAACACC --discard-untrimmed -o "Output/"$shortname"_5trim.fastq.gz" $filename >> $datum"_Nuclease_BELib_cutadaptlog.txt"
done
for  filename in *R1_001.fastq.gz; do
	shortname="${filename::-16}"
	cutadapt -j 0 -e 2 -g GCCAAGCT --discard-untrimmed -o "Output/"$shortname"_3trim.fastq.gz" $filename >> $datum"_Nuclease_BELib_cutadaptlog.txt"
done
