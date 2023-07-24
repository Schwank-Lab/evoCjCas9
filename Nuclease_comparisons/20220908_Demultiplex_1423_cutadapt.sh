#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_1423_demultiplex_cutadaptlog.txt"
for  filename in *d3d5d7d10*R1.fastq.gz; do
	shortname="${filename::-11}"
	veryshortname="${filename:28:-12}"
	cutadapt -j 0 --action=none -a d3=GTAGCCGTGGAA -a d5=AAACGCGTGGAA -a d7=CTACATGTGGAA -a d10=TCGGGTGTGGAA -o Demultiplexed/$veryshortname"_"{name}"_R1.fastq.gz" -p Demultiplexed/$veryshortname"_"{name}"_R2.fastq.gz" $filename $shortname"R2.fastq.gz" >> $datum"_1423_demultiplex_cutadaptlog.txt"
done
