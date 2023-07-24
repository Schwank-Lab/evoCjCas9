#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_cutadaptlog.txt"
for  filename in *R1.fastq.gz; do
	shortname="${filename::-11}"
	veryshortname="${filename:20:-11}"
	cutadapt -j 0 --action=none -a final=GTTCTAGTCC -o Demultiplexed/$shortname""{name}"_R1.fastq.gz" -p Demultiplexed/$shortname""{name}"_R2.fastq.gz" $filename $shortname"R2.fastq.gz"
done
