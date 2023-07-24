# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 15:26:10 2021

Step 1 after generating clusters with medaka/pomoxy.
Copy script into folder where all final cluster .fasta files are stored. (143.fasta, 2003.fasta, etc)
Script runs through all cluster .fasta files and checks which of the barcodes occurs the most.
Then the cluster number is assigned to a round and this is stored in a .csv file in the same folder.
Use this .csv file for further downstream analysis (see file Step_2).

@author: nimath
"""

from os import listdir
import pandas as pd
from Bio.Seq import Seq

clusterbarcodefile = '20211013_cluster_barcodes_LS.csv'

rounds = ['R1','R2','R3','R4','R5','R6','R7','R8','noround']
barcodes_fwd = ['tttaattccaagagttcgtacccc', 'ttccccttagtgcctacttattct', 'gaagctattttaccacgcattgtt', 'gaaaaagctttgggcataaggtcc', 'acggatttagagtggtatctccct', 'cgaaggtttagatacgcaagtatc', 'gatgtgcacccggactctggccaa', 'ggtctgtgagtgccccacgtcaca']
barcodes_fwd = [x.upper() for x in barcodes_fwd]
barcodes_rev = [str(Seq(x).reverse_complement()) for x in barcodes_fwd]
barcodes_fwd.append('noround')
barcodes_rev.append('noround')

# make dictionary of FW/RV rounds:
barcodes_fwd_dic = {}
barcodes_rev_dic = {}
allround_dic = {}
for index, evoround in enumerate(rounds):
    barcodes_fwd_dic[evoround] = barcodes_fwd[index]
    barcodes_rev_dic[evoround] = barcodes_rev[index]
    allround_dic[evoround] = []  
    
filenamelist = [f for f in listdir() if ".fasta" in f]

i = 0
counter = 0
print(len(filenamelist))
for filename in filenamelist:
    with open(filename) as currentFile:
        try:
            clusternr = filename.split('cluster_')[1]
            clusternr = clusternr.split('.fasta')[0]
        except:
            print(filename)
        text = currentFile.read()
        barcodecountlist = []
        for index, evoround in enumerate(rounds):  # loop through all barcodes and count how many times a match can be found (fw and rv)
            # if (barcodes_fwd_dic[evoround] in text) and (barcodes_rev_dic[evoround] in text):
            if (barcodes_fwd_dic[evoround] in text):
                barcodecount = text.count(barcodes_fwd_dic[evoround])
                barcodecountlist.append(barcodecount)
            elif evoround == 'noround': ## if no barcode was found, append to noRound list
                barcodecountlist.append(1)
            else:
                barcodecountlist.append(0)
        # print(barcodecountlist)
        # print()
        evoround = rounds[barcodecountlist.index(max(barcodecountlist))]
        allround_dic[evoround].append(clusternr)
    i+=1
    counter+=1
    if i == 100:
        i=0
        print(counter)
print(allround_dic[evoround])

allrounddf = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in allround_dic.items() ]))


allrounddf.to_csv(clusterbarcodefile)