# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 11:51:02 2021

Copy this file into "consensus" folder on workstation and run it with "python 2021013_Step_2_Rename_Consensus_Fasta_LS.py".
Also copy the reference sequence (including the primer sequences defined at the beginning) into the "consensus" folder as a .txt file.
Also copy the cluster.csv file from the previous "demultiplexer" script into "consensus" folder and adapt variable name below.
Create folder "allfasta" in "consensus" folder if it does not yet exist.
Final output are .fasta files which contain the sequences from the consensus folders but sorted by round. (stored in allfasta folder)


For alternative use, adapt variables in the top part.


@author: nimath
"""

from os import listdir
from Bio.Seq import Seq
from collections import Counter
import pandas as pd

path = "/media/localuser/LocalRAID/Directed_Evolution/Nanopore_20211006/data/AB20211006/p25115_26305_1/20211006_1503_X2_FAN51746_604779c9/20211010_all_finals/consensus/"

###### VARIABLES TO ADAPT FOR ALTERNATIVE USECASES ###

rounds = ['R1','R2','R3','R4','R5','R6','R7','R8','noround']
clusterbarcodename = '20211013_cluster_barcodes_LS.csv' # check the true name from previous Step 1 file and make sure you copied this file into the consensus folder

referencename = 'CjCas_Reference_Sequence.txt'  # TXT file which contains full cropped to the primers defined in previous script

domainlist = ['fullCjCas9']  # add more "domains" if necessary, but should work with only one big domain now (so it checks the whole sequence at once)
primer_fw_list = ['ATGGCAC'] # fw primer
primer_rev_list = ['CCATTAGC']  #reverse complement primer (5to3, lower strand), needs to be reverse complemented first on the next line to have it on the same strand as primer_fw_list for cropping
fw_end_seq_list = [str(Seq(x).reverse_complement()) for x in primer_rev_list]
sequencelength_list = [3241]  # length of reference sequence (see referencename)

#######


clusterbarcodedf = pd.read_csv(path+clusterbarcodename)

with open(referencename) as ref:
    referenceseq = ref.read()

with open(path+'allfasta/allsequences.fasta', 'w') as fp:
    pass

for rnd in rounds:
    with open(path+'allfasta/'+rnd+'sequences.fasta', 'w') as allfasta:
        pass


# consensusfolderlist = [f for f in listdir(path) if "consensus_cluster_" in f]  # if folders are named "consensus_cluster_XYZ"
consensusfolderlist = [f for f in listdir(path) if "consensus_assm_cluster_" in f]  # if folders are named "consensus_assm_cluster_XYZ_final"


#loop through domainlist; currently only one domain
for index, domain in enumerate(domainlist):
    i = 0
    counter = 0
    
    refseqcrop = primer_fw_list[index] + referenceseq.split(primer_fw_list[index])[1]
    refseqcrop = refseqcrop.split(fw_end_seq_list[index])[0] + fw_end_seq_list[index]
    
    print(domain)
    print(len(consensusfolderlist))
    lenlist = []
    for folder in consensusfolderlist:
        # clustername = folder[18:] # if folders are named "consensus_cluster_XYZ"
        clustername = folder[23:-6] # if folders are named "consensus_assm_cluster_XYZ_final"
        
        try:
            with open(path+folder+'/consensus.fasta') as fasta:
                content = fasta.readlines()
                sequence = content[1]
                if sequence.find(primer_rev_list[index]) == 1:
                    sequence = str(Seq(sequence).reverse_complement())
                try:
                    sequencecrop = primer_fw_list[index] + sequence.split(primer_fw_list[index])[1]
                    sequencecrop = sequencecrop.split(fw_end_seq_list[index])[0] + fw_end_seq_list[index]
                    # translation = str(Seq(sequencecrop).translate())
                    lenlist.append(len(sequencecrop))
                    print(f'Clustername: {clustername}, seqlength: {len(sequencecrop)}')
                except IndexError: #ignore folders where flanking sequences for cropping could not be found
                    # print('no flanking seq found')
                    continue
                # print(len(sequencecrop))
                
                if (len(sequencecrop) > sequencelength_list[index]-9) and (len(sequencecrop) < sequencelength_list[index]+9):
                    # print('test')
                    with open(path+'allfasta'+'/'+clustername+'.fasta','w') as finalfasta:
                        finalfasta.write(content[0])
                        finalfasta.write(sequencecrop)
                    with open(path+'allfasta/allsequences.fasta', 'a') as allfasta:
                        allfasta.write(">"+clustername+"\n")
                        allfasta.write(sequencecrop+"\n")
                    # print(clustername)
                    for rnd in rounds:  # create a .fasta file for every round (R1sequences.fasta, R2sequences.fasta etc.)
                        clusterbarcodelist = clusterbarcodedf[rnd]
                        clusterbarcodelist = clusterbarcodelist.dropna()
                        clusterbarcodelist = [int(x) for x in clusterbarcodelist]
                        if int(clustername) in clusterbarcodelist:
                            # print('success')
                            with open(path+'allfasta/'+rnd+'sequences.fasta', 'a') as allfasta:
                                allfasta.write(">"+clustername+"\n")
                                allfasta.write(sequencecrop+"\n")
                       

        except FileNotFoundError or IndexError: #ignore folders without consensus.fasta 
            # print('FileError')
            continue
        i+=1
        counter+=1
        if i == 100:
            i=0
            print(domain,counter)
    print('Total files: {}'.format(counter))
    print(Counter(lenlist))