#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 15:36:47 2021

Run this file from within consensus folder, after running files Step_1 and Step_2.

@author: nimath
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from diff_match_patch import diff_match_patch
from collections import Counter

referencename = 'CjCas_Reference_Sequence.txt'  # TXT file which contains full cropped to the primers defined in previous script
dateofrunning = '20211013'
shortname = 'LS'
startingaminoacidindex = -91 #standard is 0, adapt if the aminoacids at the beginning of the sequences should start negative
rounds = ['R1','R2','R3','R4','R5','R6','R7','R8','noround']


def aminoacidchangefunc(position, initialbp, mutatedbp):
    '''Show amino acid change from reference to query.'''
    mutposintriplet = position%3
    aastartpos = position - mutposintriplet
    originaltriplet = reference[aastartpos:aastartpos+3]
    mutatedtriplet = reference[aastartpos:aastartpos+mutposintriplet] + mutatedbp + reference[aastartpos + mutposintriplet+1:aastartpos+3]
    # print(originaltriplet)
    # print(mutatedtriplet)
    
    initialaa = str(Seq(originaltriplet).translate())
    mutatedaa = str(Seq(mutatedtriplet).translate())
    aaposition = int(aastartpos/3)
    return initialaa+str(aaposition+startingaminoacidindex)+mutatedaa
    
    


def differencestrings(ref,query):
    '''Compare reference to query and return list with mutations on nucleotide basis.'''
    dif = diff_match_patch()
    differencelist = dif.diff_main(ref,query)  # reference vs query
    
    mutations = []
    aachange = []
    aachange_nosilent = []
    positionref = 0
    positionquery = 0
    for index, stretch in enumerate(differencelist):
        if stretch[0] == 0:
            positionref+=len(stretch[1])
            positionquery+=len(stretch[1])
        elif stretch[0] == -1: 
            if differencelist[index+1][0] == 1:  # if the following tuple is insertion, convert deletion/insertion to replacement
                initialbp = str(stretch[1])
                mutatedbp = str(differencelist[index+1][1])
                mutation = str(positionref)+'_'+initialbp+'/'+mutatedbp
                mutations.append((positionref, 'Replacement', str(stretch[1]), str(differencelist[index+1][1]), mutation))
                aminoacidchange = aminoacidchangefunc(positionref, initialbp, mutatedbp)
                aachange.append(aminoacidchange)
                if aminoacidchange[0] != aminoacidchange[-1]:
                    aachange_nosilent.append(aminoacidchange)
                else:
                    print('aachange: ',aminoacidchange)
                positionref+=len(stretch[1])
                positionquery+=len(stretch[1])
            else: #deletion in query
                mutation = str(positionref)+'_'+str(stretch[1])+'/'+'-'
                mutations.append((positionref, 'Deletion', str(stretch[1]), '-', mutation))
                positionref+=len(stretch[1])
                
        elif stretch[0] == 1: #insertion in query
            if not differencelist[index-1][0] == -1: # skip if previous was deletion since then it was converted to replacement
                mutation = str(positionref)+'_'+'-'+'/'+str(stretch[1])
                mutations.append((positionref, 'Insertion', '-', str(stretch[1]),mutation))
    mutationshortlist = [x[4] for x in mutations]
    replacementonly = [x[4] for x in mutations if x[1] == 'Replacement']
    return mutations, mutationshortlist, replacementonly, aachange, aachange_nosilent


def fasta2df(infile):
    seq_df = pd.DataFrame()
    with open(infile) as fasta_file:  # Will close handle cleanly
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            seq_df.loc[int(seq_record.id),'sequence'] = str(seq_record.seq)
    return seq_df

with open(referencename) as ref:
    reference = ref.read()

# fasta = fasta2df('allfasta//allsequences.fasta')
fullaacountdf = []
fullaacombinationcountdf = []
fullmutationcountdf = []
fullmutationcombinationcountdf = []
fullaa_nosilentcombinationcountdf = []

for rnd in rounds:
    fasta = fasta2df('allfasta//'+rnd+'sequences.fasta')
    
    mutationdf = pd.DataFrame()
    mutationsummary = []
    aasummary = []
    aa_nosilentsummary = []
    for index, query in fasta.iterrows():
        mutations, mutationshortlist, replacementonly, aachange, aachange_nosilent = differencestrings(reference,query['sequence'])
        mutationsummary+=mutationshortlist
        aasummary+=aachange
        aa_nosilentsummary+=aachange_nosilent
        mutationdf.loc[index,'replacementmutations'] = str(replacementonly)
        mutationdf.loc[index,'mutationshortlist'] = str(mutationshortlist)
        mutationdf.loc[index,'mutations'] = str(mutations)
        mutationdf.loc[index,'aachangereplacement'] = str(aachange)
        mutationdf.loc[index,'aachangenosilentreplacement'] = str(aachange_nosilent)
        print(replacementonly)
        print()
    
    aacount = Counter(aasummary)
    aacountdf = pd.DataFrame.from_dict(aacount, orient='index').reset_index()
    aacountdf = aacountdf.set_index('index')
    fullaacountdf.append(aacountdf)
    
    aacombinationcount = Counter(mutationdf['aachangereplacement'])
    aacombinationcountdf = pd.DataFrame.from_dict(aacombinationcount, orient='index').reset_index()
    aacombinationcountdf = aacombinationcountdf.set_index('index')
    aacombinations = list(aacombinationcountdf.index)
    aacombinations = [x.replace('[','').replace(']','').replace(',','').replace("'","") for x in aacombinations]
    aacombinationcountdf['aacombination'] = aacombinations
    aacombinationcountdf = aacombinationcountdf.set_index('aacombination')

    fullaacombinationcountdf.append(aacombinationcountdf)    
    
    aa_nosilentcombinationcount = Counter(mutationdf['aachangenosilentreplacement'])
    aa_nosilentcombinationcountdf = pd.DataFrame.from_dict(aa_nosilentcombinationcount, orient='index').reset_index()
    aa_nosilentcombinationcountdf = aa_nosilentcombinationcountdf.set_index('index')
    aa_nosilentcombinations = list(aa_nosilentcombinationcountdf.index)
    aa_nosilentcombinations = [x.replace('[','').replace(']','').replace(',','').replace("'","") for x in aa_nosilentcombinations]
    aa_nosilentcombinationcountdf['aa_nosilentcombination'] = aa_nosilentcombinations
    aa_nosilentcombinationcountdf = aa_nosilentcombinationcountdf.set_index('aa_nosilentcombination')

    fullaa_nosilentcombinationcountdf.append(aa_nosilentcombinationcountdf)  
    
    
    
    mutationcount = Counter(mutationsummary)
    mutationcountdf = pd.DataFrame.from_dict(mutationcount, orient='index').reset_index()
    mutationcountdf = mutationcountdf.set_index('index')
    fullmutationcountdf.append(mutationcountdf) 
    
    mutationcombinationcount = Counter(mutationdf['mutationshortlist'])
    mutationcombinationcountdf = pd.DataFrame.from_dict(mutationcombinationcount, orient='index').reset_index()
    mutationcombinationcountdf = mutationcombinationcountdf.set_index('index')
    fullmutationcombinationcountdf.append(mutationcombinationcountdf) 
    

### collect all present mutations and combinations and make them in one dataframe before adding round per round
### then repeat analysis and add it to the prepared dataframe

mutationindices = []
mutationcombinationindices = []
aaindices = []
aacombinationindices = []
aanosilentcombinationindices = []
for ind, df in enumerate(fullmutationcountdf):
    mutationindices += list(fullmutationcountdf[ind].index)
    mutationcombinationindices += list(fullmutationcombinationcountdf[ind].index)
    aaindices += list(fullaacountdf[ind].index)
    aacombinationindices += list(fullaacombinationcountdf[ind].index)
    aanosilentcombinationindices += list(fullaa_nosilentcombinationcountdf[ind].index)
aaindices = list(set(aaindices))
aacombinationindices = list(set(aacombinationindices))
aanosilentcombinationindices = list(set(aanosilentcombinationindices))
mutationindices = list(set(mutationindices))
mutationcombinationindices = list(set(mutationcombinationindices))

fullaacountdf = pd.DataFrame()
fullaacountdf[1] = aaindices
fullaacountdf = fullaacountdf.set_index(1)
fullaacombinationcountdf = pd.DataFrame()
fullaacombinationcountdf[1] = aacombinationindices
fullaacombinationcountdf = fullaacombinationcountdf.set_index(1)
fullaa_nosilentcombinationcountdf = pd.DataFrame()
fullaa_nosilentcombinationcountdf[1] = aanosilentcombinationindices
fullaa_nosilentcombinationcountdf = fullaa_nosilentcombinationcountdf.set_index(1)
fullmutationcountdf = pd.DataFrame()
fullmutationcountdf[1] = mutationindices
fullmutationcountdf = fullmutationcountdf.set_index(1)
fullmutationcombinationcountdf = pd.DataFrame()
fullmutationcombinationcountdf[1] = mutationcombinationindices
fullmutationcombinationcountdf = fullmutationcombinationcountdf.set_index(1)


# repeat analysis and add to prepared dataframes with all indices already set
for rnd in rounds:
    fasta = fasta2df('allfasta//'+rnd+'sequences.fasta')
    
    mutationdf = pd.DataFrame()
    mutationsummary = []
    aasummary = []
    for index, query in fasta.iterrows():
        mutations, mutationshortlist, replacementonly, aachange, aachange_nosilent = differencestrings(reference,query['sequence'])
        mutationsummary+=mutationshortlist
        aasummary+=aachange
        aa_nosilentsummary+=aachange_nosilent
        mutationdf.loc[index,'replacementmutations'] = str(replacementonly)
        mutationdf.loc[index,'mutationshortlist'] = str(mutationshortlist)
        mutationdf.loc[index,'mutations'] = str(mutations)
        mutationdf.loc[index,'aachangereplacement'] = str(aachange)
        mutationdf.loc[index,'aachangenosilentreplacement'] = str(aachange_nosilent)
    
    aacount = Counter(aasummary)
    aacountdf = pd.DataFrame.from_dict(aacount, orient='index').reset_index()
    aacountdf = aacountdf.set_index('index')
    fullaacountdf[rnd] = aacountdf[0]
    
    aacombinationcount = Counter(mutationdf['aachangereplacement'])
    aacombinationcountdf = pd.DataFrame.from_dict(aacombinationcount, orient='index').reset_index()
    aacombinationcountdf = aacombinationcountdf.set_index('index')
    aacombinations = list(aacombinationcountdf.index)
    aacombinations = [x.replace('[','').replace(']','').replace(',','').replace("'","") for x in aacombinations]
    aacombinationcountdf['aacombination'] = aacombinations
    aacombinationcountdf = aacombinationcountdf.set_index('aacombination')
    fullaacombinationcountdf[rnd] = aacombinationcountdf[0]   
    
    aa_nosilentcombinationcount = Counter(mutationdf['aachangenosilentreplacement'])
    aa_nosilentcombinationcountdf = pd.DataFrame.from_dict(aa_nosilentcombinationcount, orient='index').reset_index()
    aa_nosilentcombinationcountdf = aa_nosilentcombinationcountdf.set_index('index')
    aa_nosilentcombinations = list(aa_nosilentcombinationcountdf.index)
    aa_nosilentcombinations = [x.replace('[','').replace(']','').replace(',','').replace("'","") for x in aa_nosilentcombinations]
    aa_nosilentcombinationcountdf['aa_nosilentcombination'] = aa_nosilentcombinations
    aa_nosilentcombinationcountdf = aa_nosilentcombinationcountdf.set_index('aa_nosilentcombination')

    fullaa_nosilentcombinationcountdf[rnd] = aa_nosilentcombinationcountdf[0]
    
    
    
    mutationcount = Counter(mutationsummary)
    mutationcountdf = pd.DataFrame.from_dict(mutationcount, orient='index').reset_index()
    mutationcountdf = mutationcountdf.set_index('index')
    fullmutationcountdf[rnd] = mutationcountdf[0]   
    
    mutationcombinationcount = Counter(mutationdf['mutationshortlist'])
    mutationcombinationcountdf = pd.DataFrame.from_dict(mutationcombinationcount, orient='index').reset_index()
    mutationcombinationcountdf = mutationcombinationcountdf.set_index('index') 
    fullmutationcombinationcountdf[rnd] = mutationcombinationcountdf[0]   


fullaacountdf = fullaacountdf.fillna(0)  # fill all NaN values with 0
fullaacountdf['Mutations'] = list(fullaacountdf.index)
fullaacountdf['Variants'] = ["Var_"+str(x) for x in range(1,len(fullaacountdf)+1)]
fullaacountdf['codonoptimization'] = fullaacountdf['Mutations'].apply(lambda x: True if x[0] == x[-1] else False)



fullaacombinationcountdf = fullaacombinationcountdf.fillna(0)
fullaacombinationcountdf['Mutations'] = list(fullaacombinationcountdf.index)
fullaacombinationcountdf['Variants'] = ["Var_"+str(x) for x in range(1,len(fullaacombinationcountdf)+1)]

fullaa_nosilentcombinationcountdf = fullaa_nosilentcombinationcountdf.fillna(0)
fullaa_nosilentcombinationcountdf['Mutations'] = list(fullaa_nosilentcombinationcountdf.index)
fullaa_nosilentcombinationcountdf['Variants'] = ["Var_"+str(x) for x in range(1,len(fullaa_nosilentcombinationcountdf)+1)]
fullaa_nosilentcombinationcountdf = fullaa_nosilentcombinationcountdf[~((fullaa_nosilentcombinationcountdf['R1']==0) & (fullaa_nosilentcombinationcountdf['R2']==0) & 
                                                                      (fullaa_nosilentcombinationcountdf['R3']==0) & (fullaa_nosilentcombinationcountdf['R4']==0) & 
                                                                      (fullaa_nosilentcombinationcountdf['R5']==0) & (fullaa_nosilentcombinationcountdf['R6']==0) & 
                                                                      (fullaa_nosilentcombinationcountdf['R7']==0) & (fullaa_nosilentcombinationcountdf['R8']==0) & 
                                                                      (fullaa_nosilentcombinationcountdf['noround']==0))]

fullmutationcountdf = fullmutationcountdf.fillna(0)
fullmutationcombinationcountdf = fullmutationcombinationcountdf.fillna(0)

fullaacountdf.to_csv(dateofrunning+'_'+shortname+'_summaryaacountdf.csv', index=False)
fullaacombinationcountdf.to_csv(dateofrunning+'_'+shortname+'_summaryaacombinationcountdf.csv', index=False)
fullaa_nosilentcombinationcountdf.to_csv(dateofrunning+'_'+shortname+'_summaryaanosilentcombinationcountdf.csv', index=False)
# fullmutationcountdf.to_csv(dateofrunning+'_'+shortname+'_summarymutationcountdf.csv')
# fullmutationcombinationcountdf.to_csv(dateofrunning+'_'+shortname+'_summarymutationcombinationcountdf.csv')