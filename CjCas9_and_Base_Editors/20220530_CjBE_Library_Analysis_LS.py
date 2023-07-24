# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:01:02 2022

@author: nimath
"""

from Bio import SeqIO
import gzip
from collections import Counter
from os import listdir
import numpy as np
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq


# Get the current working directory
cwd = os.path.join(os.getcwd(), '')
path = cwd


def lookup(prototemplate, uniquebarcodetemplate):  # generate lookup dictionaries for all protospacer and barcode
    protolookup = {}
    for i, z in enumerate(prototemplate):
            protolookup.setdefault(z, []).append(i)
    
    uniquebarcodelookup = {}
    for i, z in enumerate(uniquebarcodetemplate):
            uniquebarcodelookup.setdefault(z[-6:], []).append(i)
    
    return protolookup, uniquebarcodelookup


def importfiles(shortname, protolookup, uniquebarcodelookup):
    '''

    Parameters
    ----------
    filename : string
        Filename of the fasta file, containing the NGS reads to analyze.


    '''
    diseasedict = {}
    filename = exactpath+"Output/"+shortname+'_5trim'+'.fastq.gz'
    with gzip.open(path+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        print('Start Protospacer lookup loop...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            protospacer = seq_record.seq
            protomatch = protolookup.get(str(protospacer)[:7])
            identifier = seq_record.id
            if identifier in diseasedict:
                diseasedict[identifier]["protomatch"] = protomatch
            else:
                diseasedict[identifier] = {'protomatch':protomatch}
            count+=1
            counttemp+=1
            if counttemp == 1000000:
                print(count)
                counttemp = 0
    print('Protospacer done')
    
    filename = exactpath+"Output/"+shortname+'_3trim'+'.fastq.gz'
    with gzip.open(path+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        print('Start randombarcode lookup loop...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            barcode = str(seq_record.seq)[-6:]  # barcode is only the last 6 bases of the 3trim sequence
            barcodematch = uniquebarcodelookup.get(barcode)
            identifier = seq_record.id
            if identifier in diseasedict:
                diseasedict[identifier]["barcodematch"] = barcodematch

            count+=1
            counttemp+=1
            if counttemp == 1000000:
                print(count)
                counttemp = 0
    print('Barcode done')
    

    return diseasedict


def mergeDict(dict1, dict2):
        ''' Merge dictionaries and keep values of common keys in list'''
        dict3 = {**dict1, **dict2}
        for key, value in dict3.items():
            if key in dict1 and key in dict2:
                dict3[key] = [value , dict1[key]]
        return dict3

cjlibrary_templatedf = pd.read_csv('20220228_Library_March_combined_LS.csv')
cjlibrary_templatedf['Protospacer-Sequence_part'] = cjlibrary_templatedf['Protospacer-Sequence'].apply(lambda x: x[:7])
cjlibrary_templatedf['Protospacer_Length'] = cjlibrary_templatedf['Protospacer-Sequence'].apply(lambda x: len(x))
cjlibrary_templatedf['Protospacer-Sequence'] = cjlibrary_templatedf['Protospacer-Sequence'].apply(lambda x: x[-22:])
cjlibrary_templatedf['proto_extended'] = cjlibrary_templatedf['Target-SequenceReady'].apply(lambda x: str(Seq(x[-32:]).reverse_complement()))

cjeditornames = pd.read_csv('20220523_CjBE_library_names.csv', index_col='variant')

def list_files1(directory):
    return [f for f in listdir(directory) if ('_5trim' in f)]

exactpath = ''

filelist = list_files1('Output/')

shortnamelist = []

averageeditingdict = {}

for filename in filelist:  # loop through all SAMPLES in a directory with "protospacer" in the name (listfiles1 function)
    print(filename)
    shortname = filename[:-15]
    editorname = shortname
    print(editorname)
    editortype = cjeditornames.at[shortname,'editortype']
    print(editortype)
    prototemplate = cjlibrary_templatedf['Protospacer-Sequence_part'].values.flatten()
    uniquebarcodetemplate =cjlibrary_templatedf['6nt barcode'].values.flatten()

    protolookup, uniquebarcodelookup = lookup(prototemplate, uniquebarcodetemplate)
    del prototemplate, uniquebarcodetemplate
    diseasedict = importfiles(shortname,protolookup, uniquebarcodelookup)
    diseasedf = pd.DataFrame.from_dict(diseasedict,orient='index')  # make dataframe from dict
    del diseasedict
    
    replace_NaN = pd.isnull(diseasedf['barcodematch']) # make Series with True/False if barcodematch is NaN or not
    diseasedf.loc[replace_NaN,'barcodematch'] = -1 # replace all "NaN" with "-1" which will not match with others
    diseasedf['barcodematch'] = diseasedf['barcodematch'].apply(lambda x: [x] if type(x) == int else x) # put -1 into list for performing matching algorithm below
    diseasedf.dropna(subset = ['protomatch', 'barcodematch'], inplace=True)
    diseasedf['match'] = [list(set(a).intersection(set(b))) for a, b in zip(diseasedf.protomatch, diseasedf.barcodematch)]
    diseasedf['recombination'] = diseasedf.apply(lambda x: True if (x['protomatch'] != [-1]) and (x['barcodematch'] != [-1]) and (x['barcodematch'][0] not in x['protomatch']) else False,axis=1)
    
    
    final_diseasedf = diseasedf[diseasedf['match'].map(lambda d: len(d)) > 0][['match']]  # all reads which have at least one match; those with multiple matches will be assigned by barcodematch
    final_diseasedf['match'] = final_diseasedf['match'].apply(lambda x: x[0])
    
    
    for x in range(-10,22): # create empty columns which will be filled afterwards
        cjlibrary_templatedf[x+1] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
        cjlibrary_templatedf[str(x+1)+"_percent"] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
    
    # analyse different editings for ABE or CBE editors
    if (editortype == 'ABE'):
        for x in range(-10,22):
            cjlibrary_templatedf[str(x+1)+"_AtoG_editing"] = None
            cjlibrary_templatedf["AtoG_editing_list"] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
            cjlibrary_templatedf["A_positions"] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
    elif (editortype == 'CBE'):
        for x in range(-10,22):
            cjlibrary_templatedf[str(x+1)+"_CtoT_editing"] = None
            cjlibrary_templatedf["CtoT_editing_list"] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
            cjlibrary_templatedf["C_positions"] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
            
    elif (editortype == 'control') or (editortype == 'cas9'):
        for x in range(-10,22):
            cjlibrary_templatedf[str(x+1)+"_AtoG_editing"] = None
            cjlibrary_templatedf["AtoG_editing_list"] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
            cjlibrary_templatedf["A_positions"] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
            cjlibrary_templatedf[str(x+1)+"_CtoT_editing"] = None
            cjlibrary_templatedf["CtoT_editing_list"] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
            cjlibrary_templatedf["C_positions"] = np.empty((len(cjlibrary_templatedf), 0)).tolist()
        
        
    # Analyse fastqfiles by reading the reads into memory and evaluate editing based on matched library-variants in the importfiles function above.
    readsdict = {}
    targetbasedict = {}
    n = 4  # fastq file has one read for every 4 lines
    
    with gzip.open(exactpath + "Output/"+shortname+'_3trim.fastq.gz', 'rt') as fastqfile:
        lines = []
        for line in fastqfile:
            lines.append(line.rstrip())
            if len(lines) == n:
                identifier = lines[0].split()[0][1:]
                if not identifier in final_diseasedf.index:  # skip reads which have no match (because of recombination or sequencing errors etc)
                    lines = []
                    continue
                variantindex = diseasedf.loc[identifier].match[0] # get variantindex based on previous read matching
                raw_identifier = lines[0]
                
                target = str(Seq(lines[1][-38:-6]).reverse_complement()) # Target sequence within the sequenced reads

                
                # add each base of read to list in dataframe:
                for loc, base in enumerate(target):
                    cjlibrary_templatedf.loc[variantindex,loc+1-10].append(base)
                    
                # add target to dictionary "readsdict"; all reads from the same target are hereby combined
                if variantindex in readsdict:
                    readsdict[variantindex].append(target)
                else:
                    readsdict[variantindex] = [target]

                lines = []
    
    for entry in readsdict: # analyse all the reads for each variant in the library
        
        count = Counter(readsdict[entry])

        nrunchanged = count[cjlibrary_templatedf.loc[entry,'proto_extended']]
        
        cjlibrary_templatedf.loc[entry,'percent_unedited'] = (nrunchanged/len(readsdict[entry]))*100

        cjlibrary_templatedf.loc[entry,'nrreads'] = len(readsdict[entry])
        
        
        if (editortype == 'ABE') or (editortype == 'control') or (editortype == 'cas9'):  # ABE editing analysis:
            for loc in range(-10,22): # analyse each "A" in targetsequence/targetprotospacer
                # print(loc)
                cjlibrary_templatedf.loc[entry,str(loc+1)+"_percent"] = [Counter(cjlibrary_templatedf.loc[entry,loc+1])]
                if cjlibrary_templatedf.loc[entry,'proto_extended'][loc+10] == 'A':
                    numberofAtoG = cjlibrary_templatedf.loc[entry,str(loc+1)+"_percent"]['G']
                    try:
                        totalnumberofreads = len(cjlibrary_templatedf.loc[entry,loc+1])
                        cjlibrary_templatedf.loc[entry,str(loc+1)+"_AtoG_editing"] = numberofAtoG/totalnumberofreads*100
                        cjlibrary_templatedf.loc[entry,"AtoG_editing_list"].append(numberofAtoG/totalnumberofreads*100)
                        cjlibrary_templatedf.loc[entry,"A_positions"].append(loc+1)
                    except ZeroDivisionError:  # skip variants without any reads
                        cjlibrary_templatedf.loc[entry,str(loc+1)+"_AtoG_editing"] = None    
                else:
                    cjlibrary_templatedf.loc[entry,str(loc+1)+"_AtoG_editing"] = None
            
            averageAtoG = []
            for baseposition in range(-9,23):
                meanA = cjlibrary_templatedf[str(baseposition)+"_AtoG_editing"].mean()
                averageAtoG.append(meanA)
                
                
        if (editortype == 'CBE') or (editortype == 'control') or (editortype == 'cas9'):  # CBE editing analysis:
            for loc in range(-10,22): # analyse each "C" in targetsequence/targetprotospacer
                cjlibrary_templatedf.loc[entry,str(loc+1)+"_percent"] = [Counter(cjlibrary_templatedf.loc[entry,loc+1])]
                if cjlibrary_templatedf.loc[entry,'proto_extended'][loc+10] == 'C':
                    numberofCtoT = cjlibrary_templatedf.loc[entry,str(loc+1)+"_percent"]['T']
                    try:
                        totalnumberofreads = len(cjlibrary_templatedf.loc[entry,loc+1])
                        cjlibrary_templatedf.loc[entry,str(loc+1)+"_CtoT_editing"] = numberofCtoT/totalnumberofreads*100
                        cjlibrary_templatedf.loc[entry,"CtoT_editing_list"].append(numberofCtoT/totalnumberofreads*100)
                        cjlibrary_templatedf.loc[entry,"C_positions"].append(loc+1)
                    except ZeroDivisionError:  # skip variants without any reads
                        cjlibrary_templatedf.loc[entry,str(loc+1)+"_CtoT_editing"] = None    
                else:
                    cjlibrary_templatedf.loc[entry,str(loc+1)+"_CtoT_editing"] = None
            
            averageCtoT = []
            for baseposition in range(-9,23):
                meanC = cjlibrary_templatedf[str(baseposition)+"_CtoT_editing"].mean()
                averageCtoT.append(meanC)
            
    averageeditingdict[filename] = {'overall_modified':100-cjlibrary_templatedf['percent_unedited'].mean()}
    
    #store final dataframe containing sequence analysis
    cjlibrary_templatedf.to_csv('AnalysisFiles/20220530_CjBELibraryNova_dataframe_'+editorname+'.csv')
    
    print(len(final_diseasedf['match'].unique()))
    print('Total reads:',len(diseasedf))
    print('Successful assignment:',len(final_diseasedf))
    print('Recombined:',diseasedf.recombination.sum())
    print('Percentage recombined:',round((diseasedf.recombination.sum()/(diseasedf.recombination.sum()+len(final_diseasedf)))*100,2))
    