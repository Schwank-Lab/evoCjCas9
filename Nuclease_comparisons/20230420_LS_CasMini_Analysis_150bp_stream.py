# -*- coding: utf-8 -*-
"""

Python script for the analysis of NGS sequencing reads of the LS-CasMini library.

@author: nimath
"""
import pandas as pd
from Bio import SeqIO
import gzip
from os import listdir
import numpy as np
import os
import time
from Bio.Seq import Seq

# Get the current working directory
cwd = os.path.join(os.getcwd(), '')


def lookup(prototemplate, barcodetemplate):  # generate lookup dictionaries for all three relevant regions (protospacer, barcode)
    protolookup = {}
    for i, z in enumerate(prototemplate):
            protolookup.setdefault(z, []).append(i)
    
    barcodelookup = {}
    for i, z in enumerate(barcodetemplate):
            barcodelookup.setdefault(z, []).append(i)
    
    return protolookup, barcodelookup


def importfiles(filename, protolookup, barcodelookup):
    '''

    Parameters
    ----------
    filename : string
        Filename of the fasta file, containing the NGS reads to analyze.

    Returns
    -------
    diseasedict : dictionary
        Dictionary containing read identifier as key and protospacer/barcode/targetend/fwread as sub-keys (sequence as value of those) .
    diseasedf : pandas dataframe
        pandas dataframe made from diseasedict.
    diseasedf_filtered : pandas dataframe
        filtered diseasedf without whole sequence, but only sequence parts needed to identify index.
    full_adapter_read_percent : float
        Float which is the amount of sequences (from 0 to 1) which contain all subsequences without NaN.

    '''
    diseasedict = {}
    filename = shortname+'_spacer'+filtered+'.fastq.gz'
    with gzip.open(path+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        print('Start omegaRNA lookup loop...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            protospacer = seq_record.seq
            protomatch = protolookup.get(str(protospacer))
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
    print('Spacer done')

    filename = shortname+'_target'+filtered+'.fastq.gz'
    with gzip.open(path+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        print('Start barcode lookup loop...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            barcode = str(seq_record.seq.reverse_complement())[-6:]  # 6-bp barcode at the end of the read

            barcodematch = barcodelookup.get(barcode)
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


def barcodelookupfunc(barcodetemplate):
    ''' Convert series of barcode matches to lookup dictionary'''
    barcodelookup = {}
    for index, valuelist in barcodetemplate.items():
        tempdict = {}
        for value in valuelist:
            tempdict[value]=index
        barcodelookup = mergeDict(barcodelookup,tempdict)
        
    for key in barcodelookup:
        if type(barcodelookup[key]) == list:
            templist = list(map(int, str(barcodelookup[key]).replace("[","").strip("]").split(', ')))
            barcodelookup[key] = templist
        if type(barcodelookup[key]) == int:
            barcodelookup[key] = [barcodelookup[key]]
    return barcodelookup


# scaffold = 'GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'
templatedforiginal = pd.read_csv('CasMINI_Lib.csv')
# templatedforiginal = templatedforiginal.drop(['Unnamed: 0'], axis=1)  # remove column with old index
# templatedforiginal['finalname'] = templatedforiginal.apply(lambda x: x.Name+str(x.RToverhanglength) ,axis=1)
templatedf = templatedforiginal.copy()


def list_files1(directory):
    return [f for f in listdir(directory) if '_spacer' in f]

# adapt path to the location of the fastq files (change to "\\Fastq-Files\\"" when downloading this repo):    
path = cwd + "Output/"
filelist = list_files1(path)
filtered = ''

for filename in filelist:  # loop through all SAMPLES in a directory with "R1" in the name (listfiles1 function)
    print(filename)
    templatedf = templatedforiginal.copy()
    shortname = filename[0:-16]
    prototemplate = templatedf['spacer'].values.flatten()
    barcodetemplate =templatedf['Barcode'].values.flatten()

    protolookup, barcodelookup = lookup(prototemplate,barcodetemplate)
    del prototemplate, barcodetemplate
    diseasedict = importfiles(filename,protolookup,barcodelookup)
    diseasedf = pd.DataFrame.from_dict(diseasedict,orient='index')  # make dataframe from dict
    #del diseasedict

    diseasedf.dropna(subset = ['protomatch', 'barcodematch'], inplace=True)
    diseasedf['match'] = [list(set(a).intersection(set(b))) for a, b in zip(diseasedf.protomatch, diseasedf.barcodematch)]

    final_diseasedf = diseasedf[diseasedf['match'].map(lambda d: len(d)) == 1][['match']]  # all reads which have a match;
    final_diseasedf['match'] = final_diseasedf['match'].apply(lambda x: x[0]) #convert list to integer
    
    print(f"Skew-score: {np.std(final_diseasedf['match'].value_counts()) / np.mean(final_diseasedf['match'].value_counts()):.2f}")

    
    print("calculating...")
    templatedf["uneditedcount"] = 0
    templatedf["indelcount"] = 0
    templatedf["totalreads"] = 0
    filename = shortname+'_target'+filtered+'.fastq.gz'
    templatenumpy = templatedf.to_numpy()
    
    readcounter = 0
    readcounttemp = 0
    start = time.time()
    
    
    # define column position for used column in numpy array:
    uneditedcountnr = templatedf.columns.get_loc("uneditedcount")
    indelcountnr = templatedf.columns.get_loc("indelcount")
    targetlongnr = templatedf.columns.get_loc("target_long")
    namenr = templatedf.columns.get_loc("Identifier")
    
    indelreadsaveragequalitylist = []
    wt_editedreadsaveragequalitylist = []
    with gzip.open(path+filename, "rt") as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            if len(seq_record.letter_annotations["phred_quality"]) > 0:
                if sum(seq_record.letter_annotations["phred_quality"])/len(seq_record.letter_annotations["phred_quality"]) > 27.5:
                    readseq = str(seq_record.seq)
                    
                    identifier = seq_record.id
                    readcounter+=1
                    readcounttemp+=1
                                        
                    if readcounttemp == 100000:
                        print(readcounter)
                        readcounttemp = 0
                    if not identifier in final_diseasedf.index:
                            
                        continue
                    variantindex = int(final_diseasedf.at[identifier,'match'])

                    target_long = str(Seq(templatenumpy[variantindex,targetlongnr]).reverse_complement())
                    
                    nickposition = 39
                    window_beforenick = 5
                    window_afternick = 5
                    
                    targetwindow = target_long[nickposition-window_beforenick:nickposition+window_afternick]
                    
                    actualreadwindow = readseq[nickposition-window_beforenick:nickposition+window_afternick]
                    
                    if targetwindow == actualreadwindow:
                        templatenumpy[variantindex,uneditedcountnr] += 1
                    else:
                        templatenumpy[variantindex,indelcountnr] += 1
                    
                    
    end = time.time()
    print('Time for loop:',end-start)
    print()
    
    #make again dataframe out of numpy array for easier saving as csv:
    templatedf = pd.DataFrame(data = templatenumpy, 
                      index = templatedf.index.tolist(), 
                      columns = templatedf.columns.tolist())
      
    
    # del final_diseasedf # remove final_diseasedf from memory
    totalreads = templatedf["uneditedcount"] + templatedf["indelcount"]
    templatedf['totalreads'] = totalreads
    templatedf['totalreads'].replace(0, np.nan, inplace=True)
    percentageindels = (templatedf["indelcount"]/templatedf['totalreads'])*100
    templatedf['percentageindel'] = percentageindels
    percentageunedited = (templatedf["uneditedcount"]/templatedf['totalreads'])*100
    templatedf['percentageunedited'] = percentageunedited

    

    templatedf.to_csv('Output/Analysis/20230420_'+shortname+'_analysisdf_focused_qualityfilter.csv')