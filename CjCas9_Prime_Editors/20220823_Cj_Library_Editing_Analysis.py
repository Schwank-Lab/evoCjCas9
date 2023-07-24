# -*- coding: utf-8 -*-
"""


@author: nimath
"""
import pandas as pd
from Bio import SeqIO
import gzip
from os import listdir
import numpy as np
import os
import time

# Get the current working directory
cwd = os.path.join(os.getcwd(), '')


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
    with gzip.open(filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        print('Start Protospacer lookup loop...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            protospacer = seq_record.seq
            protomatch = protolookup.get(str(protospacer)[:5])
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
    with gzip.open(filename, "rt") as fasta_file:
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


templatedforiginal = pd.read_csv('20220331_CjCas_PE_Library_Final_Table.csv')
templatedf = templatedforiginal.copy()



def list_files1(directory):
    return [f for f in listdir(directory) if ('_5trim' in f)]

# adapt path to the location of the fastq files (change to "/Fastq-Files/"" when downloading this repo):    
path = cwd + "Output/"
exactpath = cwd

filelist = list_files1(path)
filtered = ''

for filename in filelist:  # loop through all SAMPLES in a directory with "_5trim" in the name (listfiles1 function)
    print(filename)
    templatedf = templatedforiginal.copy()
    templatedf['protoshort'] = templatedf['Protospacer-Sequence'].apply(lambda x: x[:5])
    shortname = filename[0:-15]
    prototemplate = templatedf['protoshort'].values.flatten()
    uniquebarcodetemplate =templatedf['barcode'].values.flatten()

    protolookup, uniquebarcodelookup = lookup(prototemplate,uniquebarcodetemplate)
    del prototemplate, uniquebarcodetemplate
    diseasedict = importfiles(shortname,protolookup,uniquebarcodelookup)
    diseasedf = pd.DataFrame.from_dict(diseasedict,orient='index')  # make dataframe from dict
    del diseasedict

    initialdiseasedf = diseasedf.copy()
    diseasedf.dropna(subset = ['protomatch', 'barcodematch'], inplace=True) # do not drop elements without barcode match, since barcode match only includes barcodes from ambiguous sequences
    diseasedf['match'] = [list(set(a).intersection(set(b))) for a, b in zip(diseasedf.protomatch, diseasedf.barcodematch)]
    diseasedf['numbermatch'] = diseasedf['match'].apply(lambda x: len(x))
    final_pure_diseasedf = diseasedf.copy()
    final_pure_diseasedf = final_pure_diseasedf[final_pure_diseasedf['numbermatch'] == 1]
    
    print("calculating...")
    templatedf["editedcount"] = 0
    templatedf["uneditedcount"] = 0
    templatedf["indelcount"] = 0
    templatedf["nickindelcount"] = 0
    templatedf["beforeflapindelcount"] = 0
    templatedf["totalreads"] = 0
    filename = shortname+'_3trim.fastq.gz'
    templatenumpy = templatedf.to_numpy()
    
    readcounter = 0
    readcounttemp = 0
    start = time.time()
    
    
    # define column position for used column in numpy array:
    uneditedcountnr = templatedf.columns.get_loc("uneditedcount")
    editedcountnr = templatedf.columns.get_loc("editedcount")
    indelcountnr = templatedf.columns.get_loc("indelcount")
    nickindelcountnr = templatedf.columns.get_loc("nickindelcount")
    beforeflapindelcountnr = templatedf.columns.get_loc("beforeflapindelcount")
    ampliconnr = templatedf.columns.get_loc("wide_initial_target")
    mutated_ampliconnr = templatedf.columns.get_loc("wide_mutated_target")
    editingpositionnr = templatedf.columns.get_loc("Editing_Position")
    RTlengthpositionnr = templatedf.columns.get_loc("RTlength")
    RToverhanglengthpositionnr = templatedf.columns.get_loc("RToverhanglength")
    correction_typepositionnr = templatedf.columns.get_loc("Correction_Type")
    correction_lengthpositionnr = templatedf.columns.get_loc("Correction_Length")
    namenr = templatedf.columns.get_loc("name")
    
    with gzip.open(path+filename, "rt") as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            targetend = str(seq_record.seq.reverse_complement())[1:]  # start at pos 1 to shift into same frame as reference sequence
            
            identifier = seq_record.id
            readcounter+=1
            readcounttemp+=1
            
            startofwindow = 27  # bases from the left, 2 bases before nick
            
            if readcounttemp == 100000:
                print(readcounter)
                readcounttemp = 0
            if not identifier in final_pure_diseasedf.index:
                # if identifier in multiple_diseasedf.index:
                    
                continue
            variantindex = final_pure_diseasedf.at[identifier,'match'][0]
            RTlength = int(templatenumpy[variantindex,RTlengthpositionnr])
            RToverhanglength = int(templatenumpy[variantindex,RToverhanglengthpositionnr])
            correction_type = templatenumpy[variantindex,correction_typepositionnr]
            correction_length = int(templatenumpy[variantindex,correction_lengthpositionnr])
            
            # analyze window from 2bp left of nick and 5bp after flap
            if correction_type == 'Replacement':
                sequenceWT = targetend[startofwindow:startofwindow+4+3+RTlength]
                sequenceWT_beforeflap = targetend[startofwindow:startofwindow+RTlength]  # full sequence until 2bp before flap
                sequenceWT_nick = targetend[startofwindow:startofwindow+4]
                
                sequenceMUT = targetend[startofwindow:startofwindow+4+3+RTlength]
                sequenceMUT_beforeflap = targetend[startofwindow:startofwindow+RTlength]  # full sequence until 2bp before flap
                sequenceMUT_nick = targetend[startofwindow:startofwindow+4]
                
                
                controlWT = templatenumpy[variantindex,ampliconnr][startofwindow:startofwindow+4+3+RTlength]
                controlWT_beforeflap = templatenumpy[variantindex,ampliconnr][startofwindow:startofwindow+RTlength]
                controlWT_nick = templatenumpy[variantindex,ampliconnr][startofwindow:startofwindow+4]
                
                controlMUT = templatenumpy[variantindex,mutated_ampliconnr][startofwindow:startofwindow+4+3+RTlength]
                controlMUT_beforeflap = templatenumpy[variantindex,mutated_ampliconnr][startofwindow:startofwindow+RTlength]
                controlMUT_nick = templatenumpy[variantindex,mutated_ampliconnr][startofwindow:startofwindow+4]
                
                
                if sequenceWT == controlWT:
                    templatenumpy[variantindex,uneditedcountnr] += 1 # uneditedcount is column 42
                elif sequenceMUT == controlMUT:
                    templatenumpy[variantindex,editedcountnr] += 1
                else:

                    templatenumpy[variantindex,indelcountnr] += 1
                    if (sequenceWT_nick != controlWT_nick) and (sequenceMUT_nick != controlMUT_nick):  # check if 4bp window around nick has unintended edits
                        templatenumpy[variantindex,nickindelcountnr] += 1
                    
                    if (sequenceWT_beforeflap != controlWT_beforeflap) and (sequenceMUT_beforeflap != controlMUT_beforeflap):  # check if 4bp window around nick has unintended edits
                        templatenumpy[variantindex,beforeflapindelcountnr] += 1

                        
            elif correction_type == 'Deletion':
                sequenceWT = targetend[startofwindow:startofwindow+4+3+RTlength+correction_length]
                sequenceWT_beforeflap = targetend[startofwindow:startofwindow+RTlength+correction_length]  # full sequence until 2bp before flap
                sequenceWT_nick = targetend[startofwindow:startofwindow+4]
                
                controlWT = templatenumpy[variantindex,ampliconnr][startofwindow:startofwindow+4+3+RTlength+correction_length]
                controlWT_beforeflap = templatenumpy[variantindex,ampliconnr][startofwindow:startofwindow+RTlength+correction_length]
                controlWT_nick = templatenumpy[variantindex,ampliconnr][startofwindow:startofwindow+4]
                
                sequenceMUT = targetend[startofwindow:startofwindow+4+3+RTlength]
                sequenceMUT_beforeflap = targetend[startofwindow:startofwindow+RTlength]
                sequenceMUT_nick = targetend[startofwindow:startofwindow+4]
                
                controlMUT = templatenumpy[variantindex,mutated_ampliconnr][startofwindow:startofwindow+4+3+RTlength]
                controlMUT_beforeflap = templatenumpy[variantindex,mutated_ampliconnr][startofwindow:startofwindow+RTlength]
                controlMUT_nick = templatenumpy[variantindex,mutated_ampliconnr][startofwindow:startofwindow+4]
                
               
                if sequenceWT == controlWT:
                    templatenumpy[variantindex,uneditedcountnr] += 1 # uneditedcount is column 42
                elif sequenceMUT == controlMUT:
                    templatenumpy[variantindex,editedcountnr] += 1
                else:

                    templatenumpy[variantindex,indelcountnr] += 1
                    if (sequenceWT_nick != controlWT_nick) and (sequenceMUT_nick != controlMUT_nick):  # check if 4bp window around nick has unintended edits
                        templatenumpy[variantindex,nickindelcountnr] += 1
                    
                    if (sequenceWT_beforeflap != controlWT_beforeflap) and (sequenceMUT_beforeflap != controlMUT_beforeflap):  # check if 4bp window around nick has unintended edits
                        templatenumpy[variantindex,beforeflapindelcountnr] += 1
                        
            elif correction_type == 'Insertion':
                sequenceWT = targetend[startofwindow:startofwindow+4+3+RTlength-correction_length]
                sequenceWT_beforeflap = targetend[startofwindow:startofwindow+RTlength-correction_length]  # full sequence until 2bp before flap
                sequenceWT_nick = targetend[startofwindow:startofwindow+4]
                
                controlWT = templatenumpy[variantindex,ampliconnr][startofwindow:startofwindow+4+3+RTlength-correction_length]
                controlWT_beforeflap = templatenumpy[variantindex,ampliconnr][startofwindow:startofwindow+RTlength-correction_length]
                controlWT_nick = templatenumpy[variantindex,ampliconnr][startofwindow:startofwindow+4]
                
                sequenceMUT = targetend[startofwindow:startofwindow+4+3+RTlength]
                sequenceMUT_beforeflap = targetend[startofwindow:startofwindow+RTlength]
                sequenceMUT_nick = targetend[startofwindow:startofwindow+4]
                
                controlMUT = templatenumpy[variantindex,mutated_ampliconnr][startofwindow:startofwindow+4+3+RTlength]
                controlMUT_beforeflap = templatenumpy[variantindex,mutated_ampliconnr][startofwindow:startofwindow+RTlength]
                controlMUT_nick = templatenumpy[variantindex,mutated_ampliconnr][startofwindow:startofwindow+4]
                
                
            
               
                if sequenceWT == controlWT:
                    templatenumpy[variantindex,uneditedcountnr] += 1 # uneditedcount is column 42
                elif sequenceMUT == controlMUT:
                    templatenumpy[variantindex,editedcountnr] += 1
                else:
                    templatenumpy[variantindex,indelcountnr] += 1
                    if (sequenceWT_nick != controlWT_nick) and (sequenceMUT_nick != controlMUT_nick):  # check if 4bp window around nick has unintended edits
                        templatenumpy[variantindex,nickindelcountnr] += 1
                    
                    if (sequenceWT_beforeflap != controlWT_beforeflap) and (sequenceMUT_beforeflap != controlMUT_beforeflap):  # check if 4bp window around nick has unintended edits
                        templatenumpy[variantindex,beforeflapindelcountnr] += 1
                
    end = time.time()
    print('Time for loop:',end-start)

    print()
    
    #make again dataframe out of numpy array for easier saving as csv:
    templatedf = pd.DataFrame(data = templatenumpy, 
                      index = templatedf.index.tolist(), 
                      columns = templatedf.columns.tolist())
      
    
    # del final_diseasedf # remove final_diseasedf from memory
    totalreads = templatedf["uneditedcount"] + templatedf["editedcount"] + templatedf["indelcount"]
    templatedf['totalreads'] = totalreads
    templatedf['totalreads'].replace(0, np.nan, inplace=True)
    percentageedited = (templatedf["editedcount"]/templatedf['totalreads'])*100
    templatedf['percentageedited'] = percentageedited
    percentageunedited = (templatedf["uneditedcount"]/templatedf['totalreads'])*100
    templatedf['percentageunedited'] = percentageunedited
    percentageindels = (templatedf["indelcount"]/templatedf['totalreads'])*100
    templatedf['percentageindel'] = percentageindels
    print('percentageindel',templatedf['percentageindel'].mean())
    print('percentageedited',templatedf['percentageedited'].mean())
    print('percentageindelmedian',templatedf['percentageindel'].median())
    print('percentageeditedmedian',templatedf['percentageedited'].median())
    print()
    
    unintendfilteredtemplatedf = templatedf[templatedf['percentageindel'] < 20]
    print('removed {} pegRNAs due to too high background unintended editing'.format(len(templatedf)-len(unintendfilteredtemplatedf)))
    finalfilteredtemplatedf = unintendfilteredtemplatedf[unintendfilteredtemplatedf['percentageedited'] < 5]
    print('removed {} pegRNAs due to too high background editing'.format(len(unintendfilteredtemplatedf)-len(finalfilteredtemplatedf)))
    print('{} % of initial reads could be correctly assigned.'.format(round(len(final_pure_diseasedf)/len(initialdiseasedf)*100,2)))
    print()
    
    
    templatedf.to_csv('Output/Analysis/20220823_'+shortname+'_analysisdf_focused.csv')
