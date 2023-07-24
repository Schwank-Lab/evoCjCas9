# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:01:02 2022

@author: nimath
"""

from Bio import SeqIO
import gzip
from os import listdir
import numpy as np
import os
import pandas as pd
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
            barcode = str(seq_record.seq)[:6]  # barcode is only the last 6 bases of the 3trim sequence
            # print(barcode)
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

cjlibrary_templatedf = pd.read_csv('20230201_Cj_PAM_library_500.csv')
cjlibrary_templatedf['spacer_part'] = cjlibrary_templatedf['spacer'].apply(lambda x: x[:7])
cjlibrary_templatedf['Spacer_Length'] = cjlibrary_templatedf['spacer'].apply(lambda x: len(x))
#cjlibrary_templatedf['target'] = cjlibrary_templatedf['target_rv'].apply(lambda x: str(Seq(x).reverse_complement()))
# quantification window 10bp up and downstream of nick
cjlibrary_templatedf['target'] = cjlibrary_templatedf.apply(lambda x: str(Seq(x["Concat:"]).reverse_complement())[20+11+x['Spacer_Length']-3-10:20+11+x['Spacer_Length']-3+10],axis=1)
cjlibrary_templatedf['barcode_rv'] = cjlibrary_templatedf['Barcode'].apply(lambda x: str(Seq(x).reverse_complement()))



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
    prototemplate = cjlibrary_templatedf['spacer_part'].values.flatten()
    uniquebarcodetemplate =cjlibrary_templatedf['barcode_rv'].values.flatten()

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
        
        
    # Analyse fastqfiles by reading the reads into memory and evaluate editing based on matched library-variants in the importfiles function above.
    readsdict = {}
    targetbasedict = {}
    n = 4  # fastq file has one read for every 4 lines
    cjlibrary_templatedf["uneditedcount"] = 0
    cjlibrary_templatedf["indelcount"] = 0
    cjlibrary_templatedf["totalreads"] = 0

    templatenumpy = cjlibrary_templatedf.to_numpy()
    # define column position for used column in numpy array:
    uneditedcountnr = cjlibrary_templatedf.columns.get_loc("uneditedcount")
    editedcountnr = cjlibrary_templatedf.columns.get_loc("indelcount")
    libtargetnr = cjlibrary_templatedf.columns.get_loc("target")
    spacerlengthnr = cjlibrary_templatedf.columns.get_loc("Spacer_Length")

    namenr = cjlibrary_templatedf.columns.get_loc("Identifier")
    
    i = 0
    with gzip.open(exactpath + "Output/"+shortname+'_3trim.fastq.gz', 'rt') as fastqfile:
        print('Calculating editing...')
        for seq_record in SeqIO.parse(fastqfile, 'fastq'):
            identifier = seq_record.id
            if not identifier in final_diseasedf.index:  # skip reads which have no match (because of recombination or sequencing errors etc)
                lines = []
                continue
            variantindex = diseasedf.loc[identifier].match[0] # get variantindex based on previous read matching
            spacerlength = templatenumpy[variantindex,spacerlengthnr]

            str(seq_record.seq.reverse_complement())[1:]
            target = str(seq_record.seq)[11+spacerlength-3-10:11+spacerlength-3+10] # Target sequence within the sequenced reads 10bp up and downstream of nick
            #print(target)
            
            libraryseqtarget = templatenumpy[variantindex,libtargetnr]
            #print(libraryseqtarget)
            #print()

            if target == libraryseqtarget:
                templatenumpy[variantindex,uneditedcountnr] += 1 # uneditedcount is column 42
            else:
                templatenumpy[variantindex,editedcountnr] += 1
            i += 1
            if i % 100000 == 0:
                 print(i)
    
    #make again dataframe out of numpy array for easier saving as csv:
    cjlibrary_templatedf = pd.DataFrame(data = templatenumpy, 
                      index = cjlibrary_templatedf.index.tolist(), 
                      columns = cjlibrary_templatedf.columns.tolist())
    
    totalreads = cjlibrary_templatedf["uneditedcount"] +  cjlibrary_templatedf["indelcount"]
    cjlibrary_templatedf['totalreads'] = totalreads
    cjlibrary_templatedf['totalreads'].replace(0, np.nan, inplace=True)
    
    #store final dataframe containing sequence analysis
    cjlibrary_templatedf.to_csv('Output/Analysis/20230423_'+shortname+'.csv')
    
    print(len(final_diseasedf['match'].unique()))
    print('Total reads:',len(diseasedf))
    print('Successful assignment:',len(final_diseasedf))
    print('Recombined:',diseasedf.recombination.sum())
    print('Percentage recombined:',round((diseasedf.recombination.sum()/(diseasedf.recombination.sum()+len(final_diseasedf)))*100,2))
    