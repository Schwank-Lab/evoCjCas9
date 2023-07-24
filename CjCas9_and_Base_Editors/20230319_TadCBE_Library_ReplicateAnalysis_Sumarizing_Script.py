# -*- coding: utf-8 -*-
"""
Created on Wed May 25 12:40:24 2022
Analyze replicate analysis files from CjBE library screens.


@author: nimath
"""

import pandas as pd
from os import listdir

def list_files1(directory):
    return [f for f in listdir(directory) if ('.csv' in f) and ('20221003' in f)]

filelist = list_files1('AnalysisFiles/')

sampleoverview = pd.read_csv('20230319_TadCBE_sampleoverview.csv')
controloverviewdf = sampleoverview[sampleoverview['variant'] == 'ctrl']

# fileprefix = '20220530_CjBELibraryNova_dataframe_LS_'
fileprefix = '20230314_CjBE_dataframe_LS-'


CtoTcolumnlist = ['-9_CtoT_editing',
  '-8_CtoT_editing',
  '-7_CtoT_editing',
  '-6_CtoT_editing',
  '-5_CtoT_editing',
  '-4_CtoT_editing',
  '-3_CtoT_editing',
  '-2_CtoT_editing',
  '-1_CtoT_editing',
  '0_CtoT_editing',
  '1_CtoT_editing',
  '2_CtoT_editing',
  '3_CtoT_editing',
  '4_CtoT_editing',
  '5_CtoT_editing',
  '6_CtoT_editing',
  '7_CtoT_editing',
  '8_CtoT_editing',
  '9_CtoT_editing',
  '10_CtoT_editing',
  '11_CtoT_editing',
  '12_CtoT_editing',
  '13_CtoT_editing',
  '14_CtoT_editing',
  '15_CtoT_editing',
  '16_CtoT_editing',
  '17_CtoT_editing',
  '18_CtoT_editing',
  '19_CtoT_editing',
  '20_CtoT_editing',
  '21_CtoT_editing',
  '22_CtoT_editing']

AtoGcolumnlist = ['-9_AtoG_editing',
  '-8_AtoG_editing',
  '-7_AtoG_editing',
  '-6_AtoG_editing',
  '-5_AtoG_editing',
  '-4_AtoG_editing',
  '-3_AtoG_editing',
  '-2_AtoG_editing',
  '-1_AtoG_editing',
  '0_AtoG_editing',
  '1_AtoG_editing',
  '2_AtoG_editing',
  '3_AtoG_editing',
  '4_AtoG_editing',
  '5_AtoG_editing',
  '6_AtoG_editing',
  '7_AtoG_editing',
  '8_AtoG_editing',
  '9_AtoG_editing',
  '10_AtoG_editing',
  '11_AtoG_editing',
  '12_AtoG_editing',
  '13_AtoG_editing',
  '14_AtoG_editing',
  '15_AtoG_editing',
  '16_AtoG_editing',
  '17_AtoG_editing',
  '18_AtoG_editing',
  '19_AtoG_editing',
  '20_AtoG_editing',
  '21_AtoG_editing',
  '22_AtoG_editing']



collist_general = ['Unnamed: 0',
 'Unnamed: 0.1',
 'Name',
 'Gene',
 'Phenotype',
 'PAM',
 'ABE pos',
 'CBE pos',
 'Sequence',
 'Reverse-Complement',
 'Mutated-FW-Sequences',
 'Mutated-RV-Sequences',
 'TAAC Target-Strand',
 'ReferenceAllele',
 'TAAC PAM-Position',
 "5'Overhang",
 'Protospacer-Sequence',
 'Scaffold',
 'Target-SequenceReady',
 'TAAC Halffullblock-Sequence',
 '6nt barcode',
 'threeprimeoverhang',
 'Unnamed: 21',
 'FinalLibrary 01.03.2022',
 'Unnamed: 23',
 'Protospacer-Sequence_part',
 'Protospacer_Length']

# calculate average values for control replicates
for controlsample in controloverviewdf.variant:
    # if sample == 'control': # skip control here since it is done before already
    #     continue
    rep1 = pd.read_csv('AnalysisFiles/'+fileprefix+controlsample+'-R1-all.csv')
    rep2 = pd.read_csv('AnalysisFiles/'+fileprefix+controlsample+'-R2-all.csv')
    rep3 = pd.read_csv('AnalysisFiles/'+fileprefix+controlsample+'-R3-all.csv')
    
    controldf = rep1[collist_general].copy()
    for column in CtoTcolumnlist:
        try:
            controldf[column+'_rep1'] = rep1[column]
            controldf[column+'_rep2'] = rep2[column]
            controldf[column+'_rep3'] = rep3[column]
        except KeyError:  # skip this if for this editor there is no CtoT 
            break
            
        
    for column in AtoGcolumnlist:
        try:
            controldf[column+'_rep1'] = rep1[column]
            controldf[column+'_rep2'] = rep2[column]
            controldf[column+'_rep3'] = rep3[column]
        except KeyError:  # skip this if for this editor there is no AtoG
            break
        
    controldf['percent_unedited'+'_rep1'] = rep1['percent_unedited']
    controldf['percent_unedited'+'_rep2'] = rep2['percent_unedited']
    controldf['percent_unedited'+'_rep3'] = rep3['percent_unedited']
        
    controldf['nrreads'+'_rep1'] = rep1['nrreads']
    controldf['nrreads'+'_rep2'] = rep2['nrreads']
    controldf['nrreads'+'_rep3'] = rep3['nrreads']
    
    for column in CtoTcolumnlist:
        try:
            controldf[column+'_average_raw'] = controldf[[column+'_rep1',column+'_rep2',column+'_rep3']].mean(axis=1)
        except KeyError:  # skip this if for this editor there is no CtoT 
            break
    for column in AtoGcolumnlist:
        try:
            controldf[column+'_average_raw'] = controldf[[column+'_rep1',column+'_rep2',column+'_rep3']].mean(axis=1)
        except KeyError:  # skip this if for this editor there is no AtoG
            break
    controldf['percent_unedited'+'_average_raw'] = controldf[['percent_unedited'+'_rep1','percent_unedited'+'_rep2','percent_unedited'+'_rep3']].mean(axis=1)
    controldf['percent_modified'+'_average_raw'] = controldf['percent_unedited'+'_average_raw'].apply(lambda x: 100-x)
    
    for index, row in controldf.iterrows():
        highcoveragelist = []
        for readnr in ['nrreads_rep1','nrreads_rep2','nrreads_rep3']:
            if row[readnr] > 100:
                highcoveragelist.append(1)
            else:
                highcoveragelist.append(0)
        controldf.at[index,'highcoveragereps'] = sum(highcoveragelist)

controldf.to_csv('20230319_Cj_TadCBE_control_analysis_average.csv')





for sample in sampleoverview.variant:
    # if sample == 'control': # skip control here since it is done before already
    #     continue
    rep1 = pd.read_csv('AnalysisFiles/'+fileprefix+sample+'-R1-all.csv')
    rep2 = pd.read_csv('AnalysisFiles/'+fileprefix+sample+'-R2-all.csv')
    rep3 = pd.read_csv('AnalysisFiles/'+fileprefix+sample+'-R3-all.csv')
    
    summarydf = rep1[collist_general].copy()
    for column in CtoTcolumnlist:
        try:
            summarydf[column+'_rep1'] = rep1[column]
            summarydf[column+'_rep2'] = rep2[column]
            summarydf[column+'_rep3'] = rep3[column]
        except KeyError:  # skip this if for this editor there is no CtoT 
            break
            
        
    for column in AtoGcolumnlist:
        try:
            summarydf[column+'_rep1'] = rep1[column]
            summarydf[column+'_rep2'] = rep2[column]
            summarydf[column+'_rep3'] = rep3[column]
        except KeyError:  # skip this if for this editor there is no AtoG
            break
        
    summarydf['percent_unedited'+'_rep1'] = rep1['percent_unedited']
    summarydf['percent_unedited'+'_rep2'] = rep2['percent_unedited']
    summarydf['percent_unedited'+'_rep3'] = rep3['percent_unedited']
        
    summarydf['nrreads'+'_rep1'] = rep1['nrreads']
    summarydf['nrreads'+'_rep2'] = rep2['nrreads']
    summarydf['nrreads'+'_rep3'] = rep3['nrreads']
    
    for column in CtoTcolumnlist:
        try:
            summarydf[column+'_average_raw'] = summarydf[[column+'_rep1',column+'_rep2',column+'_rep3']].mean(axis=1)
            for index, value in summarydf[column+'_average_raw'].items():
                # if value in controls is larger than 20%, set to None (ignore for following analysis)
                if controldf.at[index,column+'_average_raw'] > 20:
                    summarydf.at[index,column+'_average_ctrladjusted'] = None
                else:
                    summarydf.at[index,column+'_average_ctrladjusted'] = (value-controldf.at[index,column+'_average_raw'])/(100-controldf.at[index,column+'_average_raw'])*100
            for replicatename in ['_rep1','_rep2','_rep3']: 
                for index, value in summarydf[column+replicatename].items():
                # if value in controls is larger than 20%, set to None (ignore for following analysis)
                    if controldf.at[index,column+'_average_raw'] > 20:
                        summarydf.at[index,column+replicatename+'_ctrladjusted'] = None
                    else:
                        summarydf.at[index,column+replicatename+'_ctrladjusted'] = (value-controldf.at[index,column+'_average_raw'])/(100-controldf.at[index,column+'_average_raw'])*100
        
        except KeyError:  # skip this if for this editor there is no CtoT 
            break
    for column in AtoGcolumnlist:
        try:
            summarydf[column+'_average_raw'] = summarydf[[column+'_rep1',column+'_rep2',column+'_rep3']].mean(axis=1)
            for index, value in summarydf[column+'_average_raw'].items():
                if controldf.at[index,column+'_average_raw'] > 20:
                    summarydf.at[index,column+'_average_ctrladjusted'] = None
                else:
                    summarydf.at[index,column+'_average_ctrladjusted'] = (value-controldf.at[index,column+'_average_raw'])/(100-controldf.at[index,column+'_average_raw'])*100
            for replicatename in ['_rep1','_rep2','_rep3']: 
                for index, value in summarydf[column+replicatename].items():
                # if value in controls is larger than 20%, set to None (ignore for following analysis)
                    if controldf.at[index,column+'_average_raw'] > 20:
                        summarydf.at[index,column+replicatename+'_ctrladjusted'] = None
                    else:
                        summarydf.at[index,column+replicatename+'_ctrladjusted'] = (value-controldf.at[index,column+'_average_raw'])/(100-controldf.at[index,column+'_average_raw'])*100
        
            
        except KeyError:  # skip this if for this editor there is no AtoG
            break
    summarydf['percent_unedited'+'_average_raw'] = summarydf[['percent_unedited'+'_rep1','percent_unedited'+'_rep2','percent_unedited'+'_rep3']].mean(axis=1)
    summarydf['percent_modified'+'_average_raw'] = summarydf['percent_unedited'+'_average_raw'].apply(lambda x: 100-x)

    for index, value in summarydf['percent_modified'+'_average_raw'].items():
        if controldf.at[index,'percent_modified'+'_average_raw'] > 20:
            summarydf.at[index,column+'_average_ctrladjusted'] = None
        else:
            summarydf.at[index,'percent_modified'+'_average_ctrladjusted'] = (value-controldf.at[index,'percent_modified'+'_average_raw'])/(100-controldf.at[index,'percent_modified'+'_average_raw'])*100
    
    for index, row in summarydf.iterrows():
        highcoveragelist = []
        for readnr in ['nrreads_rep1','nrreads_rep2','nrreads_rep3']:
            if row[readnr] > 100:
                highcoveragelist.append(1)
            else:
                highcoveragelist.append(0)
        
        summarydf.at[index,'highcoveragereps'] = sum(highcoveragelist)
        summarydf.at[index,'controlhighcoveragereps'] = controldf.at[index,'highcoveragereps']
    
    collist = list(summarydf.columns)
    
    rawlist = [x for x in collist if 'raw' in x]
    adjustlist = [x for x in collist if 'adjust' in x]
    
    collist = [x for x in collist if (x not in rawlist)]
    collist = [x for x in collist if (x not in adjustlist)]
    
    filtercols = collist + rawlist + adjustlist
    
    summarydf = summarydf[filtercols]
    summarydf.to_csv('AnalysisFiles/AverageFiles/'+'20230319_'+sample+'_averagedf.csv', index=None)
    # shortdfcols = [x for x in filtercols if 'raw' not in x]
    # shortdfcols = [x for x in shortdfcols if 'rep' not in x]
    # shortdf = summarydf[shortdfcols]
    # shortdf.to_csv('AverageFiles/'+'20220530_'+sample+'_averagedf_short.csv', index=None)



