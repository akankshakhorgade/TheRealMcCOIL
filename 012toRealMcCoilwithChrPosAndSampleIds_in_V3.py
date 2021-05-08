#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 15:21:37 2020

@author: Akanksha Khorgade
@contact: akhorgad@broadinstitute.org
"""

#Creating a RealMcCoil input format from 012 file with positional pruning 
#or selecting on chromosomal positions 
#or selecting on individuals

######### Changes made from vcftools012toRealMcCoil_in.py 
#------takes a list of chromosomal positions to filter the o/p on
#------additionally, takes a list of sample ids to filter th o/p on
#------return a file with the rest of the logic as it is but slected on postions and individuals


import os
import argparse
import random
import numpy as np
import pandas as pd

#process 1.1
#reads the 012 file into a int8 dataframe
def o12toDataframe(o12_f):
    o12_df=pd.read_csv(o12_f, sep='\t', header=None, dtype=np.int8)
    o12_df=o12_df.drop(columns=[0])
    o12_df.columns=range(o12_df.shape[1])
    return(o12_df)

#process 1.2        
#computes minor allele frequencies - takes the 012 and the pos file from the arguments and returns a df of positions and mafs
def computeMaf(o12_df, pos_f):
    colnames=['Pos','refc','altc','missingc','hetc','ref_freq','alt_freq','maj_allele']
    dtypes={'Pos':'str', 'refc':np.int32,'altc':np.int32,'missingc':np.int32,'hetc':np.int32,'ref_freq':np.float32,'alt_freq':np.float32,'maj_allele':np.int8}
    print("Reading pos file..")
    freq_df=pd.read_csv(pos_f,names=colnames, header=None)
    freq_df=freq_df.apply(lambda y : np.nan_to_num(y))
    print(freq_df.head())
    freq_df=freq_df.astype(dtypes)
    for col in o12_df.columns:
        colarr=np.array(o12_df[col])
        allele, counts=np.unique(colarr, return_counts=True)
        allele_dict=dict(zip(allele, counts)) 
        print (allele_dict)
        #dict with count of alleles ex: {-1: 6, 0: 3, 2: 1}
        freq_df.loc[col,'refc']=np.int32(allele_dict.get(0,0))
        freq_df.loc[col,'altc']=np.int32(allele_dict.get(2,0))
        freq_df.loc[col,'missingc']=np.int32(allele_dict.get(-1,0))
        freq_df.loc[col,'hetc']=np.int32(allele_dict.get(1,0))
        homo_allele_count=freq_df.loc[col,'refc']+freq_df.loc[col,'altc']
        freq_df.loc[col,'ref_freq']=0 if (homo_allele_count==0) else np.float32(freq_df.loc[col,'refc'])/np.float32(homo_allele_count)
        freq_df.loc[col,'alt_freq']=0 if (homo_allele_count==0) else np.float32(freq_df.loc[col,'altc'])/np.float32(homo_allele_count)
        freq_df.loc[col,'maj_allele']=np.int8(0) if freq_df.loc[col,'refc']>freq_df.loc[col,'altc'] else np.int8(2) if freq_df.loc[col,'refc']<freq_df.loc[col,'altc'] else np.int8(3) if (freq_df.loc[col,'refc']==freq_df.loc[col,'altc'] and freq_df.loc[col,'refc']!=0) else np.int8(4);
    maf05_df=freq_df.loc[freq_df['maj_allele']==3]
    maf05_df.to_csv('./maf05_freqs.csv', sep='\t',index=False)
    #randomly drawing out of 0 and 1 for the positions with equal frequencies for ref and alt that is maf of 0.5 
    freq_df['maj_allele']=freq_df['maj_allele'].replace([3], random.choice([0, 2]))
    #print(freq_df.dtypes)
    return(freq_df)

##process 2.1.1
##trimming variant positons from arg.window size and arg 012.pos or optionally args.selectedSNP
def trimPositions(select_pos_f, window_size): 
    selectSNP_list=[]
    with open(select_pos_f) as f:
        selectSNP_list = [line.rstrip() for line in f]
    pos_df=pd.DataFrame(selectSNP_list, columns=['chrPos'])
    print(pos_df.head())
    pos_df['chr']=pos_df.chrPos.apply(lambda x: x.split(":")[0])
    pos_df['pos']=pos_df.chrPos.apply(lambda x: x.split(":")[1])
    pos_df=pos_df.drop(columns=['chrPos'])
    print(pos_df.head())
    filtered_chr_pos_df=pos_df.groupby('chr').agg(lambda chr_series: pruningWindow(chr_series,window_size))
    filtered_chr_pos_df['npos']=filtered_chr_pos_df['pos'].apply(len)
    #filtered_chr_pos_df.to_csv('./selected_positions.csv', sep='\t',index=True)
    print("Total number of positions selected: "+str(filtered_chr_pos_df['npos'].sum()))
    return(filtered_chr_pos_df, pos_df) 


#process 2.1.2  
#group by  will take array of positions prune it by window size 
def pruningWindow(pos_arr_per_chr, window_size):
    pos_arr=np.array(pos_arr_per_chr, dtype=np.int32)
    pos_min=pos_arr.min()
    pos_max=pos_arr.max() 
    pos_list=[]
    for i in range(pos_min, pos_max+1, window_size):
        random_pos_arr = list(filter(lambda x: i <= x < i+window_size, pos_arr))
        if len(random_pos_arr)!=0:
            random_pos = np.random.choice(random_pos_arr,1, replace=False)
            pos_list.append(random_pos[0]) 
    return(pos_list) 
 

#process 3.1
#adds the sample ids and the positions to the 012 file
def addPositionsAndSampleIds(o12_df, pos_f, indv_f):
    colnames=['chromosome']
    sample_ids=['indv']
    pos_all_df=pd.read_csv(pos_f,names=colnames, header=None)
    indv_df=pd.read_csv(indv_f,names=sample_ids, header=None)
    o12_df.columns=list(pos_all_df['chromosome'])
    o12_df.index=list(indv_df['indv'])
    return(o12_df)

#process 3.2- if window_size in args
#prunes the o12 dataframe for LD to keep the selected SNP's 
def selectSNP(filtered_chr_pos_df, o12_df, freq_df):
    o12_transposed=o12_df.transpose()
    selectSNP_list=[]
    for index, row in filtered_chr_pos_df.iterrows():
        for i in range(0,row['npos']):
            selectSNP_list.append(index+':'+str(row['pos'][i]))
    selectSNP_df=pd.DataFrame(selectSNP_list)
    o12_selectSNP_df=pd.merge(selectSNP_df,o12_transposed,left_on=0,right_index=True,how='left')
    o12_selectSNP_df=o12_selectSNP_df.dropna()
    freq_selected_df=pd.merge(selectSNP_df,freq_df,left_on=0,right_on='Pos',how='left')
    freq_selected_df=freq_selected_df.dropna()
    o12_selectSNP_df.index=list(o12_selectSNP_df[0])
    o12_selectSNP_df=o12_selectSNP_df.drop(columns=[0])
    o12_selectSNP_df=o12_selectSNP_df.transpose()
    freq_selected_df=freq_selected_df.drop(columns=[0])
    return(o12_selectSNP_df,freq_selected_df)

#process 3.2- if selectSNP in args
#optionally - prunes the o12 dataframe to keep the selected SNP's 
def selectSNPfromFile(select_pos_f, o12_df, freq_df):
    o12_transposed=o12_df.transpose()
    with open(select_pos_f) as f:
        selectSNP_list = [line.rstrip() for line in f]
    selectSNP_df=pd.DataFrame(selectSNP_list, columns=['Pos'])
    print(selectSNP_df.head())
    o12_selectSNP_df=pd.merge(selectSNP_df,o12_transposed, left_on='Pos',right_index=True,how='left')
    o12_selectSNP_df=o12_selectSNP_df.dropna()
    freq_selected_df=pd.merge(selectSNP_df,freq_df, on='Pos',how='left')
    freq_selected_df=freq_selected_df.dropna()
    print(freq_selected_df.head())
    o12_selectSNP_df.index = list(o12_selectSNP_df['Pos'])
    o12_selectSNP_df=o12_selectSNP_df.drop(columns=['Pos'])
    o12_selectSNP_df=o12_selectSNP_df.transpose()
    print(o12_selectSNP_df.head())
    return(o12_selectSNP_df, freq_selected_df)


#process 3.2- if selectIndv in args
#optionally - prunes the o12 dataframe to keep the selected individuals 
def selectIndv(select_indv_f, o12_selected_df): 
    with open(select_indv_f) as f:
        selectIndv_list = [line.rstrip() for line in f]
    selectIndv_df=pd.DataFrame(selectIndv_list, columns=['Indv'])
    print(selectIndv_df.head())
    o12_selectIndv_df=pd.merge(selectIndv_df,o12_selected_df,left_on='Indv',right_index=True,how='left')
    o12_selectIndv_df=o12_selectIndv_df.dropna()
    o12_selectIndv_df.index = list(o12_selectIndv_df['Indv'])
    o12_selectIndv_df=o12_selectIndv_df.drop(columns=['Indv'])
    print(o12_selectIndv_df.head())
    return(o12_selectIndv_df)

#process 3.3
#takes in the 012 file of selected snp's and converts ref-alt to minor-major allele frequencies
#----------------------------------
# vcftools 012 format:
# 0 Homozygous Reference 
# 1 Heterozygous   (coded to 0.5)
# 2 Homozygous Alternate
# -1 Missing
#----------------------------------    

def recodeDf(o12_selected_df, freq_selected_df):
    colarr=np.array(freq_selected_df['maj_allele'], dtype=np.int8)
    final_recoded_df=pd.DataFrame(index=o12_selected_df.index)
    for i in range(0,len(colarr)):
        final_recoded_df[i]=o12_selected_df.apply(lambda x: 1 if x[i]==colarr[i] else x[i] if x[i]==-1 else 0.5 if x[i]==1 else 0, axis=1)
    final_recoded_df.columns=o12_selected_df.columns
    final_recoded_df.insert(0,'sample',final_recoded_df.index)
    return(final_recoded_df)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o12",required=True, help="path to 012 file")
    parser.add_argument("-indv",required=True, help="path to 012.indv file")
    parser.add_argument("-pos", required=True, help="path to 012.pos file")
    parser.add_argument("-out", required=True, help="path to create out folder")
    parser.add_argument("-window_size", required=False, type=int, help="window size to prune variants")
    parser.add_argument('-allSNPs', action='store_true', help="set to have all SNPs i.e. window size=1 to select variants")
    parser.add_argument("-selectSnp", required=False, help="path to chr:pos file for selected positions")
    parser.add_argument("-selectIndv", required=False, help="path to file for selected individuals/sample ids")
    args = parser.parse_args()
    
    print("Initialising..")
    if os.path.isfile(args.o12 and args.pos and args.indv):
        print("Input files..")
        opath=os.path.join(os.getcwd(),args.out)
        os.mkdir(opath)
        os.chdir(opath)
        o12_df=o12toDataframe(args.o12)
        print(o12_df.head())
        allele_freq_df=computeMaf(o12_df, args.pos)
        print(allele_freq_df.head())
        allele_freq_df.to_csv("./allele_freqs.csv", sep='\t')
        o12_df=addPositionsAndSampleIds(o12_df, args.pos, args.indv)
        print(o12_df.head())
        if args.selectSnp and os.path.isfile(args.selectSnp):
            print("Selecting positions from selectSNP file..")
            if args.window_size:
                filtered_chr_pos_df, pos_df=trimPositions(args.selectSnp, args.window_size)
                pos_df.to_csv("./position_list.csv", sep='\t')
                #filtered_chr_pos_df('./selected_pos_per_chr.csv', sep='\t')
                o12_selectedSNPs_df, freq_selectedSNPs_df=selectSNP(filtered_chr_pos_df, o12_df, allele_freq_df)
                freq_selectedSNPs_df.to_csv("./freq_selectedSNPs.csv", sep='\t')
                final_recoded_df=recodeDf(o12_selectedSNPs_df, freq_selectedSNPs_df)
            elif args.allSNPs:
                print("Selecting all SNP positions..")
                o12_selectedSNPs_df, freq_selectedSNPs_df=selectSNPfromFile(args.selectSnp, o12_df, allele_freq_df)
                print(o12_selectedSNPs_df.head())
                print(freq_selectedSNPs_df.head())
                freq_selectedSNPs_df.to_csv("./freq_selectedSNPs.csv", sep='\t')
                final_recoded_df=recodeDf(o12_selectedSNPs_df, freq_selectedSNPs_df)
                print(final_recoded_df.head())
        elif args.selectIndv and os.path.isfile(args.selectIndv):
            if args.window_size:
                filtered_chr_pos_df, pos_df=trimPositions(args.pos, args.window_size)
                pos_df.to_csv("./position_list.csv", sep='\t')
                #filtered_chr_pos_df('./selected_pos_per_chr.csv', sep='\t')
                o12_selectedSNPs_df, freq_selectedSNPs_df=selectSNP(filtered_chr_pos_df, o12_df, allele_freq_df)
                freq_selectedSNPs_df.to_csv("./freq_selectedSNPs.csv", sep='\t')
                o12_selectedIndv_df=selectIndv(args.selectIndv, o12_selectedSNPs_df)
                final_recoded_df=recodeDf(o12_selectedIndv_df, freq_selectedSNPs_df)
            elif args.allSNPs:
                o12_selectedIndv_df=selectIndv(args.selectIndv, o12_df)
                final_recoded_df=recodeDf(o12_selectedIndv_df, allele_freq_df)            
        else:
            print("ERROR: Invalid file path provided for one of the arguments.\n Check the filepaths provided")
        
        final_recoded_df.to_csv("./final_recoded_matrix.csv", sep='\t')
    else:
        print("ERROR: Invalid file path provided for one of the arguments.\n Check the filepaths provided")
        print("012 filepath: "+args.o12)
        print("pos filepath: "+args.pos)
        print("indv filepath: "+args.pos)
  
if __name__ == "__main__":
    main()



     
