#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 15:21:37 2020

@author: akhorgad
"""

import os
import argparse
import random
import numpy as np
import pandas as pd

#process 1.1
#reads the 012 file into a int8 dataframew
def o12toDataframe(o12):
    o12_df=pd.read_csv(o12, sep='\t', header=None, dtype=np.int8)
    o12_df=o12_df.drop(columns=[0])
    o12_df.columns=range(o12_df.shape[1])
    #nrow=len(o12_df)
    #ncol=len(o12_df.columns)
    return(o12_df)

#process 1.2        
#computes minor allele frequencies - takes the 012 and the pos file from the arguments and returns a df of positions and mafs

def computeMaf(o12_df,pos_f):
    colnames=['Pos','refc','altc','missingc','hetc','ref_freq','alt_freq','maj_allele']
    dtypes={'Pos':'str', 'refc':np.int32,'altc':np.int32,'missingc':np.int32,'hetc':np.int32,'ref_freq':np.float32,'alt_freq':np.float32,'maj_allele':np.int8}
    freq_df=pd.read_csv(pos_f,names=colnames, header=None)
    freq_df=freq_df.apply(lambda y : np.nan_to_num(y))
    freq_df=freq_df.astype(dtypes)
    #print(freq_df.dtypes)
    for col in o12_df.columns:
        colarr=np.array(o12_df[col])
        allele, counts=np.unique(colarr, return_counts=True)
        allele_dict=dict(zip(allele, counts)) 
        #print (allele_dict)
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
    maf05_df.to_csv('./out/maf05_freqs.tsv', sep='\t',index=False)
    #randomly drawing out of 0 and 1 for the positions with equal frequencies for ref and alt that is maf of 0.5 
    freq_df['maj_allele']=freq_df['maj_allele'].replace([3], random.choice([0, 2]))
    #print(freq_df.dtypes)
    return(freq_df)

#process 2.1
#trimming variant positons from arg.window size and arg 012.pos
def trimPositions(pos, window_size):
    colnames=['chromosome']
    pos_df=pd.read_csv(pos,names=colnames, header=None)
    pos_df['chr']=pos_df.chromosome.apply(lambda x: x.split(":")[0])
    pos_df['pos']=pos_df.chromosome.apply(lambda x: x.split(":")[1])
    pos_df=pos_df.drop(columns=['chromosome'])
    filtered_chr_pos_df=pos_df.groupby('chr').agg(lambda chr_series: pruningWindow(chr_series,window_size))
    filtered_chr_pos_df['npos']=filtered_chr_pos_df['pos'].apply(len)
    #filtered_chr_pos_df.to_csv('./out/selected_positions.tsv', sep='\t',index=True)
    print("Total number of positions selected: "+str(filtered_chr_pos_df['npos'].sum()))
    return(filtered_chr_pos_df, pos_df) 

#process 2.1.2
#group by  will take array of positions prune it by window size 
def pruningWindow(pos_arr_per_chr,window_size):
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
def addPositionsAndSampleIdsAndTranspose(o12_df, pos_f, indv_f):
    colnames=['chromosome']
    sample_ids=['indv']
    pos_df=pd.read_csv(pos_f,names=colnames, header=None)
    indv_df=pd.read_csv(indv_f,names=sample_ids, header=None)
    o12_df.columns=list(pos_df['chromosome'])
    o12_df.index = list(indv_df['indv'])
    o12_transposed=o12_df.transpose()
    return(o12_transposed)


#process 3.2
#prunes the o12 dataframe to keep the selected SNP's 
def selectSNP(filtered_chr_pos_df,o12_transposed, freq_df):
    selectSNP_list=[]
    for index, row in filtered_chr_pos_df.iterrows():
        for i in range(0,row['npos']):
            selectSNP_list.append(index+':'+str(row['pos'][i]))
    selectSNP_df=pd.DataFrame(selectSNP_list)
    o12_selectSNP_df=pd.merge(selectSNP_df,o12_transposed,left_on=0,right_index=True,how='left')
    freq_selected_df=pd.merge(selectSNP_df,freq_df,left_on=0,right_on='Pos',how='left')
    o12_selectSNP_df.index = list(o12_selectSNP_df[0])
    o12_selectSNP_df=o12_selectSNP_df.drop(columns=[0])
    o12_selectSNP_df=o12_selectSNP_df.transpose()
    freq_selected_df=freq_selected_df.drop(columns=[0])
    return(o12_selectSNP_df,freq_selected_df)


#process 3.3
#takes in the 012 file of selected snp's and converts ref-alt to minor-major allele frequencies
#----------------------------------
# vcftools 012 format:
# 0 Homozygous Reference 
# 1 Heterozygous   (coded to 0.5)
# 2 Homozygous Alternate
# -1 Missing
#----------------------------------    


def recodeDf(o12_selectSNP_df, freq_selected_df):
    colarr=np.array(freq_selected_df['maj_allele'], dtype=np.int8)
    final_recoded_df=pd.DataFrame(index=o12_selectSNP_df.index)
    for i in range(0,len(colarr)):
        final_recoded_df[i]=o12_selectSNP_df.apply(lambda x: 1 if x[i]==colarr[i] else x[i] if x[i]==-1 else 0.5 if x[i]==1 else 0, axis=1)
    final_recoded_df.columns=o12_selectSNP_df.columns
    final_recoded_df.insert(0,'sample',final_recoded_df.index)
    return(final_recoded_df)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o12",required=True, help="path to 012 file.")
    parser.add_argument("-indv",required=True, help="path to 012.indv file")
    parser.add_argument("-pos", required=True, help="path to 012.pos file")
    parser.add_argument("-window_size", required=True, type=int, help="window size to prune variants")
    args = parser.parse_args()
    
    if os.path.isfile(args.o12 and args.pos and args.indv):
        opath=os.path.join( os.getcwd(),'out')
        os.mkdir(opath)
        o12_df=o12toDataframe(args.o12)
        allele_freq_df=computeMaf(o12_df, args.pos)
        allele_freq_df.to_csv('./out/allele_freqs.tsv', sep='\t')
        filtered_chr_pos_df, pos_df =trimPositions(args.pos, args.window_size)
        filtered_chr_pos_df.to_csv('./out/position_list.tsv', sep='\t')
        o12_transposed = addPositionsAndSampleIdsAndTranspose(o12_df, args.pos, args.indv)
        o12_selectedSNPs_df, freq_selectedSNPs_df = selectSNP(filtered_chr_pos_df, o12_transposed, allele_freq_df)
        final_recoded_df=recodeDf(o12_selectedSNPs_df, freq_selectedSNPs_df)
        final_recoded_df.to_csv('./out/final_recoded_matrix.tsv', sep='\t')
    else:
        print("ERROR: Invalid file path provided for one of the arguments.\n Check the filepaths privided" )
        print("012 filepath: "+args.o12)
        print("pos filepath: "+args.pos)
        print("indv filepath: "+args.pos)
  
if __name__ == "__main__":
    main()



     
