#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 07:38:58 2021

@author: bioinfo
"""

import numpy as np
import pandas as pd
import glob
import os
import csv

##define path
path1 = "/home/bioinfo/Documents/Valeria/NAFuse/files_bed/"
path2 = "/home/bioinfo/Documents/Valeria/NAFuse/files_bed/"
out_path = "/home/bioinfo/Documents/Valeria/NAFuse/output/"


###Extract file names
all_files = glob.glob(os.path.join(path1, "*_read1.bed"))

for file in all_files:
    name = os.path.splitext(os.path.basename(file))[0]
    sample_name = name.replace('_NAFuse_read1', '')
    filenames1 = glob.glob(os.path.join(path1+sample_name+ "*NAFuse_read1.bed"))
    filenames2 = glob.glob(os.path.join(path2+sample_name+ "*NAFuse_read2.bed"))
    print("Doing sample", sample_name)
    ###Fusion
    fusion_path = "/home/bioinfo/Documents/Valeria/NAFuse/Fusions/"
    fusion_output = pd.read_csv(fusion_path+sample_name+".txt", header=None, delimiter= " ", low_memory=False, names=("Gene1","Gene2"))
    gene1 = fusion_output['Gene1'].tolist()
    fusions1= [f'/b{x}/b' for x in gene1]
    gene2 = fusion_output['Gene2'].tolist()
    fusions2= [f'/b{x}/b' for x in gene2]
    print("Upload fusions")
    matches1= []
    for file1 in filenames1:
        outfile1 = open(file1,'r')
        data1 = outfile1.readlines()
        outfile1.close()
        for item1 in data1:
            lin = item1.replace('\t','/b')
            line1 = lin.replace('\n','/b')
            for f in fusions1:
                if f in line1:
                    match_line1 = line1
                    match1 = match_line1.replace('/b',',')
                    matches1.append(match1)
    final1 = pd.DataFrame(matches1)
    final1 = final1[0].str.split(',',expand=True)
    print("Read1 Done")
    matches2= []
    for file2 in filenames2:
        outfile2 = open(file2,'r')
        data2 = outfile2.readlines()
        outfile2.close()
        for item2 in data2:
            lin = item2.replace('\t','/b')
            line2 = lin.replace('\n','/b')
            for d in fusions2:
                if d in line2:
                    match_line2 = line2
                    match2 = match_line2.replace('/b',',')
                    matches2.append(match2)
    final2 = pd.DataFrame(matches2)
    final2 = final2[0].str.split(',',expand=True)
    print("Read2 Done")
    sample_merged= pd.merge(final1,final2,how='inner',on=0,suffixes=('_1', '_2'))
    l = sample_merged.values.tolist()
    print("Starting NAFuse core")
    NAFuse_out= []
    for x in l:
        n='2'
        ind = x.index(n)
        list1_result = x[0:ind]
        list2_result = x[ind:-1]
        for f in range(len(fusion_output)):
            gene1= fusion_output["Gene1"][f]
            gene2= fusion_output["Gene2"][f]
            if gene1 in list1_result:
                if gene2 in list2_result:
                    NAFuse_out.append([gene1, gene2])
    out = pd.DataFrame(NAFuse_out).drop_duplicates()
    out.to_csv(out_path+sample_name+'.txt', index = False,  sep= "\t", header= ['5_end', '3_end'])
    print("Done")


