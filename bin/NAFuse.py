#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 10:25:44 2021
@author: Valeria Difilippo
"""
import numpy as np
import pandas as pd
import glob
import os
import csv
import sys
import argparse


parser = argparse.ArgumentParser(prog='NAFuse.py',description='The tool allows for the detection of breakpoints in both partner genes/regions on both the RNA and DNA levels without having to manually search output files generated during downstream processing of BAM-files.')
parser.add_argument('Read1', nargs='+', help='Read1 of a case')
parser.add_argument('Read2', nargs='+', help='Read2 of a case')
parser.add_argument('Fusion', nargs='+', help='Fusion list')
parser.add_argument('Output', nargs='+', help='output')

args = parser.parse_args()
print(args.accumulate(args.integers))

filenames1 = glob.glob(os.path.join(sys.argv[1]+ "*NAFuse_read1.bed"))
filenames2 = glob.glob(os.path.join(sys.argv[2]+ "*NAFuse_read2.bed"))


###Fusion
fusion_output = pd.read_csv(sys.argv[3]+".txt", header=None, delimiter= " ", low_memory=False, names=("Gene1","Gene2"))

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
out.to_csv(sys.argv[4]+'.txt', index = False,  sep= "\t", header= ['5_end', '3_end'])

print("Done")
