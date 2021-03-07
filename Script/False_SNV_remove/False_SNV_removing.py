#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 14:18:08 2020

@author: Chong Li

"""
import sys
#This is the file to show the nucleotide differences between challenge strains by segment
inputfile_1= sys.argv[1]
Ref_file=open(inputfile_1,'r',encoding="ISO-8859-1")
Ref_file.readline()
Ref_pos=[]
Ref_snp=[]
for line_1 in Ref_file:
    line_1=line_1.strip()
    split_file_1=line_1.split(sep=',')
    Ref_pos.append(split_file_1[1])
    snp=split_file_1[0]
    Ref_snp.append(snp[0]+snp[2])

Ref_file.close()
# The input flle which need to be corrected.
inputfile_2= sys.argv[2]
SNP_file=open(inputfile_2,'r')

Chrom=[]
Position=[]
Ref=[]
Cons=[]
Reads1=[]
Reads2=[]
VarFreq=[]
Strands1=[]
Strands2=[]
Qual1=[]
Qual2=[]
Pvalue=[]
MapQual1=[]
MapQual2=[]
Reads1Plus=[]
Reads1Minus=[]
Reads2Plus=[]
Reads2Minus=[]
VarAllele=[]
SNP_C=[]
SNP_R=[]
title=SNP_file.readline()

for line_2 in SNP_file:
    line_2=line_2.strip()
    split_file_2=line_2.split(sep='\t')
    Chrom.append(split_file_2[0])
    Position.append(split_file_2[1])
    Ref.append(split_file_2[2])
    Cons.append(split_file_2[3])
    Reads1.append(split_file_2[4])
    Reads2.append(split_file_2[5])
    VarFreq.append(split_file_2[6])
    Strands1.append(split_file_2[7])
    Strands2.append(split_file_2[8])
    Qual1.append(split_file_2[9])
    Qual2.append(split_file_2[10])
    Pvalue.append(split_file_2[11])
    MapQual1.append(split_file_2[12])
    MapQual2.append(split_file_2[13])
    Reads1Plus.append(split_file_2[14])
    Reads1Minus.append(split_file_2[15])
    Reads2Plus.append(split_file_2[16])
    Reads2Minus.append(split_file_2[17])
    VarAllele.append(split_file_2[18])
    SNP_C.append(split_file_2[2]+split_file_2[18])
    SNP_R.append(split_file_2[18]+split_file_2[2])
SNP_file.close()
index=[]   
for index_1, words_1 in enumerate (Position):
    if not words_1 in Ref_pos:
       index.append(int(index_1))
    else: 
       Index_2= Ref_pos.index(str(words_1))
       Ref_snp_2=Ref_snp[Index_2]
       if str(Ref_snp_2)!=str(SNP_C[int(index_1)]) and  str(Ref_snp_2)!=str(SNP_R[int(index_1)]):
           index.append(int(index_1))
           
outfile_1= sys.argv[3]
out_file_1=open(outfile_1,'w')
out_file_1.write(title)

for i in range(0,len(index),1):
    number=int(index[i])
    print(str(Chrom[number]),str(Position[number]),str(Ref[number]),str(Cons[number]),str(Reads1[number]),str(Reads2[number]),str(VarFreq[number]),str(Strands1[number]),str(Strands2[number]),str(Qual1[number]),str(Qual2[number]),str(Pvalue[number]),str(MapQual1[number]),str(MapQual2[number]),str(Reads1Plus[number]),str(Reads1Minus[number]),str(Reads2Plus[number]),str(Reads2Minus[number]),str(VarAllele[number]),sep='\t', end='\n', file= out_file_1)
out_file_1.close()
