#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 23:14:46 2020

@author: apple
"""
import random
import collections
import pandas as pd
import seaborn as sns
from random import sample
from collections import Counter

    
#Generate the snp file for H1 virus.    

with open("SNP_H1_protein.csv", "r") as SNP_H1_protein:
    

    Sample_H1=[]
    Group_H1=[]
    PB2_H1=[]
    PB1_H1=[]
    PB1_F2_H1=[]
    PA_H1=[]
    PA_X_H1=[]
    NP_H1=[]
    HA_H1=[]
    NA_H1=[]
    M1_H1=[]
    M2_H1=[]
    NS1_H1=[]
    NS2_H1=[]
    title_H1=SNP_H1_protein.readline()

    for line_H1 in SNP_H1_protein:
        line_H1=line_H1.strip()
        split_file_H1=line_H1.split(sep=',')
        Sample_H1.append(split_file_H1[0])
        Group_H1.append(split_file_H1[2])
        PB2_H1.append(split_file_H1[5])
        PB1_H1.append(split_file_H1[6])
        PB1_F2_H1.append(split_file_H1[7])
        PA_H1.append(split_file_H1[8])
        PA_X_H1.append(split_file_H1[9])
        HA_H1.append(split_file_H1[10])
        NP_H1.append(split_file_H1[11])
        NA_H1.append(split_file_H1[12])
        M1_H1.append(split_file_H1[13])
        M2_H1.append(split_file_H1[14])
        NS1_H1.append(split_file_H1[15])
        NS2_H1.append(split_file_H1[16])

COM_COM_H1={}
AUT_AUT_H1={}
LAIV_NONE_H1={}
LAIV_COM_H1={}
NO_VAC_H1={}

for index_H1, words_H1 in enumerate (Group_H1):
    if words_H1 == "1":
       COM_COM_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]} 
    if words_H1 == "2":
       AUT_AUT_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]}
    if words_H1 == "5":
       LAIV_NONE_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]}
    if words_H1 == "6":
       LAIV_COM_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]} 
    if words_H1 == "7":
       NO_VAC_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]}  
    
# Generate the snp files for H3 virus.
with open("SNP_H3_protein.csv", "r") as SNP_H3_protein:
    

    Sample_H3=[]
    Group_H3=[]
    PB2_H3=[]
    PB1_H3=[]
    PB1_F2_H3=[]
    PA_H3=[]
    PA_X_H3=[]
    NP_H3=[]
    HA_H3=[]
    NA_H3=[]
    M1_H3=[]
    M2_H3=[]
    NS1_H3=[]
    NS2_H3=[]
    title_H3=SNP_H3_protein.readline()

    for line_H3 in SNP_H3_protein:
        line_H3=line_H3.strip()
        split_file_H3=line_H3.split(sep=',')
        Sample_H3.append(split_file_H3[0])
        Group_H3.append(split_file_H3[2])
        PB2_H3.append(split_file_H3[5])
        PB1_H3.append(split_file_H3[6])
        PB1_F2_H3.append(split_file_H3[7])
        PA_H3.append(split_file_H3[8])
        PA_X_H3.append(split_file_H3[9])
        HA_H3.append(split_file_H3[10])
        NP_H3.append(split_file_H3[11])
        NA_H3.append(split_file_H3[12])
        M1_H3.append(split_file_H3[13])
        M2_H3.append(split_file_H3[14])
        NS1_H3.append(split_file_H3[15])
        NS2_H3.append(split_file_H3[16])

COM_COM_H3={}
AUT_AUT_H3={}
LAIV_NONE_H3={}
LAIV_COM_H3={}
NO_VAC_H3={}
COM_AUT_H3={}

for index_H3, words_H3 in enumerate (Group_H3):
    if words_H3 == "1":
       COM_COM_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]} 
    if words_H3 == "2":
       AUT_AUT_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]} 
    if words_H3 == "5":
       LAIV_NONE_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]} 
    if words_H3 == "6":
       LAIV_COM_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]} 
    if words_H3 == "7":
       NO_VAC_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]} 
    if words_H3 == "4":
       COM_AUT_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]}  
    
#Generate the set of the numbers of nucleotides in each gene segment.
Nucleotide_region_H1={}
Nucleotide_region_H3={}

for num_H1 in range(0,len(Sample_H1),1):
    Nucleotide_region_H1[Sample_H1[int(num_H1)]]={"PB2":759, "PB1":757, "PB1-F2":79, "PA":716, "PA-X":232, "HA":566, "NP":498, "NA":469, "M1":252, "M2":97, "NS1":219, "NS2":121}
for num_H3 in range(0,len(Sample_H3),1):
    Nucleotide_region_H3[Sample_H3[int(num_H3)]]={"PB2":759, "PB1":757, "PB1-F2":79, "PA":716, "PA-X":232, "HA":566, "NP":498, "NA":469, "M1":252, "M2":97, "NS1":219, "NS2":121}

#Generate the set of the 70% numbers of nucleotides in each gene segment.
Nucleotide_region_70_H1={}
Nucleotide_region_70_H3={}

for num_H1_70 in range(0,len(Sample_H1),1):
    Nucleotide_region_70_H1[Sample_H1[int(num_H1_70)]]={"PB2":531, "PB1":530, "PB1-F2":55, "PA":501, "PA-X":162, "HA":396, "NP":349, "NA":328, "M1":176, "M2":68, "NS1":153, "NS2":85}
for num_H3_70 in range(0,len(Sample_H3),1):
    Nucleotide_region_70_H3[Sample_H3[int(num_H3_70)]]={"PB2":531, "PB1":530, "PB1-F2":55, "PA":501, "PA-X":162, "HA":396, "NP":349, "NA":328, "M1":176, "M2":68, "NS1":153, "NS2":85}

#build the function
def permutation_test(SNP_array, Nucleotide_region, stimulation):
    
    total_shared_sites=[]
    iteration_results = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
    for a in range(0, stimulation, 1):
        temporary_results = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[],"NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
        for samples in SNP_array:
            for gene in SNP_array[samples]:
                Nucleotides_range = range(0,Nucleotide_region[samples][gene])
                
                
                Num_draw_list = SNP_array[samples][gene]
                Num_draw = int(Num_draw_list[0])
                
            
                random_draw = (random.sample(Nucleotides_range, Num_draw))
                    
                for r in random_draw:
                    temporary_results[gene].append(r)
                    
        for gene in temporary_results:
            count=collections.Counter(temporary_results[gene]).values()
            
            more_than_once = 0
            for i in count:
                if i > 1:
                   more_than_once += 1
            
            iteration_results[gene].append(more_than_once)

    for c in range (0, stimulation, 1):
        shared_sites = iteration_results["PB2"][c]+iteration_results["PB1"][c]+ iteration_results["PB1-F2"][c]+iteration_results["PA"][c]+iteration_results["PA-X"][c]+iteration_results["HA"][c]+iteration_results["NP"][c]+iteration_results["NA"][c]+iteration_results["M1"][c]+iteration_results["M2"][c]+iteration_results["NS1"][c]+iteration_results["NS2"][c]
        total_shared_sites.append(shared_sites)

    return(total_shared_sites) 

#run the method.
stimulation = 100000
H1_COM_COM= permutation_test(COM_COM_H1, Nucleotide_region_H1, stimulation)    
H1_AUT_AUT= permutation_test(AUT_AUT_H1, Nucleotide_region_H1, stimulation)
H1_LAIV_NONE= permutation_test(LAIV_NONE_H1, Nucleotide_region_H1, stimulation)    
H1_LAIV_COM= permutation_test(LAIV_COM_H1, Nucleotide_region_H1, stimulation)
H1_NO_VAC= permutation_test(NO_VAC_H1, Nucleotide_region_H1, stimulation)

H3_COM_COM= permutation_test(COM_COM_H3, Nucleotide_region_H3, stimulation)
H3_AUT_AUT= permutation_test(AUT_AUT_H3, Nucleotide_region_H3, stimulation)
H3_LAIV_NONE= permutation_test(LAIV_NONE_H3, Nucleotide_region_H3, stimulation)
H3_LAIV_COM= permutation_test(LAIV_COM_H3, Nucleotide_region_H3, stimulation)
H3_NO_VAC= permutation_test(NO_VAC_H3, Nucleotide_region_H3, stimulation)
H3_COM_AUT= permutation_test(COM_AUT_H3, Nucleotide_region_H3, stimulation)
# run the method for 70% amino acid sites.
stimulation = 100000
H1_70_COM_COM= permutation_test(COM_COM_H1, Nucleotide_region_70_H1, stimulation)    
H1_70_AUT_AUT= permutation_test(AUT_AUT_H1, Nucleotide_region_70_H1, stimulation)
H1_70_LAIV_NONE= permutation_test(LAIV_NONE_H1, Nucleotide_region_70_H1, stimulation)    
H1_70_LAIV_COM= permutation_test(LAIV_COM_H1, Nucleotide_region_70_H1, stimulation)
H1_70_NO_VAC= permutation_test(NO_VAC_H1, Nucleotide_region_70_H1, stimulation)

H3_70_COM_COM= permutation_test(COM_COM_H3, Nucleotide_region_70_H3, stimulation)
H3_70_AUT_AUT= permutation_test(AUT_AUT_H3, Nucleotide_region_70_H3, stimulation)
H3_70_LAIV_NONE= permutation_test(LAIV_NONE_H3, Nucleotide_region_70_H3, stimulation)
H3_70_LAIV_COM= permutation_test(LAIV_COM_H3, Nucleotide_region_70_H3, stimulation)
H3_70_NO_VAC= permutation_test(NO_VAC_H3, Nucleotide_region_70_H3, stimulation)
H3_70_COM_AUT= permutation_test(COM_AUT_H3, Nucleotide_region_70_H3, stimulation)

# H1 convert dataset.
df_H1_NO_VAC = pd.DataFrame({"whole_genome":H1_NO_VAC})
df_H1_NO_VAC = df_H1_NO_VAC.reset_index()
df_H1_NO_VAC ["Group"]= "NO_VAC"  

df_H1_LAIV_NONE = pd.DataFrame({"whole_genome":H1_LAIV_NONE})
df_H1_LAIV_NONE = df_H1_LAIV_NONE.reset_index()
df_H1_LAIV_NONE ["Group"]= "LAIV_NONE"  

df_H1_protein=pd.concat([df_H1_LAIV_NONE,df_H1_NO_VAC])
#H3 convert dataset.
df_H3_COM_COM = pd.DataFrame({"whole_genome":H3_COM_COM})
df_H3_COM_COM = df_H3_COM_COM.reset_index()
df_H3_COM_COM ["Group"]= "COM_COM"
df_H3_COM_AUT = pd.DataFrame({"whole_genome":H3_COM_AUT})
df_H3_COM_AUT = df_H3_COM_AUT.reset_index()
df_H3_COM_AUT ["Group"]= "COM_AUT"
df_H3_LAIV_NONE = pd.DataFrame({"whole_genome":H3_LAIV_NONE})
df_H3_LAIV_NONE = df_H3_LAIV_NONE.reset_index()
df_H3_LAIV_NONE ["Group"]= "LAIV_NONE"
df_H3_LAIV_COM = pd.DataFrame({"whole_genome":H3_LAIV_COM})
df_H3_LAIV_COM = df_H3_LAIV_COM.reset_index()
df_H3_LAIV_COM ["Group"]= "LAIV_COM"
df_H3_NO_VAC = pd.DataFrame({"whole_genome":H3_NO_VAC})
df_H3_NO_VAC = df_H3_NO_VAC.reset_index()
df_H3_NO_VAC ["Group"]= "NO_VAC"

df_H3_protein=pd.concat([df_H3_COM_COM,df_H3_COM_AUT,df_H3_LAIV_NONE,df_H3_LAIV_COM,df_H3_NO_VAC])

# H1 convert dataset for 70% sites.
df_H1_70_NO_VAC = pd.DataFrame({"whole_genome":H1_70_NO_VAC})
df_H1_70_NO_VAC = df_H1_70_NO_VAC.reset_index()
df_H1_70_NO_VAC ["Group"]= "NO_VAC"  

df_H1_70_LAIV_NONE = pd.DataFrame({"whole_genome":H1_70_LAIV_NONE})
df_H1_70_LAIV_NONE = df_H1_70_LAIV_NONE.reset_index()
df_H1_70_LAIV_NONE ["Group"]= "LAIV_NONE"  

df_H1_70_protein=pd.concat([df_H1_70_LAIV_NONE,df_H1_70_NO_VAC])
#H3 convert dataset for 70% sites.
df_H3_70_COM_COM = pd.DataFrame({"whole_genome":H3_70_COM_COM})
df_H3_70_COM_COM = df_H3_70_COM_COM.reset_index()
df_H3_70_COM_COM ["Group"]= "COM_COM"
df_H3_70_COM_AUT = pd.DataFrame({"whole_genome":H3_70_COM_AUT})
df_H3_70_COM_AUT = df_H3_70_COM_AUT.reset_index()
df_H3_70_COM_AUT ["Group"]= "COM_AUT"
df_H3_70_LAIV_NONE = pd.DataFrame({"whole_genome":H3_70_LAIV_NONE})
df_H3_70_LAIV_NONE = df_H3_70_LAIV_NONE.reset_index()
df_H3_70_LAIV_NONE ["Group"]= "LAIV_NONE"
df_H3_70_LAIV_COM = pd.DataFrame({"whole_genome":H3_70_LAIV_COM})
df_H3_70_LAIV_COM = df_H3_70_LAIV_COM.reset_index()
df_H3_70_LAIV_COM ["Group"]= "LAIV_COM"
df_H3_70_NO_VAC = pd.DataFrame({"whole_genome":H3_70_NO_VAC})
df_H3_70_NO_VAC = df_H3_70_NO_VAC.reset_index()
df_H3_70_NO_VAC ["Group"]= "NO_VAC"

df_H3_70_protein=pd.concat([df_H3_70_COM_COM,df_H3_70_COM_AUT,df_H3_70_LAIV_NONE,df_H3_70_LAIV_COM,df_H3_70_NO_VAC])

#Export data.
df_H3_protein.to_csv(r'/Users/apple/Desktop/H3_Permutation_protein.csv',index=False)
df_H1_protein.to_csv(r'/Users/apple/Desktop/H1_Permutation_protein.csv',index=False)

#Export data for 70% sites.
df_H3_70_protein.to_csv(r'/Users/apple/Desktop/H3_Permutation_70%_protein.csv',index=False)
df_H1_70_protein.to_csv(r'/Users/apple/Desktop/H1_Permutation_70%_protein.csv',index=False)
