#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 16:41:40 2020

@author: Chong Li
"""

import random
import collections
import pandas as pd
import seaborn as sns
from random import sample
from collections import Counter

    
#Generate the snp file for H1 virus.    

with open("SNP_H1_nucleotide.csv", "r") as SNP_H1_nucleotide:
    

    Sample_H1=[]
    Group_H1=[]
    PB2_H1=[]
    PB1_H1=[]
    PA_H1=[]
    NP_H1=[]
    HA_H1=[]
    NA_H1=[]
    M_H1=[]
    NS_H1=[]
    title_H1=SNP_H1_nucleotide.readline()

    for line_H1 in SNP_H1_nucleotide:
        line_H1=line_H1.strip()
        split_file_H1=line_H1.split(sep=',')
        Sample_H1.append(split_file_H1[0])
        Group_H1.append(split_file_H1[2])
        PB2_H1.append(split_file_H1[5])
        PB1_H1.append(split_file_H1[6])
        PA_H1.append(split_file_H1[7])
        HA_H1.append(split_file_H1[8])
        NP_H1.append(split_file_H1[9])
        NA_H1.append(split_file_H1[10])
        M_H1.append(split_file_H1[11])
        NS_H1.append(split_file_H1[12])

COM_COM_H1={}
AUT_AUT_H1={}
LAIV_NONE_H1={}
LAIV_COM_H1={}
NO_VAC_H1={}

for index_H1, words_H1 in enumerate (Group_H1):
    if words_H1 == "1":
       COM_COM_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M":[M_H1[int(index_H1)]], "NS":[NS_H1[int(index_H1)]]} 
    if words_H1 == "2":
       AUT_AUT_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M":[M_H1[int(index_H1)]], "NS":[NS_H1[int(index_H1)]]} 
    if words_H1 == "5":
       LAIV_NONE_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M":[M_H1[int(index_H1)]], "NS":[NS_H1[int(index_H1)]]} 
    if words_H1 == "6":
       LAIV_COM_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M":[M_H1[int(index_H1)]], "NS":[NS_H1[int(index_H1)]]} 
    if words_H1 == "7":
       NO_VAC_H1[Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M":[M_H1[int(index_H1)]], "NS":[NS_H1[int(index_H1)]]}  
    
# Generate the snp files for H3 virus.
with open("SNP_H3_nucleotide.csv", "r") as SNP_H3_nucleotide:
    

    Sample_H3=[]
    Group_H3=[]
    PB2_H3=[]
    PB1_H3=[]
    PA_H3=[]
    NP_H3=[]
    HA_H3=[]
    NA_H3=[]
    M_H3=[]
    NS_H3=[]
    title_H3=SNP_H3_nucleotide.readline()

    for line_H3 in SNP_H3_nucleotide:
        line_H3=line_H3.strip()
        split_file_H3=line_H3.split(sep=',')
        Sample_H3.append(split_file_H3[0])
        Group_H3.append(split_file_H3[2])
        PB2_H3.append(split_file_H3[5])
        PB1_H3.append(split_file_H3[6])
        PA_H3.append(split_file_H3[7])
        HA_H3.append(split_file_H3[8])
        NP_H3.append(split_file_H3[9])
        NA_H3.append(split_file_H3[10])
        M_H3.append(split_file_H3[11])
        NS_H3.append(split_file_H3[12])

COM_COM_H3={}
AUT_AUT_H3={}
LAIV_NONE_H3={}
LAIV_COM_H3={}
NO_VAC_H3={}
COM_AUT_H3={}

for index_H3, words_H3 in enumerate (Group_H3):
    if words_H3 == "1":
       COM_COM_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M":[M_H3[int(index_H3)]], "NS":[NS_H3[int(index_H3)]]} 
    if words_H3 == "2":
       AUT_AUT_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M":[M_H3[int(index_H3)]], "NS":[NS_H3[int(index_H3)]]}  
    if words_H3 == "5":
       LAIV_NONE_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M":[M_H3[int(index_H3)]], "NS":[NS_H3[int(index_H3)]]}  
    if words_H3 == "6":
       LAIV_COM_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M":[M_H3[int(index_H3)]], "NS":[NS_H3[int(index_H3)]]}  
    if words_H3 == "7":
       NO_VAC_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M":[M_H3[int(index_H3)]], "NS":[NS_H3[int(index_H3)]]}  
    if words_H3 == "4":
       COM_AUT_H3[Sample_H3[int(index_H3)]]= {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M":[M_H3[int(index_H3)]], "NS":[NS_H3[int(index_H3)]]}  

    
#Generate the set of the numbers of nucleotides in each gene segment.
Nucleotide_region_H1={}
Nucleotide_region_H3={}

for num_H1 in range(0,len(Sample_H1),1):
    Nucleotide_region_H1[Sample_H1[int(num_H1)]]={"PB2":2280, "PB1":2274, "PA":2151, "HA":1701, "NP":1497, "NA":1410, "M":982, "NS":838}
for num_H3 in range(0,len(Sample_H3),1):
    Nucleotide_region_H3[Sample_H3[int(num_H3)]]={"PB2":2280, "PB1":2274, "PA":2151, "HA":1701, "NP":1497, "NA":1410, "M":982, "NS":838}
### With the 70% available region for mutation (about 30% mutation are lethal)

Nucleotide_70_region_H1={}
Nucleotide_70_region_H3={}

for num_H1_70 in range(0,len(Sample_H1),1):
    Nucleotide_70_region_H1[Sample_H1[int(num_H1_70)]]={"PB2":1596, "PB1":1592, "PA":1506, "HA":1191, "NP":1048, "NA":987, "M":687, "NS":587}
for num_H3_70 in range(0,len(Sample_H3),1):
    Nucleotide_70_region_H3[Sample_H3[int(num_H3_70)]]={"PB2":1596, "PB1":1592, "PA":1506, "HA":1191, "NP":1048, "NA":987, "M":687, "NS":587}


#build the function
def permutation_test(SNP_array, Nucleotide_region, stimulation):
    
    total_shared_sites=[]
    iteration_results = {"PB2":[], "PB1":[], "PA":[], "HA":[], "NP":[],"NA":[], "M":[], "NS":[]}
    for a in range(0, stimulation, 1):
        temporary_results = {"PB2":[], "PB1":[], "PA":[], "HA":[], "NP":[],"NA":[], "M":[], "NS":[]}
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
        shared_sites = iteration_results["PB2"][c]+iteration_results["PB1"][c]+iteration_results["PA"][c]+iteration_results["HA"][c]+iteration_results["NP"][c]+iteration_results["NA"][c]+iteration_results["M"][c]+iteration_results["NS"][c]
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
# run the methods for 70% mutation area
stimulation = 100000
H1_70_COM_COM= permutation_test(COM_COM_H1, Nucleotide_70_region_H1, stimulation)    
H1_70_AUT_AUT= permutation_test(AUT_AUT_H1, Nucleotide_70_region_H1, stimulation)
H1_70_LAIV_NONE= permutation_test(LAIV_NONE_H1, Nucleotide_70_region_H1, stimulation)    
H1_70_LAIV_COM= permutation_test(LAIV_COM_H1, Nucleotide_70_region_H1, stimulation)
H1_70_NO_VAC= permutation_test(NO_VAC_H1, Nucleotide_70_region_H1, stimulation)

H3_70_COM_COM= permutation_test(COM_COM_H3, Nucleotide_70_region_H3, stimulation)
H3_70_AUT_AUT= permutation_test(AUT_AUT_H3, Nucleotide_70_region_H3, stimulation)
H3_70_LAIV_NONE= permutation_test(LAIV_NONE_H3, Nucleotide_70_region_H3, stimulation)
H3_70_LAIV_COM= permutation_test(LAIV_COM_H3, Nucleotide_70_region_H3, stimulation)
H3_70_NO_VAC= permutation_test(NO_VAC_H3, Nucleotide_70_region_H3, stimulation)
H3_70_COM_AUT= permutation_test(COM_AUT_H3, Nucleotide_70_region_H3, stimulation)

# H1 convert dataset.
df_H1_COM_COM = pd.DataFrame({"whole_genome":H1_COM_COM})
df_H1_COM_COM = df_H1_COM_COM.reset_index()
df_H1_COM_COM ["Group"]= "COM/COM"


df_H1_NO_VAC = pd.DataFrame({"whole_genome":H1_NO_VAC})
df_H1_NO_VAC = df_H1_NO_VAC.reset_index()
df_H1_NO_VAC ["Group"]= "NO_VAC"  

df_H1_LAIV_NONE = pd.DataFrame({"whole_genome":H1_LAIV_NONE})
df_H1_LAIV_NONE = df_H1_LAIV_NONE.reset_index()
df_H1_LAIV_NONE ["Group"]= "LAIV_NONE"  

df_H1_nucleotide=pd.concat([df_H1_COM_COM,df_H1_LAIV_NONE,df_H1_NO_VAC])
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

df_H3_nucleotide=pd.concat([df_H3_COM_COM,df_H3_COM_AUT,df_H3_LAIV_NONE,df_H3_LAIV_COM,df_H3_NO_VAC])

#Export data.
df_H3_nucleotide.to_csv(r'/Users/apple/Desktop/H3_permutation_result.csv',index=False)
df_H1_nucleotide.to_csv(r'/Users/apple/Desktop/H1_permutation_result.csv',index=False)

# H1 convert dataset (70%).
df_H1_70_COM_COM = pd.DataFrame({"whole_genome":H1_70_COM_COM})
df_H1_70_COM_COM = df_H1_70_COM_COM.reset_index()
df_H1_70_COM_COM ["Group"]= "COM/COM"


df_H1_70_NO_VAC = pd.DataFrame({"whole_genome":H1_70_NO_VAC})
df_H1_70_NO_VAC = df_H1_70_NO_VAC.reset_index()
df_H1_70_NO_VAC ["Group"]= "NO_VAC"  

df_H1_70_LAIV_NONE = pd.DataFrame({"whole_genome":H1_70_LAIV_NONE})
df_H1_70_LAIV_NONE = df_H1_70_LAIV_NONE.reset_index()
df_H1_70_LAIV_NONE ["Group"]= "LAIV_NONE"  

df_H1_70_nucleotide=pd.concat([df_H1_70_COM_COM,df_H1_70_LAIV_NONE,df_H1_70_NO_VAC])
#H3 convert dataset (70%).
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

df_H3_70_nucleotide=pd.concat([df_H3_70_COM_COM,df_H3_70_COM_AUT,df_H3_70_LAIV_NONE,df_H3_70_LAIV_COM,df_H3_70_NO_VAC])

#Export data (70%).
df_H3_70_nucleotide.to_csv(r'/Users/apple/Desktop/H3_permutation_70_result.csv',index=False)
df_H1_70_nucleotide.to_csv(r'/Users/apple/Desktop/H1_permutation_70_result.csv',index=False)



