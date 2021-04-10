#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 20:35:03 2020

@author: Chong Li
"""
import collections
from collections import Counter

with open("BALF_H1_protein.csv", "r") as BALF_H1_protein:
     
     Annotate_H1_protein=[]
     Mutation_H1_protein=[]
     Treatment_H1_protein=[]
     Protein_H1_protein=[]
     
     for line_H1 in BALF_H1_protein:
         line_H1=line_H1.strip()
         BALF_H1_protein=line_H1.split(sep=',')
         Annotate_H1_protein.append(BALF_H1_protein[22])
         Mutation_H1_protein.append(BALF_H1_protein[23])
         Treatment_H1_protein.append(BALF_H1_protein[20])
         Protein_H1_protein.append(BALF_H1_protein[24])
Total_H1={"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}     
LAIV_NONE_H1={"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
Prime_boost_H1={"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
NO_VAC_H1={"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]} 

for index_H1, words_H1 in enumerate (Treatment_H1_protein):
    for segment in LAIV_NONE_H1.keys():   
        if Protein_H1_protein[int(index_H1)] == segment and Annotate_H1_protein[int(index_H1)] != "Synonymous":
           Total_H1[segment].append(Mutation_H1_protein[int(index_H1)][1:-1]) 
           if words_H1 == "LAIV/NONE":
              LAIV_NONE_H1[segment].append(Mutation_H1_protein[int(index_H1)][1:-1])
           if words_H1 == "COM/COM":
              Prime_boost_H1[segment].append(Mutation_H1_protein[int(index_H1)][1:-1])
           if words_H1 == "AUT/AUT":
              Prime_boost_H1[segment].append(Mutation_H1_protein[int(index_H1)][1:-1])
           if words_H1 == "LAIV/COM":
              Prime_boost_H1[segment].append(Mutation_H1_protein[int(index_H1)][1:-1])
           if words_H1 == "NO VAC/CHA":
              NO_VAC_H1[segment].append(Mutation_H1_protein[int(index_H1)][1:-1])

for gene in Total_H1:
    print(gene+":"+str(Counter(Total_H1[gene])))
#Results: 13 actual shared mutations in all pigs.
for gene in Prime_boost_H1:
    print(gene+":"+str(Counter(Prime_boost_H1[gene])))
#Results: 0 actual shared mutations in Prime_boost pigs.
for gene in LAIV_NONE_H1:
    print(gene+":"+str(Counter(LAIV_NONE_H1[gene])))
#Results: 7 actual shared mutations in LAIV_NONE pigs.    
for gene in NO_VAC_H1:
    print(gene+":"+str(Counter(NO_VAC_H1[gene])))
#Results: 5 actual shared mutations in NO_VAC pigs.    
 
with open("BALF_H3_protein.csv", "r") as BALF_H3_protein:
     
     Annotate_H3_protein=[]
     Mutation_H3_protein=[]
     Treatment_H3_protein=[]
     Protein_H3_protein=[]
     
     for line_H3 in BALF_H3_protein:
         line_H3=line_H3.strip()
         BALF_H3_protein=line_H3.split(sep=',')
         Annotate_H3_protein.append(BALF_H3_protein[22])
         Mutation_H3_protein.append(BALF_H3_protein[23])
         Treatment_H3_protein.append(BALF_H3_protein[20])
         Protein_H3_protein.append(BALF_H3_protein[24])
     
Total_H3={"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]} 
LAIV_NONE_H3={"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
Prime_boost_H3={"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
NO_VAC_H3={"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]} 

for index_H3, words_H3 in enumerate (Treatment_H3_protein):
    for segment in LAIV_NONE_H3.keys():
        if Protein_H3_protein[int(index_H3)] == segment and Annotate_H3_protein[int(index_H3)] != "Synonymous":
           Total_H3[segment].append(Mutation_H3_protein[int(index_H3)][1:-1])
           if words_H3 == "LAIV/NONE":
              LAIV_NONE_H3[segment].append(Mutation_H3_protein[int(index_H3)][1:-1])
           if words_H3 == "COM/COM":
              Prime_boost_H3[segment].append(Mutation_H3_protein[int(index_H3)][1:-1])
           if words_H3 == "AUT/AUT":
              Prime_boost_H3[segment].append(Mutation_H3_protein[int(index_H3)][1:-1])
           if words_H3 == "LAIV/COM":
              Prime_boost_H3[segment].append(Mutation_H3_protein[int(index_H3)][1:-1])
           if words_H3 == "NO VAC/CHA":
              NO_VAC_H3[segment].append(Mutation_H3_protein[int(index_H3)][1:-1])
           if words_H3 == "COM/AUT":
              Prime_boost_H3[segment].append(Mutation_H3_protein[int(index_H3)][1:-1])

for gene in Total_H3:
    print(gene+":"+str(Counter(Total_H3[gene])))
#Results:8 shared mutations in all pigs.    

for gene in LAIV_NONE_H3:
    print(gene+":"+str(Counter(LAIV_NONE_H3[gene])))
#Results:4 shared mutations in LAIV_NONE pigs.

for gene in Prime_boost_H3:
    print(gene+":"+str(Counter(Prime_boost_H3[gene])))
#Results:0 shared mutations in Prime_boost pigs.

for gene in NO_VAC_H3:
    print(gene+":"+str(Counter(NO_VAC_H3[gene])))

#Results:2 shared mutations in NO_VAC pigs.


