#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 21:17:34 2020

@author: Chong Li
"""
with open("H1N1_SNV.csv", "r") as SNP_H1_nucleotide:
     VarFreq_H1=[]
     Annotate_H1=[]
     Protein_H1=[]
     Sample_H1=[]
     group_H1=[]
     title_H1=SNP_H1_nucleotide.readline()
     
     for line_H1 in SNP_H1_nucleotide:
         line_H1=line_H1.strip()
         split_file_H1=line_H1.split(sep=',')
         Protein_H1.append(split_file_H1[22])
         VarFreq_H1.append(split_file_H1[6])
         Annotate_H1.append(split_file_H1[23])
         Sample_H1.append(split_file_H1[0])
         group_H1.append(split_file_H1[20])

with open("failed_sequencing.csv", "r") as failed_sequencing:
     SampleID=[]
     PB2=[]
     PB1=[]
     PA=[]
     HA=[]
     NP=[]
     NA=[]
     M=[]
     NS=[]
     BALF=[]
     first_line=failed_sequencing.readline()
     
     for row in failed_sequencing:
         row=row.strip()
         Parse_file=row.split(sep=",")
         SampleID.append(Parse_file[0])
         PB2.append(Parse_file[1])
         PB1.append(Parse_file[2])
         PA.append(Parse_file[3])
         HA.append(Parse_file[4])
         NP.append(Parse_file[5])
         NA.append(Parse_file[6])
         M.append(Parse_file[7])
         NS.append(Parse_file[8])
         BALF.append(Parse_file[9])


H1_seq={}
H3_seq={}

for a, b in enumerate (BALF):
    if b == "H1":
       H1_seq[SampleID[a]]={"PB2":[PB2[int(a)]], "PB1":[PB1[int(a)]], "PB1-F2":[PB1[int(a)]], "PA":[PA[int(a)]], "PA-X":[PA[int(a)]], "HA":[HA[int(a)]], "NP":[NP[int(a)]], "NA":[NA[int(a)]], "M1":[M[int(a)]], "M2":[M[int(a)]], "NS1":[NS[int(a)]], "NS2":[NS[int(a)]]}
    else:
       H3_seq[SampleID[a]]={"PB2":[PB2[int(a)]], "PB1":[PB1[int(a)]], "PB1-F2":[PB1[int(a)]], "PA":[PA[int(a)]], "PA-X":[PA[int(a)]], "HA":[HA[int(a)]], "NP":[NP[int(a)]], "NA":[NA[int(a)]], "M1":[M[int(a)]], "M2":[M[int(a)]], "NS1":[NS[int(a)]], "NS2":[NS[int(a)]]} 


Syn_H1={}
Nonsyn_H1={}
Stop_H1={}


treatment={}
for num, pig in enumerate (Sample_H1):
    if pig not in treatment.keys():
       treatment[str(pig)]=(str(group_H1[num]))

Synsites={"PB2":570, "PB1":568.5, "PB1-F2":60, "PA":537.75, "PA-X":174.75, "HA":425.25, "NP":374.25, "NA":352.5, "M1":189.75, "M2":73.5, "NS1":165, "NS2":91.5}
Nonsynsites={"PB2":1641.6, "PB1":1637.28, "PB1-F2":172.8, "PA":1548.72, "PA-X":503.28, "HA":1224.72, "NP":1077.84, "NA":1015.2, "M1":546.48, "M2":211.68, "NS1":475.2, "NS2":263.52}
Stopsites={"PB2":68.4, "PB1":68.22, "PB1-F2":7.2, "PA":64.53, "PA-X":20.97, "HA":51.03, "NP":44.91, "NA":42.3, "M1":22.77, "M2":8.82, "NS1":19.8, "NS2":10.98}

for index_H1, words_H1 in enumerate (Annotate_H1):
    if words_H1 == "Synonymous":
       if Sample_H1[index_H1] not in Syn_H1.keys():
          Syn_H1[Sample_H1[index_H1]]= {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
       Syn_H1[Sample_H1[index_H1]][Protein_H1[index_H1]].append(VarFreq_H1[index_H1])
    
    elif words_H1 == "Non_Synonymous":
         if Sample_H1[index_H1] not in Nonsyn_H1.keys():
            Nonsyn_H1[Sample_H1[index_H1]]= {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
         Nonsyn_H1[Sample_H1[index_H1]][Protein_H1[index_H1]].append(VarFreq_H1[index_H1])
            
    else:
         if Sample_H1[index_H1] not in Stop_H1.keys():
            Stop_H1[Sample_H1[index_H1]]= {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}       
         Stop_H1[Sample_H1[index_H1]][Protein_H1[index_H1]].append(VarFreq_H1[index_H1])           

for pigs in H1_seq.keys():
    if pigs not in Syn_H1.keys():
       Syn_H1[pigs] = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
    if pigs not in Nonsyn_H1.keys():
       Nonsyn_H1[pigs] = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]} 
    if pigs not in Stop_H1.keys():
       Stop_H1[pigs] = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}

H1_Pig_ID=[]
H1_Mutation_type=[]
H1_Protein=[]
H1_evolution_rate=[]
H1_group=[]

for sample_1 in Syn_H1:
    whole_freq_1 = 0
    whole_sites_1 = 0    
    for gene_1 in Syn_H1[sample_1]:
        if H1_seq[sample_1][gene_1][0] != "F": 
           freq_1=Syn_H1[sample_1][gene_1]
           total_freq_1 = 0
           for i in freq_1:
               total_freq_1 += float(i)
           whole_freq_1 += total_freq_1
           whole_sites_1 += Synsites[gene_1]
           rate_1=total_freq_1/float(Synsites[gene_1])/7
           H1_Pig_ID.append(str(sample_1))
           H1_Mutation_type.append("Synonymous")
           H1_Protein.append(str(gene_1))
           H1_evolution_rate.append(float(rate_1))
           H1_group.append(treatment[sample_1])
    
    whole_rate_1= whole_freq_1/float(whole_sites_1)/7       
    H1_Pig_ID.append(str(sample_1))
    H1_Mutation_type.append("Synonymous")
    H1_Protein.append("Whole_genome")
    H1_evolution_rate.append(float(whole_rate_1))
    H1_group.append(treatment[sample_1])
        

for sample_2 in Nonsyn_H1:
    whole_freq_2 = 0
    whole_sites_2 = 0  
    for gene_2 in Nonsyn_H1[sample_2]:
        if H1_seq[sample_2][gene_2][0] != "F":
           freq_2=Nonsyn_H1[sample_2][gene_2]
           total_freq_2 = 0
           for ii in freq_2:
               total_freq_2 += float(ii)
           whole_freq_2 += total_freq_2
           whole_sites_2 += Nonsynsites[gene_2]    
           rate_2=total_freq_2/float(Nonsynsites[gene_2])/7
           H1_Pig_ID.append(str(sample_2))
           H1_Mutation_type.append("Non-Synonymous")
           H1_Protein.append(str(gene_2))
           H1_evolution_rate.append(float(rate_2))
           H1_group.append(treatment[sample_2])
     
    whole_rate_2= whole_freq_2/float(whole_sites_2)/7       
    H1_Pig_ID.append(str(sample_2))
    H1_Mutation_type.append("Non-Synonymous")
    H1_Protein.append("Whole_genome")
    H1_evolution_rate.append(float(whole_rate_2))
    H1_group.append(treatment[sample_2])      
           
        
for sample_3 in Stop_H1:
    whole_freq_3 = 0
    whole_sites_3 = 0 
    for gene_3 in Stop_H1[sample_3]:
        if H1_seq[sample_3][gene_3][0] != "F":
           freq_3=Stop_H1[sample_3][gene_3]
           total_freq_3 = 0
           for iii in freq_3:
               total_freq_3 += float(iii)
           whole_freq_3 += total_freq_3
           whole_sites_3 += Stopsites[gene_3]    
           rate_3=total_freq_3/float(Stopsites[gene_3])/7
           H1_Pig_ID.append(str(sample_3))
           H1_Mutation_type.append("Stop-gained")
           H1_Protein.append(str(gene_3))
           H1_evolution_rate.append(float(rate_3))
           H1_group.append(treatment[sample_3])
    whole_rate_3= whole_freq_3/float(whole_sites_3)/7       
    H1_Pig_ID.append(str(sample_3))
    H1_Mutation_type.append("Stop-gained")
    H1_Protein.append("Whole_genome")
    H1_evolution_rate.append(float(whole_rate_3))
    H1_group.append(treatment[sample_3])            
           
        
outfile_1="BALF_H1_evolution_rate.csv"
out_file_1=open(outfile_1,'w')

print("Pig_ID", "Group", "Mutation_type", "Protein", "Evolution_rate", sep=',', end='\n', file= out_file_1)

for i in range(0,len(H1_Protein),1):
    print(str(H1_Pig_ID[i]), str(H1_group[i]), str(H1_Mutation_type[i]), str(H1_Protein[i]), str(H1_evolution_rate[i]), sep=',', end='\n', file= out_file_1)      
out_file_1.close()        
        
# Same script applied on BALF H3 data.        
with open("H3N2_SNV.csv", "r") as SNP_H3_nucleotide:
     VarFreq_H3=[]
     Annotate_H3=[]
     Protein_H3=[]
     Sample_H3=[]
     group_H3=[]
     title_H3=SNP_H3_nucleotide.readline()
     
     for line_H3 in SNP_H3_nucleotide:
         line_H3=line_H3.strip()
         split_file_H3=line_H3.split(sep=',')
         Protein_H3.append(split_file_H3[22])
         VarFreq_H3.append(split_file_H3[6])
         Annotate_H3.append(split_file_H3[23])
         Sample_H3.append(split_file_H3[0])
         group_H3.append(split_file_H3[20])        
        
Syn_H3={}
Nonsyn_H3={}
Stop_H3={}

treatmentH3={}
for num_2, pig_2 in enumerate (Sample_H3):
    if pig_2 not in treatmentH3.keys():
       treatmentH3[str(pig_2)]=(str(group_H3[num_2]))

for index_H3, words_H3 in enumerate (Annotate_H3):
    if words_H3 == "Synonymous":
       if Sample_H3[index_H3] not in Syn_H3.keys():
          Syn_H3[Sample_H3[index_H3]]= {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
       Syn_H3[Sample_H3[index_H3]][Protein_H3[index_H3]].append(VarFreq_H3[index_H3])
    
    elif words_H3 == "Non_Synonymous":
         if Sample_H3[index_H3] not in Nonsyn_H3.keys():
            Nonsyn_H3[Sample_H3[index_H3]]= {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
         Nonsyn_H3[Sample_H3[index_H3]][Protein_H3[index_H3]].append(VarFreq_H3[index_H3])
            
    else:
         if Sample_H3[index_H3] not in Stop_H3.keys():
            Stop_H3[Sample_H3[index_H3]]= {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}       
         Stop_H3[Sample_H3[index_H3]][Protein_H3[index_H3]].append(VarFreq_H3[index_H3])           

for pigs2 in H3_seq.keys():
    if pigs2 not in Syn_H3.keys():
       Syn_H3[pigs2] = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
    if pigs2 not in Nonsyn_H3.keys():
       Nonsyn_H3[pigs2] = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]} 
    if pigs2 not in Stop_H3.keys():
       Stop_H3[pigs2] = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
        
H3_Pig_ID=[]
H3_Mutation_type=[]
H3_Protein=[]
H3_evolution_rate=[]
H3_group=[]

for sample_4 in Syn_H3:
    whole_freq_4 = 0
    whole_sites_4 = 0 
    for gene_4 in Syn_H3[sample_4]:
        if H3_seq[sample_4][gene_4][0] != "F":
           freq_4=Syn_H3[sample_4][gene_4]
           total_freq_4 = 0
           for iv in freq_4:
               total_freq_4 += float(iv)
           whole_freq_4 += total_freq_4
           whole_sites_4 += Synsites[gene_4]    
           rate_4=total_freq_4/float(Synsites[gene_4])/7
           H3_Pig_ID.append(str(sample_4))
           H3_Mutation_type.append("Synonymous")
           H3_Protein.append(str(gene_4))
           H3_evolution_rate.append(float(rate_4))
           H3_group.append(treatmentH3[sample_4])
    whole_rate_4= whole_freq_4/float(whole_sites_4)/7       
    H3_Pig_ID.append(str(sample_4))
    H3_Mutation_type.append("Synonymous")
    H3_Protein.append("Whole_genome")
    H3_evolution_rate.append(float(whole_rate_4))
    H3_group.append(treatmentH3[sample_4])      
        

for sample_5 in Nonsyn_H3:
    whole_freq_5 = 0
    whole_sites_5 = 0 
    for gene_5 in Nonsyn_H3[sample_5]:
        if H3_seq[sample_5][gene_5][0] != "F":
           freq_5=Nonsyn_H3[sample_5][gene_5]
           total_freq_5 = 0
           for v in freq_5:
               total_freq_5 += float(v)
           whole_freq_5 += total_freq_5
           whole_sites_5 += Nonsynsites[gene_5]    
           rate_5=total_freq_5/float(Nonsynsites[gene_5])/7
           H3_Pig_ID.append(str(sample_5))
           H3_Mutation_type.append("Non-Synonymous")
           H3_Protein.append(str(gene_5))
           H3_evolution_rate.append(float(rate_5))
           H3_group.append(treatmentH3[sample_5])
    whole_rate_5= whole_freq_5/float(whole_sites_5)/7       
    H3_Pig_ID.append(str(sample_5))
    H3_Mutation_type.append("Non-Synonymous")
    H3_Protein.append("Whole_genome")
    H3_evolution_rate.append(float(whole_rate_5))
    H3_group.append(treatmentH3[sample_5])        
        
for sample_6 in Stop_H3:
    whole_freq_6 = 0
    whole_sites_6 = 0 
    for gene_6 in Stop_H3[sample_6]:
        if H3_seq[sample_6][gene_6][0] != "F":
           freq_6=Stop_H3[sample_6][gene_6]
           total_freq_6 = 0
           for vi in freq_6:
               total_freq_6 += float(vi)
           whole_freq_6 += total_freq_6
           whole_sites_6 += Stopsites[gene_6]     
           rate_6=total_freq_6/float(Stopsites[gene_6])/7
           H3_Pig_ID.append(str(sample_6))
           H3_Mutation_type.append("Stop-gained")
           H3_Protein.append(str(gene_6))
           H3_evolution_rate.append(float(rate_6))
           H3_group.append(treatmentH3[sample_6])
    whole_rate_6= whole_freq_6/float(whole_sites_6)/7       
    H3_Pig_ID.append(str(sample_6))
    H3_Mutation_type.append("Stop-gained")
    H3_Protein.append("Whole_genome")
    H3_evolution_rate.append(float(whole_rate_6))
    H3_group.append(treatmentH3[sample_6])          

outfile_2="BALF_H3_evolution_rate.csv"
out_file_2=open(outfile_2,'w')

print("Pig_ID", "Group", "Mutation_type", "Protein", "Evolution_rate", sep=',', end='\n', file= out_file_2)

for ii in range(0,len(H3_Protein),1):
    print(str(H3_Pig_ID[ii]), str(H3_group[ii]), str(H3_Mutation_type[ii]), str(H3_Protein[ii]), str(H3_evolution_rate[ii]), sep=',', end='\n', file= out_file_2)      
out_file_2.close()        
     
