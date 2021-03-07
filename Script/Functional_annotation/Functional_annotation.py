#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 15:28:06 2020

@author: Chong Li
"""
import sys
#command line: python script.name inputfile outputfile
#eg: python Functional_annotation.py Summary_BALF_NS_H1_Annotated.csv Summary_BALF_NS_H1_Annotated_functional.csv
inputfile= sys.argv[1]
SNP_file=open(inputfile,'r', encoding="ISO-8859-1")

Sample=[]
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
Day=[]
Treatment=[]
Segment=[]
Annotate=[]
Annotate_S=[]
Mutation=[]
Mutation_S=[]
title=SNP_file.readline()

#Name the output file
outfile= sys.argv[2]
out_file=open(outfile,'w')
out_file.write(title)

for line_2 in SNP_file:
    line_2=line_2.strip()
    split_file_2=line_2.split(sep=',')
    Sample.append(split_file_2[0])
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
    Day.append(split_file_2[19])
    Treatment.append(split_file_2[20])
    Segment.append(split_file_2[21])
    Annotate.append(split_file_2[22])
    Annotate_S.append(split_file_2[23])
    Mutation.append(split_file_2[24])
    Mutation_S.append(split_file_2[25])
SNP_file.close()
name=str(Segment[0])
if "PB2" in name:
    Functional_relevant=[]
    Feature=[]
    Comment=[]
    Reference=[]
    Functional_relevant_S=[]
    Feature_S=[]
    Comment_S=[]
    Reference_S=[]
    for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Functional_relevant_S.append("-")
        Feature_S.append("-")
        Comment_S.append("-")
        Reference_S.append("-")
        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):
            
            if amino_site > 0 and amino_site < 36:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_PB1-interacting-region_1")
               Comment.append("N-terminal 35 amino acid region of PB2 forms an interface with the C-terminal three helix bundle of PB1 (aa 685-757)")
               Reference.append("Structural Insight Into the Essential PB1-PB2 Subunit Contact of the Influenza Virus RNA Polymerase")
            
            elif amino_site > 735 and amino_site < 740:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_nuclear-localization-motif_736")
               Comment.append("This region is required for nuclear localization of PB2")
               Reference.append("Two Signals Mediate Nuclear Localization of Influenza Virus (A/WSN/33) Polymerase Basic Protein 2/UniProtKB - P03428 (PB2_I34A1)") 
            
            elif amino_site == 735 or amino_site == 740 or amino_site ==741:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_nuclear-localization-motif_735")
               Comment.append("NA")
               Reference.append("Two Signals Mediate Nuclear Localization of Influenza Virus (A/WSN/33) Polymerase Basic Protein 2")
        
            elif amino_site > 483 and amino_site < 497:
               Functional_relevant.append("Yes") 
               Feature.append("Influenza A_PB2_nuclear-localization-motif_448")
               Comment.append("This region is required for nuclear localization of PB2 and deletion at this site exhibits a phenotype that remains bound to perinuclear membrane")
               Reference.append("Two Signals Mediate Nuclear Localization of Influenza Virus (A/WSN/33) Polymerase Basic Protein 2")
            
            elif amino_site > 751 and amino_site < 756:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_nuclear-localization-motif_752")
               Comment.append("This region is required for nuclear localization of PB2 and forms a classibal bipartite NLS with amino acids 736-739.")
               Reference.append("Structure and Nuclear Import Function of the C-terminal Domain of Influenza Virus Polymerase PB2 Subunit")
               
            elif amino_site > 319 and amino_site < 448 and amino_site!=361 and amino_site!=376:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_cap-binding-site_320")
               Comment.append("This fragment is a domain co-crystalized with m7GTP.")
               Reference.append("The Structural Basis for Cap Binding by Influenza Virus Polymerase Subunit PB2")
            
            elif amino_site > 447 and amino_site < 484:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_cap-binding-site_320/Influenza A_PB2_nuclear-localization-motif_448")
               Comment.append("This fragment is a domain co-crystalized with m7GTP./This region is required for nuclear localization of PB2 and deletion at this site exhibits a phenotype that remains bound to perinuclear membrane.")
               Reference.append("The Structural Basis for Cap Binding by Influenza Virus Polymerase Subunit PB2/Two Signals Mediate Nuclear Localization of Influenza Virus (A/WSN/33) Polymerase Basic Protein 2 ")

            elif amino_site == 361 or amino_site == 376:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_m7GTP-binding-sites_361")
               Comment.append("Residues E361 and K376 directly interact with the guanine base of the 5'-m7G cap of host mRNAs. The aromatic residues H357 and F404 support this interaction.")
               Reference.append("The Structural Basis for Cap Binding by Influenza Virus Polymerase Subunit PB2")
               
            elif amino_site == 714:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_determinants-of-host-range_701") 
               Comment.append("The 701N and 714R; together with PA 615N; and NP 319K; PB1 13P and 678N cause increase in polymerase activity and confers adaptation of avian influenza virus to the mammalian host.")
               Reference.append("16339318")
               
            elif amino_site == 271:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_Polymerase-activity-in-mammalian-cells_271")
               Comment.append("This residue is responsible for the enhanced polymerase activity in mammalian cells.")
               Reference.append("PB2 Residue 271 Plays a Key Role in Enhanced Polymerase Activity of Influenza A Viruses in Mammalian Host Cells")
               
            elif amino_site == 28 or amino_site == 274 or amino_site == 526 or amino_site == 553 or amino_site == 607:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_polymerase-activity_28")
               Comment.append("A/duck/Guangxi/53/2002 differed from A/duck/Fujian/01/2002 by Met28Ile; Ala274Thr; Lys526Arg; Ile553Val; Leu607Val mutations. A/duck/Guangxi/53/2002 showed reduced polymerase activity.")
               Reference.append("Correlation Between Polymerase Activity and Pathogenicity in Two Duck H5N1 Influenza Viruses Suggests That the Polymerase Contributes to Pathogenicity")
               
            elif amino_site == 256:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_polymerase-activity_256")
               Comment.append("Introduction of an Asp256Gly substitution in the A/chicken/Yamaguchi/7/2004 backbone conferred increased replication efficiency in pigs as measured by viral titers from nasal swabs 1 day pi and increased polymerase activity of the RNP expressed.")
               Reference.append("19052090")
               
            elif amino_site == 627:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_determinant-of-virulence_627/replication-efficiency_627/polymerase-activity_627/tissue-tropism_627/transmissibility_627")
               Comment.append("Introduction of the Lys627Glu substitution in the A/swan/Germany/R65/2006 backbone conferred decreased replication in mammalian cells using growth kinetics; decreased virulence in mice using lethal dose./ A/Vietnam/1203/2004 isolate possessing 627Lys compared to A/Vietnam/1204/2004 with 627Glu increased replicated systemically in mice. Introduction of the Glu627Lys substitution in the A/chicken/Yamaguchi/7/2004 backbone conferred increased polymerase activity of RNP expressed. Introduction of the Glu627Lys substitution in the A/Vietnam/1204/2004 backbone conferred increased replication in MDCK cells. A/Hong Kong/483/97 with the same mutation showed increased replication efficiency in cultured mouse astrocytes and LA-4 mouse lung adenoma cells. Introduction of the Glu627Lys substitution in the A/Indonesia/05/2005 backbone conferred increased airborne transmission in ferrets.")
               Reference.append("21849466/20016035/17922570/17521765/11546875/15016548/17098982/19052090/21846828/22723413/19119420/16533883/19393699")
               
            elif amino_site == 701:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_replication-efficiency_701/Influenza A_PB2_tissue-tropism_701/Influenza A_PB2_Transmissibility_627/Influenza A_PB2_determinants-of-host-range_701")
               Comment.append("Introduction of Asp701Asn substitution in the A/duck/Guangxi/22/2001 backbone conferred efficient replication in the nose; trachea and lung of guinea pigs at titer levels comparable to A/duck/Guangxi/35/2001./ Introduction Lys627Glu and Asp701Asn substitutions in the A/Vietnam/1203/2004 backbone conferred increased contact transmission amongst guinea pigs./The 701N and 714R; together with PA 615N; and NP 319K; PB1 13P and 678N cause increase in polymerase activity and confers adaptation of avian influenza virus to the mammalian host.")
               Reference.append("20041223/19264775/19119420/16140781/16339318")
            
            elif amino_site == 368 or amino_site == 391 or amino_site == 447:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_tissue-tropism_368/ determinant-of-virulence_368")
               Comment.append("Introduction of the substitutions Arg368Gln; Gln391Glu; Gln447His; Lys627Glu in the A/Vietnam/1203/2004 backbone conferred reduced virulence as indicated by lethality in mice and conferred histologic alteration in the lungs; liver and brain of ferrets.")
               Reference.append("16533883")
               
            elif amino_site == 89 or amino_site == 309 or amino_site == 339 or amino_site == 477 or amino_site == 495 or amino_site == 676:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_polymerase-activity_89")
               Comment.append("Introduction of Leu89Val; Gly309Asp; Thr339Lys; Arg477Gly; Ile495Val; Lys627Glu; Ala676Thr naturally occurring substitutions in the A/wild duck/Hunan/021/2005 backbone conferred increased polymerase activity in mouse cells.")
               Reference.append("19393699")
               
            elif amino_site == 591:
               Functional_relevant.append("Yes")
               Feature.append("Influenza A_PB2_replication-efficiency_591")
               Comment.append("Introduction of the Gln591Lys substitution from A/Indonesia/UT3006/05 into the A/chicken/Indonesia/UT3091/2005 backbone conferred increased replication in NHBE cells and increased virulence as indicated by lethality in mice.")
               Reference.append("20700447")
               
            else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA") 
            
if "PB1" in name:               
   Functional_relevant=[]
   Feature=[]
   Comment=[]
   Reference=[]             
   Functional_relevant_S=[]
   Feature_S=[]
   Comment_S=[]
   Reference_S=[]
   for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Amino_mutation_S=Mutation_S[order]
        amino_site_S=int(Amino_mutation_S[1:-1]) 
        
        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):
           
           if amino_site > 186 and amino_site < 196:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_nuclear-localization-motif_187")
              Comment.append("NA")
              Reference.append("Uniprot:P16511")
              
           elif amino_site > 202 and amino_site < 217:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_nuclear-localization-motif_203")
              Comment.append("NA")
              Reference.append("Uniprot:P16511")   
               
           elif amino_site > 249 and amino_site < 257:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_nuclear-localization-motif_249")
              Comment.append("NA")
              Reference.append("Uniprot:P16511")
              
           elif amino_site == 249:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_nuclear-localization-motif_249/Influenza A_PB1_promoter-binding-site_233")
              Comment.append("NA")
              Reference.append("Uniprot:P16511/16476991")
              
           elif amino_site == 233 or amino_site == 238 or amino_site == 239:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_promoter-binding-site_233")
              Comment.append("NA")
              Reference.append("16476991")
               
           elif amino_site > 0 and amino_site < 26 and amino_site != 3:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_PA-binding-region_1")
              Comment.append("The N-terminal region of PB1 interacts with the C-terminus of PA (residues 257-716). This subunit interface complex is essential for initiation of transcription.")
              Reference.append("18615018/8811014")
             
           elif amino_site > 684 and amino_site < 758 and amino_site != 678:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_PB2-binding-region_685")
              Comment.append("The C-terminal three helix bundle of PB1 binds to 1-37 and 1-86 fragments on the N-terminus of PB2. This interface is crucial for the regulation of overall enzyme activity")
              Reference.append("19461581/17967456/8811014")
              
           elif amino_site == 678 :   
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_determinants-of-host-range_13/Influenza A_PB1_PB2-binding-region_685")
              Comment.append("The PB1 13P and 678N; together with PB2 701N and 714R; PA 615N; and NP 319K cause a dramatic increase in polymerase activity and confer adaptation of avian influenza virus to the mammalian host./The C-terminal three helix bundle of PB1 binds to 1–37 and 1–86 fragments on the N-terminus of PB2. This interface is crucial for the regulation of overall enzyme activity")
              Reference.append("16339318/19461581/17967456/8811014")
              
           elif amino_site == 13 :   
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_determinants-of-host-range_13")
              Comment.append("The PB1 13P and 678N; together with PB2 701N and 714R; PA 615N; and NP 319K cause a dramatic increase in polymerase activity and confer adaptation of avian influenza virus to the mammalian host.")
              Reference.append("16339318")
              
           elif amino_site == 436:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_determinant-of-virulence_436")
              Comment.append("Introduction of Tyr436His substitution in the A/Vietnam/1203/2004 backbone conferred decreased virulence as indicated by the survival rate of mice.")
              Reference.append("17553873")
              
           elif amino_site == 99 or amino_site == 368:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_transmissibility_99")
              Comment.append("Introduction of His99Tyr and Ile368val naturally occurring substitutions in the A/Indonesia/5/2005 backbone conferred increased airborne transmission in mammals")
              Reference.append("22723413")
              
           elif amino_site == 207:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_clinical-symptoms-of-disease_207/Influenza A_PB1_determinant-of-virulence_207/Influenza A_PB1_polymerase-activity_207")
              Comment.append("Introduction of Lys207Arg substitution in the A/Vietnam/1203/2004 backbone conferred increased virulence as indicated by mortality in mallards. Clinical signs of disease observed in mallards: cloudy eyes; appeared depressed; neurological signs. Introduction of Lys207Arg substitution in the A/Vietnam/1203/2004 backbone conferred decreased polymerase activity as indicated by the luciferase activity.")
              Reference.append("17553873")
              
           elif amino_site == 328 or amino_site == 375:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_determinant-of-virulence_3/Influenza A_PB1_tissue-tropism_3")
              Comment.append("Introduction of Val3Ala; Asn328Lys; Asn375Ser substitutions in the A/Vietnam/1203/2004 backbone conferred increased virulence as indicated by lethality in mice.")
              Reference.append("16533883")
            
           elif amino_site == 473 or amino_site == 598:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_polymerase-activity_473/Influenza A_PB1_replication-efficiency_473")
              Comment.append("Introduction of Val473Leu and Pro598Leu substitutions in the recombinant virus A/Cambodia/P0322095/2005 (PB1; PB2; PA; NP) x WSN conferred decreased replication in MDCK and A549 cells")
              Reference.append("22090209")
              
           elif amino_site == 3 :
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_PB1_determinant-of-virulence_3/Influenza A_PB1_tissue-tropism_3/Influenza A_PB1_PA-binding-region_1") 
              Comment.append("Introduction of Val3Ala; Asn328Lys; Asn375Ser substitutions in the A/Vietnam/1203/2004 backbone conferred increased virulence as indicated by lethality in mice./The N-terminal region of PB1 interacts with the C-terminus of PA (residues 257-716). This subunit interface complex is essential for initiation of transcription.")
              Reference.append("16533883/18615018/8811014")
          
           else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA")   
    
   
        
        if str(Amino_mutation_S[0]) != str(Amino_mutation_S[-1]):
           if amino_site_S == 35:
             Functional_relevant_S.append("Yes")
             Feature_S.append("Influenza A_PB1-F2_PKC-phosphorylation-site_35") 
             Comment_S.append("PKC phosphorylation site is related to proapoptotic activity of the PB1-F2 in PR8 strains.")
             Reference_S.append("19523156")
              
           elif amino_site_S > 61 and amino_site_S < 71 and amino_site_S != 66:
             Functional_relevant_S.append("Yes")
             Feature_S.append("Influenza A_PB1-F2_determinants-of-disease-progression_62/Influenza A_PB1-F2_determinant-of-secondary-bacterial-infections_61")
             Comment_S.append("Immunizing mice with a synthetic peptide corresponding to positions 62-70 of PB1-F2 is recognized by CD8+ T-lymphocytes that lyse virus-infected cells./Enhancement of secondary bacterial pneumonia in mice is seen by pretreatment with PB1-F2 C-terminus derived peptides ")
             Reference_S.append("11726970/18005742")
              
           elif amino_site_S > 70 and amino_site_S < 87 :
             Functional_relevant_S.append("Yes")
             Feature_S.append("Influenza A_PB1-F2_determinant-of-secondary-bacterial-infections_61")
             Comment_S.append("Enhancement of secondary bacterial pneumonia in mice is seen by pretreatment with PB1-F2 C-terminus derived peptides.")
             Reference_S.append("18005742")
              
           elif amino_site_S == 61:
             Functional_relevant_S.append("Yes")
             Feature_S.append("Influenza A_PB1-F2_determinant-of-secondary-bacterial-infections_61")
             Comment_S.append("Enhancement of secondary bacterial pneumonia in mice is seen by pretreatment with PB1-F2 C-terminus derived peptides.")
             Reference_S.append("18005742")
              
           elif amino_site_S == 87:
             Functional_relevant_S.append("Yes")
             Feature_S.append("Influenza A_PB1-F2_determinant-of-secondary-bacterial-infections_61/Influenza A_PB1-F2_putative-determinants-of-pathogenicity_51")
             Comment_S.append("Changes in these three residues of pathogenic H5N1 viruses affect pathogenicity in mallard ducks. A majority of H5N1 viruses have M51; V56; and E87 in PB1-F2./Enhancement of secondary bacterial pneumonia in mice is seen by pretreatment with PB1-F2 C-terminus derived peptides")
             Reference_S.append("18005742/20383540")
              
           elif amino_site_S == 51 or amino_site_S == 56:
             Functional_relevant_S.append("Yes")
             Feature_S.append("Influenza A_PB1-F2_putative-determinants-of-pathogenicity_51")
             Comment_S.append("Changes in these three residues of pathogenic H5N1 viruses affect pathogenicity in mallard ducks. A majority of H5N1 viruses have M51; V56; and E87 in PB1-F2.")
             Reference_S.append("20383540")
              
           elif amino_site_S == 66:
             Functional_relevant_S.append("Yes")
             Feature_S.append("Influenza A_PB1-F2_replication-efficiency_66/Influenza A_PB1-F2_clinical-symptoms-of-disease_66/Influenza A_PB1-F2_determinant-of-secondary-bacterial-infections_61/Influenza A_PB1-F2_determinants-of-disease-progression_62")
             Comment_S.append("Introduction of Asn66Ser substitution in the A/Hong Kong/156/1997 backbone conferred increased replication efficiency as indicated by growth kinetics in MDCK and lungs of mice. Mice inoculated with the mutant virus showed significant weight loss. Introduction of Asn66Ser substitution in the A/Vietnam/1203/2004 backbone conferred increased replication in CNS 8 days post infection using plaque assay./Introduction of Asn66Ser substitution in the A/Hong Kong/156/1997 backbone conferred increased replication efficiency as indicated by growth kinetics in MDCK and lungs of mice. Mice inoculated with the mutant virus showed significant weight loss. Introduction of Asn66Ser substitution in the A/Vietnam/1203/2004 backbone conferred increased replication in CNS 8 days post infection using plaque assay./Immunizing mice with a synthetic peptide corresponding to positions 62-70 of PB1-F2 is recognized by CD8+ T-lymphocytes that lyse virus-infected cells.")
             Reference_S.append("21852950/18005742/11726970")

           else:
              if "*" in Amino_mutation_S :
                 Functional_relevant_S.append("Yes")
                 Feature_S.append("NA")
                 Comment_S.append("NA")
                 Reference_S.append("NA")
                  
              else:
                 Functional_relevant_S.append("No")
                 Feature_S.append("NA")
                 Comment_S.append("NA")
                 Reference_S.append("NA")
                  
               
        if str(Amino_mutation_S[0]) == str(Amino_mutation_S[-1]):
           Functional_relevant_S.append("No")
           Feature_S.append("NA")
           Comment_S.append("NA")
           Reference_S.append("NA")                   
 
if "PA" in name:               
   Functional_relevant=[]
   Feature=[]
   Comment=[]
   Reference=[]             
   Functional_relevant_S=[]
   Feature_S=[]
   Comment_S=[]
   Reference_S=[]
   for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Amino_mutation_S=Mutation_S[order]
        amino_site_S=int(Amino_mutation_S[1:-1]) 

        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):
           if amino_site > 256 and amino_site < 717:
              
              if amino_site == 631:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_PA_determinant-of-virulence_631")
                 Comment.append("Passaging of HK156 viruses in mouse brain and embryonated eggs led to the selection of high and low virulent variants in the mice model respectively. These phenotypic changes are confered by changes in amino acids in HA residues (211); PB1 (456 and 712); NP (127) and NS1 (101) proteins; together with the PA (631).")
                 Reference.append("10873787")
                 
              elif amino_site == 356:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_PA_determinant-of-pathogenicity and replication")
                 Comment.append("Prevailing PA Mutation K356R in Avian Influenza H9N2 Virus Increases Mammalian Replication and Pathogenicity.")
                 Reference.append("27384648")
                 
              elif amino_site == 615:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_PA_determinant-of-host-range_615")
                 Comment.append("The PA 615N; together with PB2 701N; 714R; NP 319K; PB1 13P and 678N causes increase in polymerase activity and confers adaptation of avian influenza virus to the mammalian host.")
                 Reference.append("16339318")
              
              elif amino_site == 336:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_PA_Polymerase-activity-in-mammalian-cells_85")
                 Comment.append("These residues are responsible for the enhanced polymerase activity in mammalian cells.")
                 Reference.append("21561908")
              
              elif amino_site == 515:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_PA_determinant-of-virulence_515/Influenza A_PA_polymerase-activity_515")
                 Comment.append("Introduction of Thr515Ala substitutions in the A/Vietnam/1203/2004 backbone conferred decreased polymerase activity as indicated by the luciferase activity; caused no mortality in ducks.")
                 Reference.append("17553873")
                 
              elif amino_site == 266 or amino_site == 357: 
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_PA_determinant-of-virulence_149")
                 Comment.append("A/duck/Guangxi/53/2002 differed from duck/Fujian/01/2002 by Pro149Ser; Arg266His; Lys357Ile; Thr515Ser mutations. A/duck/Guangxi/53/2002 had limited lethality in mice.")
                 Reference.append("20211480")
                 
              elif amino_site > 258 and amino_site < 262:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_PA_Inflammatory-response_259")
                 Comment.append("Among 1918 WT and 1918 FS there is no apparent change in histopathological changes. Gene ontology analysis indicated that sequences showing differential expression between 1918 WT and 1918 FS infection were associated predominantly with inflammation or immune response; apoptosis; cell differentiation; tissue remodeling. Genes in CCR5 signaling in macrophages is more highly induced in 1918 FS and 1918 PTC mutants relative to WT.")
                 Reference.append("22745253")
                 
              else:
                 if "*" in Amino_mutation :
                    Functional_relevant.append("Yes")
                    Feature.append("NA")
                    Comment.append("NA")
                    Reference.append("NA")
                  
                 
                 else:   
                    Functional_relevant.append("Yes")
                    Feature.append("Influenza A_PA_PB1-binding-region_257")
                    Comment.append("The C-terminus of PA interacts with the N-terminal region of PB1 (residues 1-25). This subunit interface complex is essential for initiation of transcription.")
                    Reference.append("18615018")
                 
           elif amino_site == 85 or amino_site == 186:
                Functional_relevant.append("Yes")
                Feature.append("Influenza A_PA_Polymerase-activity-in-mammalian-cells_85/Influenza A_PA_Determinants-of-replication_57")
                Comment.append("These residues are responsible for the enhanced polymerase activity in mammalian cells./Multiple residues on PA of avian-origin pH1N1 influenza viruses suppress host cell protein synthesis during infection; allowing for preferential production of viral proteins.")
                Reference.append("21561908/23283952")
                
           elif amino_site == 57 or amino_site == 62 or amino_site == 65 or amino_site == 86 or amino_site == 91 or amino_site == 100 or amino_site == 114:
                Functional_relevant.append("Yes")     
                Feature.append("Influenza A_PA_Determinants-of-replication_57")
                Comment.append("Multiple residues on PA of avian-origin pH1N1 influenza viruses suppress host cell protein synthesis during infection; allowing for preferential production of viral proteins.")
                Reference.append("23283952")
           
           elif amino_site == 149:
                Functional_relevant.append("Yes")
                Feature.append("Influenza A_PA_determinant-of-virulence_149")
                Comment.append("A/duck/Guangxi/53/2002 differed from duck/Fujian/01/2002 by Pro149Ser; Arg266His; Lys357Ile; Thr515Ser mutations. A/duck/Guangxi/53/2002 had limited lethality in mice.")
                Reference.append("20211480")
                
           else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA")
           
        if str(Amino_mutation_S[0]) != str(Amino_mutation_S[-1]):
           if amino_site_S > 232 and amino_site_S < 253:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_PA-X_host adaption_virulence and transmission")
              Comment_S.append("Truncation of PA-X contributes to virulence and transmission of H3N8 and H3N2 canine influenza viruses in dogs/Truncation of C-terminal 20 amino acids in PA-X contributes to adaptation of swine influenza virus in pigs")
              Reference_S.append("32461313/26912401")
              
           else:
               if "*" in Amino_mutation_S :
                  Functional_relevant_S.append("Yes")
                  Feature_S.append("NA")
                  Comment_S.append("NA")
                  Reference_S.append("NA")
                  
               else:
                  Functional_relevant_S.append("No")
                  Feature_S.append("NA")
                  Comment_S.append("NA")
                  Reference_S.append("NA")
                  
        if str(Amino_mutation_S[0]) == str(Amino_mutation_S[-1]):
           Functional_relevant_S.append("No")
           Feature_S.append("NA")
           Comment_S.append("NA")
           Reference_S.append("NA")      

if "HA_H1" in name:
    Functional_relevant=[]
    Feature=[]
    Comment=[]
    Reference=[]
    Functional_relevant_S=[]
    Feature_S=[]
    Comment_S=[]
    Reference_S=[]
    for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Functional_relevant_S.append("-")
        Feature_S.append("-")
        Comment_S.append("-")
        Reference_S.append("-")
        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):
           if amino_site == 114 or amino_site == 152 or amino_site == 153 or amino_site == 154 or amino_site == 170:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H1_sialic-acid-binding-site_114")
              Comment.append("This region is commonly refered to as the 'receptor binding site' (RBS) and is the region of HA that determines its receptor (sialic acid) binding specificity. It is a shallow surface pocket made of conserved amino acid residues.")
              Reference.append("10.3923/jbs.2007.113.122 ")
           
           elif amino_site == 200 or amino_site == 207 or amino_site == 211 or amino_site == 212 or amino_site == 243 or amino_site == 244 or amino_site == 245 or amino_site == 246 or amino_site == 247:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H1_sialic-acid-binding-site_114")
              Comment.append("This region is commonly refered to as the 'receptor binding site' (RBS) and is the region of HA that determines its receptor (sialic acid) binding specificity. It is a shallow surface pocket made of conserved amino acid residues.")
              Reference.append("10.3923/jbs.2007.113.122 ") 
              
           elif amino_site == 159 or amino_site == 147:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H1_lysine-fence_147")
              Comment.append("Region at the base of HA binding site that helps anchor the N-acetylneuraminic acid (Neu5Ac) and galactose sugars of both alpha 2-3 and alpha 2-6 glycans on the host cell. Neu5Ac acts as a receptor for influenza viruses; allowing attachment to mucous cells via hemagglutinin.")
              Reference.append("19513050")
            
           elif amino_site > 0 and amino_site < 18:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H1_signal-peptide_1")
              Comment.append("NA")
              Reference.append("Uniprot:P03452")
              
           elif amino_site == 204:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H1_determinants-of-receptor-specificity_204")
              Comment.append("These residues are shown to confer specificity of HA to the human alpha2-6 sialylated glycan receptors.")
              Reference.append("19513050/14764886")
           
           elif amino_site == 239:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H1_determinants-of-receptor-specificity_204/Influenza A_H1_determinant-of-virulence_239")
              Comment.append("These residues are shown to confer specificity of HA to the human alpha2-6 sialylated glycan receptors.")
              Reference.append("19513050/14764886/20844044/21966421")
              
           elif amino_site == 391:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H1_determinant-of-virulence_391")
              Comment.append("NA")
              Reference.append("10426210")
              
           elif amino_site == 172:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H1_determinant-of-receptor-binding_172/Influenza A_H1_sialic-acid-binding-site_114")
              Comment.append("This region is commonly refered to as the 'receptor binding site' (RBS) and is the region of HA that determines its receptor (sialic acid) binding specificity. It is a shallow surface pocket made of conserved amino acid residues.")
              Reference.append("2033664/10.3923/jbs.2007.113.122")
               
           elif amino_site == 236:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H1_determinant-of-receptor-binding_236/Influenza A_H1_determinant-of-transmission_236/Influenza A_H1_lysine-fence_147")
              Comment.append("Introduction of Ile219Lys substitution in the glycan receptor-binding site (RBS) of HA in the A/California/04/09 backbone increases its human receptor-binding affinity on the glycan array and confers airborne transmission in ferrets./Region at the base of HA binding site that helps anchor the N-acetylneuraminic acid (Neu5Ac) and galactose sugars of both alpha 2-3 and alpha 2-6 glycans on the host cell. Neu5Ac acts as a receptor for influenza viruses; allowing attachment to mucous cells via hemagglutinin.")
              Reference.append("21407805/19513050")
              
           else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA")
                 
if "HA_H3" in name:
    Functional_relevant=[]
    Feature=[]
    Comment=[]
    Reference=[]
    Functional_relevant_S=[]
    Feature_S=[]
    Comment_S=[]
    Reference_S=[]
    for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Functional_relevant_S.append("-")
        Feature_S.append("-")
        Comment_S.append("-")
        Reference_S.append("-")
        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):
           if amino_site > 0 and amino_site < 16:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H3_signal-peptide_1")
              Comment.append("NA")
              Reference.append("Uniprot:P03437")
              
           elif amino_site == 114 or amino_site == 150 or amino_site == 151 or amino_site == 152 or amino_site == 153 or amino_site == 169:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H3_sialic-acid-binding-site_114")
              Comment.append("Site responsible for binding of virus to sialic acid receptors on the host cell surface")
              Reference.append("9454721/3374584")
              
           elif amino_site == 171 or amino_site == 199 or amino_site == 206 or amino_site == 210 or amino_site == 211:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H3_sialic-acid-binding-site_114")
              Comment.append("Site responsible for binding of virus to sialic acid receptors on the host cell surface")
              Reference.append("9454721/3374584")
              
           elif amino_site > 239 and amino_site < 245:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H3_sialic-acid-binding-site_114")
              Comment.append("Site responsible for binding of virus to sialic acid receptors on the host cell surface")
              Reference.append("9454721/3374584")
     
           elif amino_site == 154:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H3_sialic-acid-binding-site_114/Influenza A_H3_determinant-of-infectivity_154")
              Comment.append("Site responsible for binding of virus to sialic acid receptors on the host cell surface/A serine at amino acid 138 is essential but not solely sufficient for a high infectivity phenotype in primary swine respiratory epithelial cells (SRECs).")
              Reference.append("9454721/3374584/22353399/18329747") 

           elif amino_site == 140:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H3_determinant-of-infectivity_140")
              Comment.append("G124D exhibits high infectivity phenotype in swine respiratory epithelial cells (SRECs); preferentially binds polylactosamine glycans.")
              Reference.append("22353399")
              
           elif amino_site == 158:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H3_determinant-of-infectivity_158")
              Comment.append("G142E exhibits high infectivity phenotype in swine respiratory epithelial cells (SRECs).")
              Reference.append("22353399")
              
           elif amino_site == 501:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H3_determinants-of-virulence_234")
              Comment.append("The G218W (HA1) as well as the T156N(HA2) mutation enhance replication and virulence in vivo in mouse model. In addition to these; mutations of the PB2 (D701N) contributes to virulence and host range.")
              Reference.append("20702632")
              
           elif amino_site == 234:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_H3_determinants-of-virulence_234/Influenza A_H3_determinant-of-pathogenicity_234/Influenza A_H3_determinant-of-receptor-binding_234")
              Comment.append("The single mutation Gly218Glu causes high pathogenicity in mice at late passage 10 and also causes enhanced binding abilities to α2-3 and α2-6 sialic acid-linked receptors./The G218W (HA1) as well as the T156N(HA2) mutation enhance replication and virulence in vivo in mouse model. In addition to these; mutations of the PB2 (D701N) contributes to virulence and host range.")
              Reference.append("18983930/20702632")

           else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA")

if "NP" in name:
    Functional_relevant=[]
    Feature=[]
    Comment=[]
    Reference=[]
    Functional_relevant_S=[]
    Feature_S=[]
    Comment_S=[]
    Reference_S=[]
    for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Functional_relevant_S.append("-")
        Feature_S.append("-")
        Comment_S.append("-")
        Reference_S.append("-")
        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):
           if amino_site > 0 and amino_site < 181:
              if amino_site > 0 and amino_site < 162:
                 if amino_site > 0 and amino_site < 19:
                    if amino_site == 16 :
                       Functional_relevant.append("Yes")
                       Feature.append("Influenza A_NP_nuclear-localization-signal1_1/Influenza A_NP_RNA-binding-domain_1/Influenza A_NP_PB2-interaction-domain_1/Influenza A_NP_determinant-of-pathogenicity_16")
                       Comment.append("Wild-type A/PR/8 is highly pathogenic in mice; whereas the PR/8 with D16G is less lethal confirming that a single mutation in the N terminus of NP of the human PR/8 virus ignificantly decreases the pathogenicity of the virus in mice. Similarly; introduction the human-like G16D substitution into the NP of highly pathogenic A/Vietnam/1203/04 (H5N1) virus decreases lethality in mice./Also called unconventional nuclear localization signal motif. This site mediates nuclear import of the NP of influenza A virus.")
                       Reference.append("9770415/Uniprot:P03466/7745727/10438825/9621005/18058063")
                       
                    else:
                       if "*" in Amino_mutation :
                          Functional_relevant.append("Yes")
                          Feature.append("NA")
                          Comment.append("NA")
                          Reference.append("NA")
                          
                       else:   
                          Functional_relevant.append("Yes")
                          Feature.append("Influenza A_NP_nuclear-localization-signal1_1/Influenza A_NP_RNA-binding-domain_1/Influenza A_NP_PB2-interaction-domain_1")
                          Comment.append("Also called unconventional nuclear localization signal motif. This site mediates nuclear import of the NP of influenza A virus.")
                          Reference.append("9770415/Uniprot:P03466/7745727/10438825/9621005")
                 
                    
                 elif amino_site == 34:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_NP_determinant-of-temperature-sensitivity_34/Influenza A_NP_RNA-binding-domain_1/Influenza A_NP_PB2-interaction-domain_1")
                     Comment.append("Multiple loci confer ts phenotype to the vaccine strains from master donor virus (MDV) A/AA/6/60 used in FluMist: PB1 (K391E and E581G); PB2 (N265S); and NP (D34G). The PB1 (A661T) also contributes to the ts phenotype.")
                     Reference.append("12620793/7745727/10438825/9621005")
                     
                 elif amino_site == 99:    
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_NP_transmissibility_99/Influenza A_NP_RNA-binding-domain_1/Influenza A_NP_PB2-interaction-domain_1")
                     Comment.append("Introduction of Arg99Lys and Ser345Asn naturally occurring substitutions in the A/Indonesia/5/2005 backbone conferred airborne transmission in mammals.")
                     Reference.append("22723413/7745727/10438825/9621005")
                 
                 elif amino_site == 105:   
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_NP_determinant-of-pathogenicity_105/Influenza A_NP_determinant-of-host-specificity_105/Influenza A_NP_RNA-binding-domain_1/Influenza A_NP_PB2-interaction-domain_1")
                     Comment.append("Valine at position 105 in NP is critical for the high pathogenicity./Valine at position 105 of NP is found to be one of the determinants for adaptation of avian influenza viruses from ducks to chickens. A/chicken/Yamaguchi/7/2004 (H5N1) rapidly kills chicken without severe clinical signs while A/duck/Yokohama/aq10/2003 (H5N1) causes prologned death time with severe clinical signs.")
                     Reference.append("21123376/7745727/10438825/9621005")
                 
                 elif amino_site == 127:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_NP_determinant-of-virulence_127/Influenza A_NP_RNA-binding-domain_1/Influenza A_NP_PB2-interaction-domain_1")
                     Comment.append("Passaging of HK156 viruses in mouse brain and embryonated eggs led to the selection of high and low virulent variants in the mice model respectively. These phenotypic changes are confered by changes in amino acids in HA residues (211); PB1 (456 and 712); PA (631); NP (127) and NS1 (101) proteins.")
                     Reference.append("10873787/7745727/10438825/9621005")
                 else:
                     if "*" in Amino_mutation :
                        Functional_relevant.append("Yes")
                        Feature.append("NA")
                        Comment.append("NA")
                        Reference.append("NA")
                     
                     else:
                        Functional_relevant.append("Yes")
                        Feature.append("Influenza A_NP_RNA-binding-domain_1/Influenza A_NP_PB2-interaction-domain_1")
                        Comment.append("NA")
                        Reference.append("7745727/10438825/9621005")
                     
              else:
                 if "*" in Amino_mutation :
                    Functional_relevant.append("Yes")
                    Feature.append("NA")
                    Comment.append("NA")
                    Reference.append("NA")
                  
                 
                 else:   
                    Functional_relevant.append("Yes")
                    Feature.append("Influenza A_NP_RNA-binding-domain_1")
                    Comment.append("NA")
                    Reference.append("7745727/10438825")
            
                
                
           elif amino_site == 184:
                Functional_relevant.append("Yes")
                Feature.append("Influenza A_NP_determinant-of-replication_184")
                Comment.append("A change from alanine to a lysine at residue 184 of NP results in increased replication and pathogenicity of the viruses in chickens.")
                Reference.append("19475480")
               
           elif amino_site > 188 and amino_site < 359:
              if amino_site > 197 and amino_site < 217:
                 if amino_site == 214:
                    Functional_relevant.append("Yes")
                    Feature.append("Influenza A_NP_Determinants-of-morphology_214/Influenza A_NP_nuclear-localization-signal2_198/Influenza A_NP-NP-association-region_189")
                    Comment.append("Residues 214; 217; and 253 of Aichi NP play a critical role in influenza virus morphology; possibly through interaction with the M1 layer during virus budding/This region is a nuclear targeting motif that mediates transport of the cytoplasmic reporter protein into the nucleus.")
                    Reference.append("24335312/9770415/Uniprot:P03466/10405371")
                 else:
                    if "*" in Amino_mutation :
                       Functional_relevant.append("Yes")
                       Feature.append("NA")
                       Comment.append("NA")
                       Reference.append("NA")
                     
                    else:
                       Functional_relevant.append("Yes")
                       Feature.append("Influenza A_NP_nuclear-localization-signal2_198/Influenza A_NP-NP-association-region_189")
                       Comment.append("This region is a nuclear targeting motif that mediates transport of the cytoplasmic reporter protein into the nucleus.")
                       Reference.append("9770415/Uniprot:P03466/10405371") 
                  
              elif amino_site == 217 or amino_site == 253:
                   Functional_relevant.append("Yes")
                   Feature.append("Influenza A_NP_Determinants-of-morphology_214/Influenza A_NP-NP-association-region_189")
                   Comment.append("Residues 214; 217; and 253 of Aichi NP play a critical role in influenza virus morphology; possibly through interaction with the M1 layer during virus budding")
                   Reference.append("24335312/10405371") 
                
              elif amino_site > 254 and amino_site < 359:
                   if amino_site > 339 and amino_site < 359:
                      if amino_site == 345:
                         Functional_relevant.append("Yes")
                         Feature.append("Influenza A_NP_transmissibility_99/Influenza A_NP_PB2-binding-site_340/Influenza A_NP_PB2-interaction-domain_255/Influenza A_NP-NP-association-region_189")
                         Comment.append("Introduction of Arg99Lys and Ser345Asn naturally occurring substitutions in the A/Indonesia/5/2005 backbone conferred airborne transmission in mammals.")
                         Reference.append("22723413/9621005/10405371")
                      else:
                         if "*" in Amino_mutation :
                            Functional_relevant.append("Yes")
                            Feature.append("NA")
                            Comment.append("NA")
                            Reference.append("NA")
                          
                         
                         else:  
                            Functional_relevant.append("Yes")
                            Feature.append("IInfluenza A_NP_PB2-binding-site_340/Influenza A_NP_PB2-interaction-domain_255/Influenza A_NP-NP-association-region_189")
                            Comment.append("NA")
                            Reference.append("9621005/10405371")
                   elif amino_site == 319:
                       Functional_relevant.append("Yes")
                       Feature.append("Influenza A_NP_determinant-of-host-specificity_319/Influenza A_NP_species-adaptation_319/Influenza A_NP_PB2-interaction-domain_255/Influenza A_NP-NP-association-region_189")
                       Comment.append("Introduction of Asn319Lys naturally occurring substitution in the A/seal/Mass/1/1980(H7N7) backbone conferred increased binding to importin alpha1./The NP 319K; together with PB2 701N and 714R; PA 615N; PB1 13P and 678N causes increase in polymerase activity and confers adaptation of avian influenza virus to the mammalian host.")
                       Reference.append("16339318/18248089/10405371")
                       
                   elif amino_site == 267 or amino_site == 314 or amino_site == 330 or amino_site == 332:    
                       Functional_relevant.append("Yes")
                       Feature.append("Influenza A_NP_RNA-binding-domain_1/Influenza A_NP_PB2-interaction-domain_255/Influenza A_NP-NP-association-region_189")
                       Comment.append("NA")
                       Reference.append("7745727/10438825/10405371")
                   else:
                       if "*" in Amino_mutation :
                          Functional_relevant.append("Yes")
                          Feature.append("NA")
                          Comment.append("NA")
                          Reference.append("NA") 
                       
                       else:
                          Functional_relevant.append("Yes")
                          Feature.append("Influenza A_NP_PB2-interaction-domain_255/Influenza A_NP-NP-association-region_189")
                          Comment.append("NA")
                          Reference.append("9621005/10405371")
              else:
                  if "*" in Amino_mutation :
                     Functional_relevant.append("Yes")
                     Feature.append("NA")
                     Comment.append("NA")
                     Reference.append("NA") 
                  
                    
                  else:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_NP-NP-association-region_189")
                     Comment.append("NA")
                     Reference.append("10405371")
               
           elif amino_site > 358 and amino_site < 499:
                if amino_site == 479:
                   Functional_relevant.append("Yes")
                   Feature.append("Influenza A_NP-NP-association-region_479/Influenza A_NP_PB2-binding-site_340")
                   Comment.append("This residue is important for NP-NP interactions; as its alteration to alanine increases self-association several-fold.increased self-association")
                   Reference.append("10405371/9621005")
                    
                elif amino_site > 358 and amino_site < 466:
                   if amino_site > 370 and amino_site < 466:
                      if amino_site == 386:
                         Functional_relevant.append("Yes")
                         Feature.append("Influenza A_NP_RNA-binding-domain_1/Influenza A_NP-NP-association-region_371/Influenza A_NP_PB2-interaction-domain_255/Influenza A_NP_PB2-binding-site_340")
                         Comment.append("NA")
                         Reference.append("7745727/10438825/9621005/10405371")
                          
                      elif amino_site > 401 and amino_site < 429:
                         if amino_site == 412 or amino_site == 416:
                            Functional_relevant.append("Yes")
                            Feature.append("Influenza A_NP_RNA-binding-domain_1/Influenza A_NP-dimerization-domain_402/Influenza A_NP-NP-association-region_371/Influenza A_NP_PB2-interaction-domain_255/Influenza A_NP_PB2-binding-site_340")
                            Comment.append("This region of the tail loop is required for NP dimerization.")
                            Reference.append("7745727/10438825/17151603/9621005/10405371")
                         else:
                            if "*" in Amino_mutation :
                               Functional_relevant.append("Yes")
                               Feature.append("NA")
                               Comment.append("NA")
                               Reference.append("NA")
                            else: 
                               Functional_relevant.append("Yes")
                               Feature.append("Influenza A_NP-dimerization-domain_402/Influenza A_NP-NP-association-region_371/Influenza A_NP_PB2-interaction-domain_255/Influenza A_NP_PB2-binding-site_340")
                               Comment.append("This region of the tail loop is required for NP dimerization.")
                               Reference.append("17151603/9621005/10405371") 
                      
                      else:
                          if "*" in Amino_mutation :
                             Functional_relevant.append("Yes")
                             Feature.append("NA")
                             Comment.append("NA")
                             Reference.append("NA")
                          
                          else:
                             Functional_relevant.append("Yes")
                             Feature.append("Influenza A_NP-NP-association-region_371/Influenza A_NP_PB2-interaction-domain_255/Influenza A_NP_PB2-binding-site_340")
                             Comment.append("NA")
                             Reference.append("9621005/10405371")
                          
                   else:
                       if "*" in Amino_mutation :
                          Functional_relevant.append("Yes")
                          Feature.append("NA")
                          Comment.append("NA")
                          Reference.append("NA")
                       else:
                          Functional_relevant.append("Yes")
                          Feature.append("Influenza A_NP-NP-association-region_371/Influenza A_NP_PB2-binding-site_340")
                          Comment.append("NA")
                          Reference.append("9621005/10405371")
                        
                else:
                    if "*" in Amino_mutation :
                       Functional_relevant.append("Yes")
                       Feature.append("NA")
                       Comment.append("NA")
                       Reference.append("NA")
                    else:
                       Functional_relevant.append("Yes")
                       Feature.append("Influenza A_NP_PB2-binding-site_340")
                       Comment.append("NA")
                       Reference.append("9621005")
           else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA")        
                    

if "NA_H1" in name:
    Functional_relevant=[]
    Feature=[]
    Comment=[]
    Reference=[]
    Functional_relevant_S=[]
    Feature_S=[]
    Comment_S=[]
    Reference_S=[]
    for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Functional_relevant_S.append("-")
        Feature_S.append("-")
        Comment_S.append("-")
        Reference_S.append("-")
        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):                
           if amino_site == 324 :
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_calcium-binding-and-protein-stabilizing-site_324/Influenza A_N1_framework-region-of-active-site_149")
              Comment.append("Site at which calcium ions bind to NA and renders stability to the protein/The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure")
              Reference.append("Uniprot:Q9IGQ6/16912325/8497041")

           elif amino_site == 294 or amino_site == 298:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_calcium-binding-and-protein-stabilizing-site_294")
              Comment.append("Site at which calcium ions bind to NA and renders stability to NA")
              Reference.append("Uniprot:Q9IGQ6")
              
           elif amino_site == 344:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-binding-site/Influenza A_N1_calcium-binding-and-protein-stabilizing-site_344")
              Comment.append("Oseltamivir Binding Site/Zanamivir Binding Site/Site at which calcium ions bind to NA and renders stability to NA")
              Reference.append("Uniprot:Q9IGQ6/23028314/18480754")

           elif amino_site == 118:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-binding-site/Influenza A_N1_substrate-binding-and-hydrolysis-site_118")
              Comment.append("Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site/Region at which NA binds to its substrate and catalyzes its hydrolysis")               
              Reference.append("Uniprot:Q9IGQ6/22028647/18480754/16915235")
              
           elif amino_site == 293:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-binding-site/Influenza A_N1_substrate-binding-and-hydrolysis-site_293")
              Comment.append("Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site/Region at which NA binds to its substrate and catalyzes its hydrolysis")               
              Reference.append("Uniprot:Q9IGQ6/22028647/18480754/16915235")               

           elif amino_site == 368:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-binding-site/Influenza A_N1_substrate-binding-and-hydrolysis-site_368")
              Comment.append("Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site/Region at which NA binds to its substrate and catalyzes its hydrolysis")               
              Reference.append("Uniprot:Q9IGQ6/22028647/18480754/16915235")
  
           elif amino_site == 369:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_determinant-of-virulence_372")
              Comment.append("The five mutations (N369I in NA; K482R (silent mutation: G912A) in PB2; D538G in PB1; T139A (silent mutation: T121C) in M1 and W47G in HA2) control virulence and replicative capacity in mice.The PB1 and PB2 mutations are shown to be host restrictive in changing the virus to a mouse specific strain.")
              Reference.append("10426210")

           elif amino_site == 469:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_determinant-of-virulence_469/Influenza A_N1_determinants-of-pathogenicity_146/Influenza A_N1_determinant-of-replication_469")
              Comment.append("The C-terminal Lys at this position is critical for virulence and its ability to replicate in the mouse. It supports plasminogen-binding activity which is critical for WSN virus neurotropism./The presence of a C-terminal Lys 453 and the lack of glycosylation at position 130 (146 in N2 numbering) are required for binding of the NA to plasminogen./The C-terminal Lys at this position is critical for virulence and its ability to replicate in the mouse. It supports plasminogen-binding activity which is critical for WSN virus neurotropism.")
              Reference.append("11533192")
              
           elif amino_site == 146:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_determinant-of-virulence_146/Influenza A_N1_determinants-of-pathogenicity_146")
              Comment.append("The absence of a glycosylation site at position 130 of the NA plays a key role in the neurovirulence of WSN virus in mice./The presence of a C-terminal Lys 453 and the lack of glycosylation at position 130 (146 in N2 numbering) are required for binding of the NA to plasminogen.")
              Reference.append("8411368/11533192")
              
           elif amino_site > 48 and amino_site < 69:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_determinant-of-virulence_49")
              Comment.append("Introduction of a deletion known as the NA stalk motif in the A/chicken/Hubei/327/2004 backbone conferred increased virulence in mice and chicken indicated by the IVPI and lethality.")
              Reference.append("19225004/19609439")
              
           elif amino_site == 117:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_117")
              Comment.append("A/Chicken/Indonesia/Wates/77/2005 isolate with Ile97Val substitution conferred decreased sensitivity to oseltamivir using fluorescence based NA enzyme inhibition assay.")
              Reference.append("17112602/20523902/18836532")

           elif amino_site == 149:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_149/Influenza A_N1_antiviral-binding-site/Influenza A_N1_framework-region-of-active-site_149")
              Comment.append("Introduction of Val129Ala substitution in the A/CAM/408008/2005 backbone conferred decreased sensitivity to zanamivir using NA inhibition assay./Oseltamivir Binding Site/The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure")
              Reference.append("21343450/16912325/8497041/23028314")

           elif amino_site == 156:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_156")
              Comment.append("Introduction of Arg156Lys naturally occurring substitution in the A/Hong Kong/213/03 backbone conferred resistance to oseltamivir; peramivir and zanamivir using NA inhibition assay.")
              Reference.append("22379077")
              
           elif amino_site == 199:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_199")
              Comment.append("A/New York/4438/2009 isolate contained the Asp199Asn substitution that conferred decreased sensitivity to oseltamivir using NA inhibition assay.")
              Reference.append("21288815")
              
           elif amino_site == 278:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_278")
              Comment.append("Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site/Introduction of Glu258Gln naturally occurring substitution in the A/Vietnam/1203/2004 backbone decreased oseltamivir sensitivity using plaque reduction assay in MDCK cells.")
              Reference.append("22028647/18480754/16915235/17296744")
              
           elif amino_site == 223:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_223")
              Comment.append("Introduction of Ile203Val and His277Tyr substitutions in the A/Pennsylvania/30/2009 backbone conferred decreased sensitivity to oseltamivir, peramivir using NA inhibitors./Introduction of Ile203Met substitutions in the A/Chicken/Laos/26/2006 backbone conferred decreased sensitivity to oseltamivir using NA inhibition assay./Introduction of Ile203Arg and His277Tyr substitutions in the A/Pennsylvania/30/2009 backbone conferred decreased sensitivity to oseltamivir, zanamivir, peramivir using NA inhibitors.")
              Reference.append(":19651908/21148493/17302366/20016036/20858074")
              
           elif amino_site == 136:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_136")
              Comment.append("Introduction of Gln136Lys naturally occurring substitution in the A/Panama/1310/2008 backbone conferred reduced susceptibility to zanamivir and peramivir./Introduction of Gln116Leu naturally occurring substitution in the A/Vietnam/1203/2004 backbone conferred oseltamivir and zanamivir resistance using fluorescence based enzyme inhibition assay.")
              Reference.append("19917319/19641000/20603155")
              
           elif amino_site == 119:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_119")
              Comment.append("Introduction of Glu119Gly substitution in the A/Quebec/144147/09 backbone conferred resistance to oseltamivir; zanamivir and peramivir using NA inhibition assay./Introduction of Glu99Ala naturally occurring substitution in the A/Turkey/65 1242/06 backbone conferred resistance to oseltamivir; zanamivir and peramivir.")
              Reference.append("19651908/21148493/:22379077/17302366/20523902")
              
           elif amino_site == 275:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_275")
              Comment.append("Introduction of His255Tyr naturally occurring substitution in the A/Vietnam/1203/2004 backbone conferred decreased oseltamivir sensitivity as indicated by measuring inhibition of neuraminidase activity.")
              Reference.append("19651908/1170976/17296744/18368779/19022400/16228009/16371632")
              
           elif amino_site == 116:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_116")
              Comment.append("Introduction of Val95Ala substitution in the A/Turkey/15/2006 backbone conferred decreased sensitivity to oseltamivir and zanamivir using NA inhibition assay and measuring NA enzyme kinetics.")
              Reference.append("20016036/20523902/17112602")
              
           elif amino_site == 295:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_295")
              Comment.append("Asn275Ser substitution found in A/Egypt/1425 NAMRU3/2006 isolate conferred decreased oseltamivir sensitivity from patients treated with oseltamivir; increased replication in ferrets./Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site")
              Reference.append("20701864/21367898/19022400/21148493/17855542/22028647/18480754/16915235")
              
           elif amino_site == 247:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_247")
              Comment.append("A/chicken/Laos/13/08 isolate with Ser227Asn substitution conferred decreased oseltamivir sensitivity using NA inhibition assay/Zanamivir Binding Site/Laninamivir Binding Site")
              Reference.append("20016036/22028647/18715929")
              
           elif amino_site == 297:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-response_223")
              Comment.append("Introduction of Ile203Arg and His277Tyr substitutions in the A/Pennsylvania/30/2009 backbone conferred decreased sensitivity to oseltamivir; zanamivir; peramivir using NA inhibitors./Introduction of Ile203Met and His277Tyr substitutions in the A/Pennsylvania/30/2009 backbone conferred decreased sensitivity to oseltamivir; peramivir using NA inhibitors.")
              Reference.append("21148493/17302366/20016036/20858074")
              
           elif amino_site == 148 or amino_site == 181 or amino_site == 182 or amino_site == 254 or amino_site == 306 or amino_site == 322 or amino_site == 401 or amino_site == 436:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_catalytic-region-of-active-site_148")
              Comment.append("The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure")
              Reference.append("16912325/8497041")
            
           elif amino_site == 186 or amino_site == 208 or amino_site == 209 or amino_site == 228 or amino_site == 252 or amino_site == 257 or amino_site == 304 or amino_site == 307 or amino_site == 405:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_framework-region-of-active-site_149")
              Comment.append("The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure")
              Reference.append("16912325/8497041")
           
           elif amino_site == 151 or amino_site == 152 or amino_site == 179 or amino_site == 225 or amino_site == 277 or amino_site == 402:   
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-binding-site")
              Comment.append("Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site")
              Reference.append("22028647/18480754/16915235")
           
           elif amino_site == 228:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-binding-site")
              Comment.append("Zanamivir Binding Site/Laninamivir Binding Site")
              Reference.append("22028647/18480754/16915235")
              
           elif amino_site == 150:   
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N1_antiviral-binding-site")
              Comment.append("Oseltamivir Binding Site")
              Reference.append("22028647/18480754/16915235")
            
           else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA")     

if "NA_H3" in name:
    Functional_relevant=[]
    Feature=[]
    Comment=[]
    Reference=[]
    Functional_relevant_S=[]
    Feature_S=[]
    Comment_S=[]
    Reference_S=[]
    for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Functional_relevant_S.append("-")
        Feature_S.append("-")
        Comment_S.append("-")
        Reference_S.append("-")
        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):   
           if amino_site == 293 or amino_site == 297 or amino_site == 324:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_calcium-binding-site_293")
              Comment.append("NA")
              Reference.append("Uniprot:P06820")
           
           elif amino_site == 118 or amino_site == 292 or amino_site == 371:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_catalytic-region-of-active-site_118/Influenza A_N2_substrate-binding-site_118/Influenza A_N2_antiviral-binding-site")   
              Comment.append("The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure/Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site")
              Reference.append("Uniprot:P06820/23531861/16912325/8497041")

           elif amino_site == 116:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_antiviral-response_116")
              Comment.append("Viruses with mutation V116A show resistance to both inhibitors; zanamivir and oseltamivir carboxylate.")
              Reference.append("21253602")
            
           elif amino_site == 119:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_antiviral-response_119/Influenza A_N2_antiviral-binding-site/Influenza A_N2_framework-region-of-active-site_119")
              Comment.append("A/Texas/12/2007 isolate with substitution Glu119Ile conferred decreased sensitivity to oseltamivir; zanamivir using chemiluminescent NA inhibition assay./A double mutation which includes E119V and a change of isoleucine to valine at position 222 (E119V+I222V) confers oseltamivir resistance./The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure/Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site")
              Reference.append("16912325/8497041/23531861/18684820/17302366/17109288/16479508/20194700/21106781/21248368/16891631")
              
           elif amino_site == 136:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_antiviral-response_136")
              Comment.append("H3N2 viruses with a Q136K mutation in NA confer reduction in zanamivir susceptibility. These viruses also contained the S31N mutation in M2.")
              Reference.append("20202427")
                  
           elif amino_site == 198:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_antiviral-response_198/Influenza A_N2_framework-region-of-active-site_119")
              Comment.append("Introduction of Asp198Asn substitution in the A/New_York/1191/2009 backbone imputed with decreased sensitivity to oseltamivir as observed in N1 genome./The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure")
              Reference.append("16912325/8497041")
            
           elif amino_site == 274:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_antiviral-response_274/Influenza A_N2_framework-region-of-active-site_119")   
              Comment.append("A/Canada/83/2006 isolate with substitution His274Asn conferred mild decreased sensitivity to oseltamivir using NA inhibition assay./The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure")
              Reference.append("16912325/8497041/16891631")    

           elif amino_site == 222:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_antiviral-response_222/Influenza A_N2_framework-region-of-active-site_119")
              Comment.append("Introduction of Ile222Arg substitution in the A/British Columbia/EFA0401/2009 backbone imputed with decreased sensitivity to oseltamivir; zanamivir as observed in N1 genome./A double mutation which includes E119V and a change of isoleucine to valine at position 222 (E119V+I222V) confers oseltamivir resistance./The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure/Laninamivir Binding Site")
              Reference.append("16912325/8497041/18684820/22028647")

           elif amino_site == 151 or amino_site == 152 or amino_site == 224 or amino_site == 276 or amino_site == 406:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_catalytic-region-of-active-site_118/Influenza A_N2_antiviral-binding-site")
              Comment.append("The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure/Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site")
              Reference.append("16912325/8497041/22028647")
              
           elif amino_site == 179 or  amino_site == 425:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_framework-region-of-active-site_119")
              Comment.append("The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure")
              Reference.append("16912325/8497041")

           elif amino_site == 246:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_antiviral-binding-site")
              Comment.append("Oseltamivir Binding Site/Zanamivir Binding Site/Laninamivir Binding Site")
              Reference.append("22028647")
           elif amino_site == 156 or amino_site == 178 or amino_site == 227:   
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_framework-region-of-active-site_119/Influenza A_N2_antiviral-binding-site")
              Comment.append("The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure/Zanamivir Binding Site/Laninamivir Binding Site")
              Reference.append("16912325/8497041/22028647")
              
           elif amino_site == 277 or amino_site == 294:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_N2_framework-region-of-active-site_119/Influenza A_N2_antiviral-binding-site")
              Comment.append("The active site of Influenza NA is made of the catalytic site residues which interact with the sialic acid substrate and the framework site residues that are indirectly involved in suporting the catalytic site. (PMID: 16912325). Framework residues stabilize the active-site structure/Oseltamivir Binding Site/Laninamivir Binding Site")
              Reference.append("16912325/8497041/22028647") 
              
           elif amino_site == 347:
              Functional_relevant.append("Yes") 
              Feature.append("Influenza A_N2_antiviral-binding-site")
              Comment.append("Oseltamivir Binding Site")
              Reference.append("22028647") 
              
           else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA") 
           
if "M" in name:               
   Functional_relevant=[]
   Feature=[]
   Comment=[]
   Reference=[]             
   Functional_relevant_S=[]
   Feature_S=[]
   Comment_S=[]
   Reference_S=[]
   for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Amino_mutation_S=Mutation_S[order]
        amino_site_S=int(Amino_mutation_S[1:-1]) 
        
        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):
           
           if amino_site > 0 and amino_site < 165:
               
              if amino_site > 0 and amino_site < 77:
                  
                 if amino_site > 61 and amino_site < 69:
                     
                    Functional_relevant.append("Yes")
                    Feature.append("Influenza A_M1_lipid-binding-site_62/Influenza A_M1_lipid-binding-site_1/Influenza A_M1_RNP-binding-region_1/Influenza A_M1_membrane-binding-region_1")
                    Comment.append("This N-terminal region binds to vRNP.")
                    Reference.append("11222100/10438836/Uniprot:P03485")
                    
                 elif amino_site == 30:
                    Functional_relevant.append("Yes")
                    Feature.append("Influenza A_M1_determinants-of-virulence_30/Influenza A_M1_Determinants-of-morphology_30/Influenza A_M1_Determinants-of-budding_30/Influenza A_M1_Determinants-of-virus-production_30/Influenza A_M1_determinant-of-virulence_30/Influenza A_M1_RNP-binding-region_1/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                    Comment.append("These two M1 residues contribute to differences in virulence of H5N1 avian influenza viruses in mice. Positions 30 and 215 are found in the N-terminal and C-terminal domains of M1 respectively./Key to spherical morphology. Substitution of any one of these residues with the corresponding TR swine residues causes transition to a filamentous morphology./Making M1 protein sufficient for virus-like particle production; indicating a unique feature of the pH1N1 M1 protein to induce and complete budding at the plasma membrane by itself./Single replacement of either one of the residue in M1 protein reduces overall viral production as well as growth kinetics./Introduction of Asn30Asp and Thr215Ala substitutions in the A/duck/Guangxi/53/2002 backbone conferred increased virulence in mice indicated by survival rate./This N-terminal region binds to vRNP.")
                    Reference.append("19117585/23209789/10438836/11222100/Uniprot:P03485")
           
                 else:
                    if "*" in Amino_mutation :
                       Functional_relevant.append("Yes")
                       Feature.append("NA")
                       Comment.append("NA")
                       Reference.append("NA") 

                    else:
                       Functional_relevant.append("Yes")
                       Feature.append("Influenza A_M1_RNP-binding-region_1/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                       Comment.append("This N-terminal region binds to vRNP.")
                       Reference.append("10438836/Uniprot:P03485/11222100")
              elif amino_site == 90:
                  Functional_relevant.append("Yes")
                  Feature.append("Influenza A_M1_RNA-binding-domain_90/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                  Comment.append("NA")
                  Reference.append("10438836/Uniprot:P03485/11222100")
                  
                  
              elif amino_site > 90 and amino_site < 112:
                  
                  if amino_site > 90 and amino_site < 109:
                      
                     if amino_site > 101 and amino_site < 105:
                         
                        Functional_relevant.append("Yes")
                        Feature.append("Influenza A_M1_RNA-binding-domain_101/Influenza A_M1_NS2/NEP-interaction-site_101/Influenza A_M1_RNA-binding-domain_90/Influenza A_M1_transcription-inhibition-site_91/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                        Comment.append("This stretch of basic amino acids can bind to vRNA./This nuclear localization signal interacts with the C-terminus of NS2 (NEP)./This region has been found to be essential for anti-RNA synthesis activity; RNA binding; and oligomerization of M1.")
                        Reference.append("9225034/12970177/10438836/8523532/11222100/Uniprot:P03485")
                       
                     elif amino_site == 101 or amino_site == 105:
                        
                        Functional_relevant.append("Yes")
                        Feature.append("Influenza A_M1_determinants-of-temperature-sensitivity_101/Influenza A_M1_RNA-binding-domain_101/Influenza A_M1_NS2/NEP-interaction-site_101/Influenza A_M1_RNA-binding-domain_90/Influenza A_M1_transcription-inhibition-site_91/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                        Comment.append("The R101S-R105S double mutation (substitution of Ser for Arg at either position 101 or position 105 of the RKLKR domain) is synergistic and results in temperature sensitivity as seen by reduced viral replication at a restrictive temperature. The double mutation is seen to be fully attenuated in mice./This stretch of basic amino acids can bind to vRNA./This nuclear localization signal interacts with the C-terminus of NS2 (NEP)./This region has been found to be essential for anti-RNA synthesis activity; RNA binding; and oligomerization of M1.")
                        Reference.append("15650216/15331690/9225034/12970177/10438836/8523532/11222100/Uniprot:P03485")
                         
                     else:
                        if "*" in Amino_mutation :
                           Functional_relevant.append("Yes")
                           Feature.append("NA")
                           Comment.append("NA")
                           Reference.append("NA")
                           
                        else:
                           Functional_relevant.append("Yes")
                           Feature.append("Influenza A_M1_RNA-binding-domain_90/Influenza A_M1_transcription-inhibition-site_91/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                           Comment.append("This region has been found to be essential for anti-RNA synthesis activity; RNA binding; and oligomerization of M1.")
                           Reference.append("10438836/8523532/11222100/Uniprot:P03485")
                           
                           
                  else:
                      if "*" in Amino_mutation :
                           Functional_relevant.append("Yes")
                           Feature.append("NA")
                           Comment.append("NA")
                           Reference.append("NA")
                           
                      else:
                           Functional_relevant.append("Yes")
                           Feature.append("Influenza A_M1_transcription-inhibition-site_91/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                           Comment.append("This region has been found to be essential for anti-RNA synthesis activity; RNA binding; and oligomerization of M1.")
                           Reference.append("8523532/11222100/Uniprot:P03485")
              
                            
              elif amino_site > 113 and amino_site < 134:
                  Functional_relevant.append("Yes")
                  Feature.append("Influenza A_M1_lipid-binding-site_114/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                  Comment.append("NA")
                  Reference.append("7288926/11222100/11222100/Uniprot:P03485")
              
              elif amino_site > 134 and amino_site < 165:
                  
                  if amino_site >147 and amino_site < 163:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_M1_determinants-of-virulence_148/Influenza A_M1_RNA-binding-domain_135/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                     Comment.append("The putative zinc finger motif CCHH in the helix 9 of M1 plays a critical role in virulence in vivo in mice although it does not play an important role in virus growth in MDCK cells in cultures . The CCHH motif may be a contributory factor in virus virulence in a species-specific manner along with the remaining residues in the H9 region.")
                     Reference.append("16731908/10438836/11222100/Uniprot:P03485")
                     
                  elif amino_site == 139:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_M1_determinant-of-virulence_139/Influenza A_M1_determinant-of-replication_139/Influenza A_M1_RNA-binding-domain_135/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                     Comment.append("All five mutations (T139A (and silent mutation: T121C) of M1; D538G in PB1; K482R (silent mutation: G912A) in PB2; N369I in NA and W47G in HA2) control virulence and replicative capacity in mice.The PB1 and PB2 mutations are shown to be host restrictive in changing the virus to a mouse specific strain./The mouse adapted A/Fort Monmouth/1/47 virus contained Thr139Ala substitution that conferred increased virulence as indicated by measuring median lethal dose in mice and increased viral yield in lungs of mice and MDCK cells.")
                     Reference.append("10426210/8879138/10438836/11222100/Uniprot:P03485")
                     
                  elif amino_site == 142:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_M1_Determinants-of-virus-production_30/Influenza A_M1_RNA-binding-domain_135/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                     Comment.append("Single replacement of either one of the residue in M1 protein reduces overall viral production as well as growth kinetics.")
                     Reference.append("23209789/10438836/11222100/Uniprot:P03485")
                     
                  else:
                     if "*" in Amino_mutation :
                        Functional_relevant.append("Yes")
                        Feature.append("NA")
                        Comment.append("NA")
                        Reference.append("NA")
                        
                     else:
                        Functional_relevant.append("Yes")
                        Feature.append("Influenza A_M1_RNA-binding-domain_135/Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                        Comment.append("NA")
                        Reference.append("10438836/11222100/Uniprot:P03485")
                
              
              else:
                  if "*" in Amino_mutation :
                     Functional_relevant.append("Yes")
                     Feature.append("NA")
                     Comment.append("NA")
                     Reference.append("NA")
                       
                  else:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_M1_membrane-binding-region_1/Influenza A_M1_lipid-binding-site_1")
                     Comment.append("NA")
                     Reference.append("11222100/Uniprot:P03485")
                     
           elif amino_site == 165:
              Functional_relevant.append("Yes")
              Feature.append("Influenza A_M1_RNP-binding-region_165/Influenza A_M1_RNA-binding-domain_135")
              Comment.append("This C-terminal region binds to vRNP.")
              Reference.append("10438836/11222100/Uniprot:P03485")
              
           elif amino_site > 165 and amino_site < 253:
              if amino_site == 215:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_M1_determinants-of-virulence_30/Influenza A_M1_RNP-binding-region_165")
                 Comment.append("These two M1 residues contribute to differences in virulence of H5N1 avian influenza viruses in mice. Positions 30 and 215 are found in the N-terminal and C-terminal domains of M1 respectively./Introduction of Asn30Asp and Thr215Ala substitutions in the A/duck/Guangxi/53/2002 backbone conferred increased virulence in mice indicated by survival rate./This C-terminal region binds to vRNP.")
                 Reference.append("19117585/11222100/Uniprot:P03485")
                 
              elif amino_site == 207 or amino_site == 209:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_M1_Determinants-of-budding_30/Influenza A_M1_Determinants-of-virus-production_30/Influenza A_M1_determinant-of-virulence_30/Influenza A_M1_RNP-binding-region_165")
                 Comment.append("Key to spherical morphology. Substitution of any one of these residues with the corresponding TR swine residues causes transition to a filamentous morphology./Making M1 protein sufficient for virus-like particle production; indicating a unique feature of the pH1N1 M1 protein to induce and complete budding at the plasma membrane by itself./Single replacement of either one of the residue in M1 protein reduces overall viral production as well as growth kinetics./This C-terminal region binds to vRNP.")
                 Reference.append("23209789/11222100/Uniprot:P03485")
                 
              else:
                 if "*" in Amino_mutation :
                     Functional_relevant.append("Yes")
                     Feature.append("NA")
                     Comment.append("NA")
                     Reference.append("NA")
                  
                 
                 else:   
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_M1_RNP-binding-region_165")
                     Comment.append("This C-terminal region binds to vRNP.")
                     Reference.append("11222100/Uniprot:P03485")
                 
           else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA")  
        
        if str(Amino_mutation_S[0]) != str(Amino_mutation_S[-1]):
           
           if amino_site_S > 45 and amino_site_S < 61 and amino_site_S != 50:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_CRAC-motif_46")
              Comment_S.append("This is the 'cholesterol recognition consensus (CRAC) motif' found downstream of the transmembrane domain in the cytoplasmic tail region and also possess the C50 palmitoylation site. It is involved in M2 cholesterol binding.")
              Reference_S.append("15221235")
              
           elif amino_site_S == 11 or amino_site_S == 14 or amino_site_S == 16 or amino_site_S == 18 or amino_site_S == 20:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_determinant-of-host-range-specificity_11")
              Comment_S.append("At positions 11; 14; 16; 18 and 20; Weybridge possesses Thr; Gly; Glu; Ser and Ser; whereas WSN possesses Ile; Glu; Gly; Arg and Asn; respectively.")
              Reference_S.append("15777646")
              
           elif amino_site_S == 26:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_determinant-of-amantadine-resistance_26/Influenza A_M2_antiviral-activity_26")
              Comment_S.append("This residue is one of the critical amino acid occuring within the transmembrane domain of M2 protein. Substitution at this residue results in loss of sensitivity to M2 inhibitor drugs./Introduction of Leu26Phe substitution in the WSN/1933 backbone conferred decreased sensitivity to amantadine as indicated by plaque reduction assay in MDBK cells and reduced plaque sizes.")
              Reference_S.append("15673732/20834097")
              
           elif amino_site_S == 27:   
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_determinant-of-amantadine-resistance_27/Influenza A_M2_antiviral-activity_27") 
              Comment_S.append("It is shown that single-amino-acid substitution at this position within the transmembrane domain of M2 produces amantadine resistance in Influenza viruses/Introduction of Val27Ala substitution in the A/Chicken/Hong Kong/YU250/03 backbone conferred decreased sensitivity to amantadine.")
              Reference_S.append("15673732/20834097/16703504/16081121")
              
           elif amino_site_S == 30:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_determinant-of-amantadine-resistance_30/Influenza A_M2_antiviral-binding-site")
              Comment_S.append("It is shown that single-amino-acid substitution at this position within the transmembrane domain of M2 produces amantadine resistance in Influenza viruses/Rimantadine Binding Site")
              Reference_S.append("15673732/22078564")
              
           elif amino_site_S == 31:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_determinant-of-amantadine-resistance_31/Influenza A_M2_Antiviral-activity_31/Influenza A_M2_antiviral-binding-site_2LJC_30")
              Comment_S.append("It is shown that single-amino-acid substitution at this position within the transmembrane domain of M2 produces amantadine resistance in Influenza viruses/A/Chicken/Hebei/108/2002 isolate with Ser31Asn substitution conferred decreased sensitivity to amantadine and rimantadine in MDCK cells using plaque assay./Rimantadine Binding Site")
              Reference_S.append("15673732/20834097/16703504/2723453/17897729/17431677/15659762/17494553/16081121/22078564")
              
           elif amino_site_S == 34:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_determinant-of-amantadine-resistance_34/Influenza A_M2_antiviral-activity_34")
              Comment_S.append("It is shown that single-amino-acid substitution at this position within the transmembrane domain of M2 produces amantadine resistance in Influenza viruses/Introduction of Gly34Glu substitution in the WSN/1933 backbone conferred decreased sensitivity to amantadine as indicated by plaque reduction assay.")
              Reference_S.append("15673732")
              
           elif amino_site_S == 50:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_determinant-of-virulence_50/Influenza A_M2_CRAC-motif_46")
              Comment_S.append("The viruses lacking the palmitoylation site at this residue have been shown to cause a modest reduction in virulence in vivo (mouse models) although the effect is not seen tissue culture cells./This is the 'cholesterol recognition consensus (CRAC) motif' found downstream of the transmembrane domain in the cytoplasmic tail region and also possess the C50 palmitoylation site. It is involved in M2 cholesterol binding.")
              Reference_S.append("19553312/15221235")
              
           elif amino_site_S > 73 and amino_site_S < 80:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_determinants-of-infectivity_74")
              Comment_S.append("The amino acids 74 to 79 of the M2 tail play a role in virion morphogenesis and affect viral infectivity.")
              Reference_S.append("16699003")
              
           elif amino_site_S == 40 or amino_site_S == 41 or amino_site_S == 43 or amino_site_S == 45:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_M2_antiviral-binding-site")
              Comment_S.append("Rimantadine Binding Site")
              Reference_S.append("18235503")
              
           else:
               if "*" in Amino_mutation_S :
                  Functional_relevant_S.append("Yes")
                  Feature_S.append("NA")
                  Comment_S.append("NA")
                  Reference_S.append("NA")
                  
               else:
                  Functional_relevant_S.append("No")
                  Feature_S.append("NA")
                  Comment_S.append("NA")
                  Reference_S.append("NA")
                  
        if str(Amino_mutation_S[0]) == str(Amino_mutation_S[-1]):
           Functional_relevant_S.append("No")
           Feature_S.append("NA")
           Comment_S.append("NA")
           Reference_S.append("NA")     

if "NS" in name:               
   Functional_relevant=[]
   Feature=[]
   Comment=[]
   Reference=[]             
   Functional_relevant_S=[]
   Feature_S=[]
   Comment_S=[]
   Reference_S=[]
   for order in range(0,len(Annotate),1):
        Amino_mutation=Mutation[order]
        amino_site=int(Amino_mutation[1:-1])
        Amino_mutation_S=Mutation_S[order]
        amino_site_S=int(Amino_mutation_S[1:-1]) 
        
        if str(Amino_mutation[0]) != str(Amino_mutation[-1]):
           
           if amino_site > 0 and amino_site < 74:
               
              if amino_site == 35:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_NS1_RNA-binding-site_5/Influenza A_NS1_RBD-dimer-interface_12/Influenza A_NS1_nuclear-localization-signal-1_35/Influenza A_NS1_RNA-binding-domain_1")
                 Comment.append("Positions 35; 38; and 41 have been identified as critical amino acids regulating the functionality of NLS1 and are required for importin-alpha binding./Although many residues appear to contribute towards efficient RNA binding by NS1; only R38 is absolutely essential. Dimerization of the RBD is also essential for RNA binding/RBD dimerization is essential for RNA-binding activity/The N-terminal domain of NS1 binds several RNA species; including dsRNA. This domain also mediates interactions with RIG-I; possibly via dsRNA intermediates; PABPI; and importin-alpha.")
                 Reference.append("10024172/18813227/14967035/9360602/2969057/17376915/18796704/17475623/18813227")

              elif amino_site == 38 or amino_site == 41:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_NS1_nuclear-localization-signal-1_35/Influenza A_NS1_RNA-binding-site_5/Influenza A_NS1_RNA-binding-domain_1")
                 Comment.append("Positions 35; 38; and 41 have been identified as critical amino acids regulating the functionality of NLS1 and are required for importin-alpha binding./Although many residues appear to contribute towards efficient RNA binding by NS1; only R38 is absolutely essential. Dimerization of the RBD is also essential for RNA binding/The N-terminal domain of NS1 binds several RNA species; including dsRNA. This domain also mediates interactions with RIG-I; possibly via dsRNA intermediates; PABPI; and importin-alpha.")
                 Reference.append("2969057/17376915/10024172/18813227/14967035/18796704/17475623/18813227")

              elif amino_site == 29 or amino_site == 46:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_NS1_RNA-binding-site_5/Influenza A_NS1_RBD-dimer-interface_12/Influenza A_NS1_RNA-binding-domain_1")
                 Comment.append("Although many residues appear to contribute towards efficient RNA binding by NS1; only R38 is absolutely essential. Dimerization of the RBD is also essential for RNA binding/RBD dimerization is essential for RNA-binding activity/The N-terminal domain of NS1 binds several RNA species; including dsRNA. This domain also mediates interactions with RIG-I; possibly via dsRNA intermediates; PABPI; and importin-alpha.")
                 Reference.append("10024172/18813227/14967035/10024172/9360602/18796704/17475623/18813227")
              
                
              elif amino_site == 42:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_NS1_determinant-of-virulence_42/Influenza A_NS1_tissue-tropism_42/Influenza A_NS1_RNA-binding-site_5/Influenza A_NS1_RNA-binding-domain_1")
                 Comment.append("Introduction of Pro42Ser substitution from A/Duck/Guangxi/27/03 in the A/Duck/Guangxi/12/03 backbone conferred increased virulence as indicated by lethality in mice and the systemic spread of infection. This substitution also affects IFN pathway. Human epithelial lung A549 cells were infected with mutant A/Duck/Guangxi/12/03. Then supernatants from A549 cells were used to determine the levels of secreted IFN alpha/beta in bioassay. Infected cells did not inhibit viral replication./Although many residues appear to contribute towards efficient RNA binding by NS1; only R38 is absolutely essential. Dimerization of the RBD is also essential for RNA binding/The N-terminal domain of NS1 binds several RNA species; including dsRNA. This domain also mediates interactions with RIG-I; possibly via dsRNA intermediates; PABPI; and importin-alpha.")
                 Reference.append("18032512/10024172/1881322714967035/18796704/17475623/18813227") 
                 
              elif amino_site == 5 or amino_site == 31 or amino_site == 34 or amino_site == 37 or amino_site == 44 or amino_site == 45 or amino_site == 49:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_NS1_RNA-binding-site_5/Influenza A_NS1_RNA-binding-domain_1")
                 Comment.append("Although many residues appear to contribute towards efficient RNA binding by NS1; only R38 is absolutely essential. Dimerization of the RBD is also essential for RNA binding/The N-terminal domain of NS1 binds several RNA species; including dsRNA. This domain also mediates interactions with RIG-I; possibly via dsRNA intermediates; PABPI; and importin-alpha.")
                 Reference.append("10024172/1881322714967035/18796704/17475623/18813227")
                 
              elif amino_site == 12 or amino_site == 19 or amino_site == 32:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_NS1_RBD-dimer-interface_12/Influenza A_NS1_RNA-binding-domain_1")
                 Comment.append("RBD dimerization is essential for RNA-binding activity/The N-terminal domain of NS1 binds several RNA species; including dsRNA. This domain also mediates interactions with RIG-I; possibly via dsRNA intermediates; PABPI; and importin-alpha.")
                 Reference.append("10024172/9360602/18796704/17475623/18813227")

              else:
                 if "*" in Amino_mutation :
                     Functional_relevant.append("Yes")
                     Feature.append("NA")
                     Comment.append("NA")
                     Reference.append("NA")
                       
                 else:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_NS1_RNA-binding-domain_1")
                     Comment.append("The N-terminal domain of NS1 binds several RNA species; including dsRNA. This domain also mediates interactions with RIG-I; possibly via dsRNA intermediates; PABPI; and importin-alpha.")
                     Reference.append("18796704/17475623/18813227")
                  
           elif amino_site > 73 and amino_site < 87:
               
              if amino_site == 80:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_NS1_determinant-of-virulence_80/Influenza A_NS1_inter-domain-linker_74")
                 Comment.append("Flexible linker between the RNA binding and effector domains. Can be variable in length; with a 5aa deletion commonly reported in recent H5N1 isolates./Introduction of an artificial 15nt deletion in the recombinant virus A/WSN/33/(H1N1) x A/Duck/Shangdong/093/2004(H5N1) (HA/NA) conferred increased the virulence in mice using lethal dose and using intravenous pathogenicity index in chickens.")
                 Reference.append("18987632/Uniprot:P0349618317917/12195436")
                 
              elif amino_site > 80 and amino_site < 85:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_determinant-of-virulence_80/Influenza A_NS1_inter-domain-linker_74")
                 Comment.append("Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./Flexible linker between the RNA binding and effector domains. Can be variable in length; with a 5aa deletion commonly reported in recent H5N1 isolates./Introduction of an artificial 15nt deletion in the recombinant virus A/WSN/33/(H1N1) x A/Duck/Shangdong/093/2004(H5N1) (HA/NA) conferred increased the virulence in mice using lethal dose and using intravenous pathogenicity index in chickens.")
                 Reference.append("10938102/18987632/Uniprot:P0349618317917/12195436") 
               
              elif amino_site == 85 or amino_site == 86:
                 Functional_relevant.append("Yes")
                 Feature.append("Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_inter-domain-linker_74")
                 Comment.append("Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./Flexible linker between the RNA binding and effector domains. Can be variable in length; with a 5aa deletion commonly reported in recent H5N1 isolates.")
                 Reference.append("10938102/18987632/Uniprot:P03496")
                 
              else:
                 if "*" in Amino_mutation :
                     Functional_relevant.append("Yes")
                     Feature.append("NA")
                     Comment.append("NA")
                     Reference.append("NA")
                       
                 else:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_NS1_inter-domain-linker_74")
                     Comment.append("Flexible linker between the RNA binding and effector domains. Can be variable in length; with a 5aa deletion commonly reported in recent H5N1 isolates.")
                     Reference.append("18987632/Uniprot:P03496")
                     
                     
           elif amino_site > 86 and amino_site < 204:
               
              if amino_site > 86 and amino_site < 114:
                  
                 if amino_site > 89 and amino_site < 100:
                     
                    if amino_site > 89 and amino_site < 95:
                        
                       if amino_site == 91:
                          Functional_relevant.append("Yes")
                          Feature.append("Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_determinants-of-virulence-and-pathogenicity_90/Influenza A_NS1_eIF4GI-interaction-site_81")
                          Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./These 5 residues are present within the eIF4GI binding domain region of NS1 and their deletion increases the virulence and pathogenicity of H5N1 viruses./The changed length of the NS1 eIF4GI binding domain in H5N1 viruses with a 5-amino acid deletion can cause increased virulence and pathogenicity. Truncation of the eIF4GI binding domain attenuates replication invitro and invivo./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs.")
                          Reference.append("20133840/16963558/17881440/19403603/17229704/20854176/10938102/18725644/18796704/16715094/18585749")
                          
                       elif amino_site == 92:
                          Functional_relevant.append("Yes")
                          Feature.append("Influenza A_NS1_Replication-efficiency_92/Influenza A_NS1_Affect-type-I-IFN-pathway_92/Influenza A_NS1_determinants-of-virulence-and-pathogenicity_90/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                          Comment.append("Introduction of Glu92Asp in the A/HK/156/97 backbone conferred cytokine resistance using antiviral activity assay by comparing viral titers after pretreatment with IFN gamma; IFN alpha; TNF alpha. Introduction of Glu92Asp in the A/HK/156/97 backbone had viral titers similar to PR8 when inoculated pigs./These 5 residues are present within the eIF4GI binding domain region of NS1 and their deletion increases the virulence and pathogenicity of H5N1 viruses./The changed length of the NS1 eIF4GI binding domain in H5N1 viruses with a 5-amino acid deletion can cause increased virulence and pathogenicity. Truncation of the eIF4GI binding domain attenuates replication invitro and invivo./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                          Reference.append("12195436/18725644/18796704/16715094/18585749/10938102/20854176")
                        
                       else:   
                          if "*" in Amino_mutation :
                             Functional_relevant.append("Yes")
                             Feature.append("NA")
                             Comment.append("NA")
                             Reference.append("NA")
                       
                          else:
                             Functional_relevant.append("Yes")
                             Feature.append("Influenza A_NS1_determinants-of-virulence-and-pathogenicity_90/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                             Comment.append("These 5 residues are present within the eIF4GI binding domain region of NS1 and their deletion increases the virulence and pathogenicity of H5N1 viruses./The changed length of the NS1 eIF4GI binding domain in H5N1 viruses with a 5-amino acid deletion can cause increased virulence and pathogenicity. Truncation of the eIF4GI binding domain attenuates replication invitro and invivo./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                             Reference.append("18725644/18796704/16715094/18585749/10938102/20854176")
                 
                    elif amino_site == 96:
                         Functional_relevant.append("Yes")
                         Feature.append("Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_TRIM25-interaction-site_96/Influenza A_NS1_determinants-of-virulence_90/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                         Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./Necessary for interaction with TRIM25; a mechanism by which NS1 inhibits RIG-I. Same residues are also essential for PI3K activation; possibly by making direct contact with the p110 activation loop /The changed length of the NS1 eIF4GI binding domain in H5N1 viruses with a 5-amino acid deletion can cause increased virulence and pathogenicity. Truncation of the eIF4GI binding domain attenuates replication invitro and invivo./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                         Reference.append("20133840/16963558/17881440/19403603/17229704/19454348/20854176/10938102/18725644/18796704/16715094/18585749")
                    
                    elif amino_site == 97:
                         Functional_relevant.append("Yes")
                         Feature.append("Influenza A_NS1_TRIM25-interaction-site_96/Influenza A_NS1_determinants-of-virulence_90/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                         Comment.append("Necessary for interaction with TRIM25; a mechanism by which NS1 inhibits RIG-I. Same residues are also essential for PI3K activation; possibly by making direct contact with the p110 activation loop /The changed length of the NS1 eIF4GI binding domain in H5N1 viruses with a 5-amino acid deletion can cause increased virulence and pathogenicity. Truncation of the eIF4GI binding domain attenuates replication invitro and invivo./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                         Reference.append("19454348/20854176/10938102/18725644/18796704/16715094/18585749")
                    
                    elif amino_site == 95 or amino_site == 98 or amino_site == 99:
                         Functional_relevant.append("Yes")
                         Feature.append("Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_determinants-of-virulence_90/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                         Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./The changed length of the NS1 eIF4GI binding domain in H5N1 viruses with a 5-amino acid deletion can cause increased virulence and pathogenicity. Truncation of the eIF4GI binding domain attenuates replication invitro and invivo./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                         Reference.append("20133840/16963558/17881440/19403603/17229704/20854176/10938102/18725644/18796704/16715094/18585749")
                        
                    else:
                        if "*" in Amino_mutation :
                           Functional_relevant.append("Yes")
                           Feature.append("NA")
                           Comment.append("NA")
                           Reference.append("NA")
                    
                        else:
                           Functional_relevant.append("Yes")
                           Feature.append("Influenza A_NS1_determinants-of-virulence_90/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                           Comment.append("The changed length of the NS1 eIF4GI binding domain in H5N1 viruses with a 5-amino acid deletion can cause increased virulence and pathogenicity. Truncation of the eIF4GI binding domain attenuates replication invitro and invivo./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                           Reference.append("20854176/10938102/18725644/18796704/16715094/18585749")

                 elif amino_site == 87 or amino_site == 89:
                      Functional_relevant.append("Yes")
                      Feature.append("Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                      Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                      Reference.append("20133840/16963558/17881440/19403603/17229704/10938102/18725644/18796704/16715094/18585749")

                 elif amino_site == 101:
                      Functional_relevant.append("Yes")
                      Feature.append("Influenza A_NS1_determinant-of-virulence_101/Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                      Comment.append("Passaging of HK156 viruses in mouse brain and embryonated eggs led to the selection of high and low virulent variants in the mice model respectively. These phenotypic changes are confered by changes in amino acids in HA residues (211); PB1 (456 and 712); PA (631); NP (127) and NS1 (101) proteins./Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                      Reference.append("10873787/20133840/16963558/17881440/19403603/17229704/10938102/18725644/18796704/16715094/18585749")

                 elif amino_site == 103:
                      Functional_relevant.append("Yes")
                      Feature.append("Influenza A_NS1_tissue-tropism and determine of virulence_103/Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                      Comment.append("Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./Introduction of Leu103Phe and Ile106Met substitutions in the A/Hong Kong/483/1997 backbone conferred increased virulence compared to WT by measuring lethality in mice. This dual substitution also spread systemically after measuring viral titers in lung; peripheral blood; spleen and brain. The histopathological assessment of lungs in mice show lung inflammation; accumulation of neutrophils and exudate in the alveolar spaces./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                      Reference.append("17522219/12667806/20444891/18725644/17442719/10938102/18725644/18796704/16715094/18585749")

                 elif amino_site == 106:
                      Functional_relevant.append("Yes")
                      Feature.append("Influenza A_NS1_ED-helix-dimer_106/Influenza A_NS1_tissue-tropism and determine of virulence_103/Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                      Comment.append("The helix-helix ED conformation is conserved in all apo-ED crystal structures; though its physiological relevance is unknown. W187 is essential for ED dimerization in vitro./Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./Introduction of Leu103Phe and Ile106Met substitutions in the A/Hong Kong/483/1997 backbone conferred increased virulence compared to WT by measuring lethality in mice. This dual substitution also spread systemically after measuring viral titers in lung; peripheral blood; spleen and brain. The histopathological assessment of lungs in mice show lung inflammation; accumulation of neutrophils and exudate in the alveolar spaces./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                      Reference.append("18585749/20133840/17522219/12667806/20444891/18725644/17442719/10938102/18725644/18796704/16715094/18585749")


                 elif amino_site == 108 or amino_site == 109 or amino_site == 110:
                      Functional_relevant.append("Yes")
                      Feature.append("Influenza A_NS1_ED-helix-dimer_106/Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                      Comment.append("The helix-helix ED conformation is conserved in all apo-ED crystal structures; though its physiological relevance is unknown. W187 is essential for ED dimerization in vitro./Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                      Reference.append("18585749/20133840/17522219/12667806/20444891/18725644/17442719/10938102/18725644/18796704/16715094/18585749")
 
                 elif amino_site == 105 or amino_site == 107:
                      Functional_relevant.append("Yes")
                      Feature.append("Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                      Comment.append("Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                      Reference.append("17522219/12667806/20444891/18725644/17442719/10938102/18725644/18796704/16715094/18585749")

                 else:         
                     if "*" in Amino_mutation :
                        Functional_relevant.append("Yes")
                        Feature.append("NA")
                        Comment.append("NA")
                        Reference.append("NA")
                    
                     else:
                        Functional_relevant.append("Yes")
                        Feature.append("Influenza A_NS1_eIF4GI-interaction-site_81/Influenza A_NS1_effector-domain_87")
                        Comment.append("Interaction of NS1 with eIF4GI may lead to preferential translation of viral mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                        Reference.append("10938102/18725644/18796704/16715094/18585749")
              
              elif amino_site == 117:
                   Functional_relevant.append("Yes")
                   Feature.append("Influenza A_NS1_ED-helix-dimer_106/Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_effector-domain_87")
                   Comment.append("The helix-helix ED conformation is conserved in all apo-ED crystal structures; though its physiological relevance is unknown. W187 is essential for ED dimerization in vitro./Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                   Reference.append("18585749/20133840/17522219/12667806/20444891/18725644/17442719/18725644/18796704/16715094/18585749")
              
              elif amino_site == 118:
                   Functional_relevant.append("Yes")
                   Feature.append("Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_effector-domain_87")
                   Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                   Reference.append("20133840/16963558/17881440/19403603/17229704/18725644/18796704/16715094/18585749")
                   
                   
              elif amino_site > 118 and amino_site < 127:
                  
                   if amino_site == 124:
                      Functional_relevant.append("Yes")
                      Feature.append("Influenza A_NS1_ED-helix-dimer_106/Influenza A_NS1_PKR-binding-site_123/Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_effector-domain_87")
                      Comment.append("Necessary for the interaction with PKR; resulting in an inhibition of eIF2alpha phosphorylation./The helix-helix ED conformation is conserved in all apo-ED crystal structures; though its physiological relevance is unknown. W187 is essential for ED dimerization in vitro./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain./Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs.")
                      Reference.append("17320139/18585749/20133840/18725644/18796704/16715094/1858574917522219/12667806/20444891/18725644/17442719")
                  
                   elif amino_site == 123 or amino_site == 126:
                        Functional_relevant.append("Yes")
                        Feature.append("Influenza A_NS1_PKR-binding-site_123/Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_effector-domain_87")
                        Comment.append("Necessary for the interaction with PKR; resulting in an inhibition of eIF2alpha phosphorylation./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain./Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs.")
                        Reference.append("17320139/18725644/18796704/16715094/1858574917522219/12667806/20444891/18725644/17442719")
                  
                   elif amino_site == 119 or amino_site == 121:
                        Functional_relevant.append("Yes")
                        Feature.append("Influenza A_NS1_ED-helix-dimer_106/Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_effector-domain_87")
                        Comment.append("The helix-helix ED conformation is conserved in all apo-ED crystal structures; though its physiological relevance is unknown. W187 is essential for ED dimerization in vitro./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain./Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs.")
                        Reference.append("18585749/20133840/18725644/18796704/16715094/1858574917522219/12667806/20444891/18725644/17442719")
                       
                   elif amino_site == 125:
                        Functional_relevant.append("Yes")
                        Feature.append("Influenza A_NS1_determinant-of-pathogenicity_125/Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_effector-domain_87")
                        Comment.append("The single mutation Asp125Gly causes high pathogenicity in mice at late passage 10 and also causes enhanced binding abilities to Alpha-2-3 and Alpha-2-6 sialic acid-linked receptors./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain./Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs.")
                        Reference.append("18983930/18725644/18796704/16715094/1858574917522219/12667806/20444891/18725644/17442719")
                                    
                   else:         
                       if "*" in Amino_mutation :
                          Functional_relevant.append("Yes")
                          Feature.append("NA")
                          Comment.append("NA")
                          Reference.append("NA")
                    
                       else:
                          Functional_relevant.append("Yes")
                          Feature.append("Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_effector-domain_87")
                          Comment.append("The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain./Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs.")
                          Reference.append("18725644/18796704/16715094/1858574917522219/12667806/20444891/18725644/17442719")

              elif amino_site == 127:
                   Functional_relevant.append("Yes")
                   Feature.append("Influenza A_NS1_PKR-binding-site_123/Influenza A_NS1_effector-domain_87")
                   Comment.append("Necessary for the interaction with PKR; resulting in an inhibition of eIF2alpha phosphorylation./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                   Reference.append("17320139/18725644/18796704/16715094/18585749")

              elif amino_site == 133 or amino_site == 135:
                   Functional_relevant.append("Yes")
                   Feature.append("Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_effector-domain_87")
                   Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                   Reference.append("20133840/16963558/17881440/19403603/17229704/18725644/18796704/16715094/18585749") 


              elif amino_site > 136 and amino_site < 148:
                   
                   if amino_site == 142 or amino_site == 143 or amino_site == 145 or amino_site == 146:
                      Functional_relevant.append("Yes")
                      Feature.append("Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_NES-mask_148/Influenza A_NS1_p85b-binding-site_87")
                      Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity/Sequence added to a heterologous protein causes nuclear export. L144 and L146 are essential for this activity.")
                      Reference.append("20133840/16963558/17881440/19403603/17229704/9560194/Uniprot:P03496/20133840/16963558/17881440/19403603/17229704") 
                  
                  
                   else:
                       if "*" in Amino_mutation :
                          Functional_relevant.append("Yes")
                          Feature.append("NA")
                          Comment.append("NA")
                          Reference.append("NA")
                    
                       else:
                          Functional_relevant.append("Yes")
                          Feature.append("Influenza A_NS1_NES-mask_148/Influenza A_NS1_p85b-binding-site_87")
                          Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity/Sequence added to a heterologous protein causes nuclear export. L144 and L146 are essential for this activity.")
                          Reference.append("9560194/Uniprot:P03496/20133840/16963558/17881440/19403603/17229704")

              elif amino_site > 147 and amino_site < 162:
                  
                   if amino_site == 148 or amino_site == 159 or amino_site == 161:
                      Functional_relevant.append("Yes") 
                      Feature.append("Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_NES-mask_148/Influenza A_NS1_effector-domain_87")
                      Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./Sequence added to Influenza A_NS1_SF29 inhibits NES activity. R148; E152; and E153 are critical for the function of the mask; including in the context of full-length NS1./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                      Reference.append("20133840/16963558/17881440/19403603/17229704/9560194/18725644/18796704/16715094/18585749") 
                  
                   elif amino_site == 151 or amino_site == 153 or amino_site == 155 or amino_site == 156 or amino_site == 157:
                      Functional_relevant.append("Yes") 
                      Feature.append("Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_NES-mask_148/Influenza A_NS1_effector-domain_87")
                      Comment.append("Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./Sequence added to Influenza A_NS1_SF29 inhibits NES activity. R148; E152; and E153 are critical for the function of the mask; including in the context of full-length NS1./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                      Reference.append("17522219/12667806/20444891/18725644/17442719/9560194/18725644/18796704/16715094/18585749") 
                                      
                   else:
                       if "*" in Amino_mutation :
                          Functional_relevant.append("Yes")
                          Feature.append("NA")
                          Comment.append("NA")
                          Reference.append("NA")

                       else:
                          Functional_relevant.append("Yes") 
                          Feature.append("Influenza A_NS1_NES-mask_148/Influenza A_NS1_effector-domain_87")
                          Comment.append("Sequence added to Influenza A_NS1_SF29 inhibits NES activity. R148; E152; and E153 are critical for the function of the mask; including in the context of full-length NS1./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                          Reference.append("9560194/18725644/18796704/16715094/18585749")
                   
              elif amino_site == 180 or amino_site == 181 or amino_site == 183 or amino_site == 184 or amino_site == 187 or amino_site == 188 or amino_site == 189:
                   Functional_relevant.append("Yes")
                   Feature.append("Influenza A_NS1_ED-helix-dimer_106/Influenza A_NS1_CPSF30-binding-site_103/Influenza A_NS1_effector-domain_87")
                   Comment.append("The helix-helix ED conformation is conserved in all apo-ED crystal structures; though its physiological relevance is unknown. W187 is essential for ED dimerization in vitro./Necessary for binding CPSF30 and inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                   Reference.append("18585749/20133840/17522219/12667806/20444891/18725644/17442719/18725644/18796704/16715094/18585749")

              elif amino_site == 162 or amino_site == 164:
                   Functional_relevant.append("Yes")
                   Feature.append("Influenza A_NS1_p85b-binding-site_87/Influenza A_NS1_effector-domain_87")
                   Comment.append("Necessary for binding p85beta and activating PI3K signaling. Required for efficient virus replication; perhaps by delaying host-cell apoptosis or by enhancing sodium channel activity./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                   Reference.append("20133840/16963558/17881440/19403603/17229704/18725644/18796704/16715094/18585749") 

              elif amino_site == 179 or amino_site == 186:
                   Functional_relevant.append("Yes")
                   Feature.append("Influenza A_NS1_ED-helix-dimer_106/Influenza A_NS1_effector-domain_87")
                   Comment.append("The helix-helix ED conformation is conserved in all apo-ED crystal structures; though its physiological relevance is unknown. W187 is essential for ED dimerization in vitro./The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                   Reference.append("18585749/20133840/18725644/18796704/16715094/18585749")

              else:
                  if "*" in Amino_mutation :
                     Functional_relevant.append("Yes")
                     Feature.append("NA")
                     Comment.append("NA")
                     Reference.append("NA")

                  else:
                     Functional_relevant.append("Yes")
                     Feature.append("Influenza A_NS1_effector-domain_87")
                     Comment.append("The NS1 effector domain mediates interactions with several host proteins and may stabilize the N-terminal RNA-binding domain.")
                     Reference.append("18725644/18796704/16715094/18585749") 

           elif amino_site > 203 and amino_site < 231:
               
                if amino_site > 212 and amino_site < 218:
                   if amino_site == 214 or amino_site == 215 or amino_site == 217: 
                      Functional_relevant.append("Yes") 
                      Feature.append("Influenza A_NS1_Crk/CrkL-SH3-binding-site_212/Influenza A_NS1_phosphorylation-site_213/Influenza A_NS1_flexible-tail_204")
                      Comment.append("Crk/CrL-SH3 binding motif commonly found in avian influenza A virus strains as well as in the human 1918 pandemic virus. Shown to be required for Crk/CrkL binding./Threonine-215 is phosphorylated by a subset of proline-directed kinases (e.g. CDKs/ERKs) acting via a defined motif. Required for efficient virus replication in tissue-culture./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                      Reference.append("18165234/19007960/18585749/Uniprot:P03496")
                
                   else:
                       if "*" in Amino_mutation :
                          Functional_relevant.append("Yes")
                          Feature.append("NA")
                          Comment.append("NA")
                          Reference.append("NA")

                       else:
                          Functional_relevant.append("Yes") 
                          Feature.append("Influenza A_NS1_phosphorylation-site_213/Influenza A_NS1_flexible-tail_204")
                          Comment.append("Threonine-215 is phosphorylated by a subset of proline-directed kinases (e.g. CDKs/ERKs) acting via a defined motif. Required for efficient virus replication in tissue-culture./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                          Reference.append("19007960/18585749/Uniprot:P03496")
                
                elif amino_site == 219:
                     Functional_relevant.append("Yes") 
                     Feature.append("Influenza A_NS1_nuclear-localization-signal 2_219/Influenza A_NS1_phosphorylation-site_213/Influenza A_NS1_flexible-tail_204")
                     Comment.append("These basic residues are essential for NLS2 function and are required for importin-alpha binding. The same residues have been shown to form a nucleolar localization signal in some strains./Threonine-215 is phosphorylated by a subset of proline-directed kinases (e.g. CDKs/ERKs) acting via a defined motif. Required for efficient virus replication in tissue-culture./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                     Reference.append("2969057/17376915/19007960/18585749/Uniprot:P03496")
                
                elif amino_site == 220:
                     Functional_relevant.append("Yes") 
                     Feature.append("Influenza A_NS1_nuclear-localization-signal 2_219/Influenza A_NS1_flexible-tail_204")
                     Comment.append("These basic residues are essential for NLS2 function and are required for importin-alpha binding. The same residues have been shown to form a nucleolar localization signal in some strains./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                     Reference.append("2969057/17376915/18585749/Uniprot:P03496")
                     
                elif amino_site == 212:
                     Functional_relevant.append("Yes") 
                     Feature.append("Influenza A_NS1_Crk/CrkL-SH3-binding-site_212/Influenza A_NS1_flexible-tail_204")
                     Comment.append("Crk/CrL-SH3 binding motif commonly found in avian influenza A virus strains as well as in the human 1918 pandemic virus. Shown to be required for Crk/CrkL binding./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                     Reference.append("18165234/18585749/Uniprot:P03496")
                
                elif amino_site == 205 or amino_site == 210:
                     Functional_relevant.append("Yes") 
                     Feature.append("Influenza A_NS1_determinants-of-pathogenicity_205/Influenza A_NS1_flexible-tail_204")
                     Comment.append("Residues at positions 200 and 205 of NS1 contribute to enhanced type I interferon (IFN) antagonistic activity. Togehter; amino acid differences at residue 134 of HA; at 200 and 205 of NS1; and positions 47 and 51 of NS2 cause difference in virulence between high and low pathogenic H5N1 viruses./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                     Reference.append("20862325/18585749/Uniprot:P03496")
                
                elif amino_site > 222 and amino_site < 231:
                     
                     if amino_site > 226 and amino_site < 231:
                        if amino_site == 229:
                           Functional_relevant.append("Yes") 
                           Feature.append("Influenza A_NS1_nuclear-localization-signal 2_219/Influenza A_NS1_Tissue-tropism and Clinical-symptoms-of-disease_226/Influenza A_NS1_PDZ-ligand-motif_227/Influenza A_NS1_PABPII-binding-site_223/Influenza A_NS1_flexible-tail_204")
                           Comment.append("These basic residues are essential for NLS2 function and are required for importin-alpha binding. The same residues have been shown to form a nucleolar localization signal in some strains./Binds PDZ-domain containing proteins. The C-terminal motif is found primarily in avian isolates (common variants include av=ESEV/EPEV/KSEV; hu=RSKV/RSEV). Note: RSKV/RSEV is not thought to be a PDZ-ligand motif./Introduction of the PL motif at the C terminal in the virus A/WSN/33 conferred significant weight loss compared to WT. The virus variant showed severe alveolitis and hemorrhage in lung tissue of mice./Refers to poly(A)-binding protein II (PABPII)- binding region of NS1. May be involved in inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                           Reference.append("2969057/17376915/16439620/20702615/18334632/20410267/20686040/18334632/10205180/11421366/18585749/Uniprot:P03496") 
                                                        
                        else:
                            if "*" in Amino_mutation :
                               Functional_relevant.append("Yes")
                               Feature.append("NA")
                               Comment.append("NA")
                               Reference.append("NA")

                            else:
                               Functional_relevant.append("Yes") 
                               Feature.append("Influenza A_NS1_Tissue-tropism and Clinical-symptoms-of-disease_226/Influenza A_NS1_PDZ-ligand-motif_227/Influenza A_NS1_PABPII-binding-site_223/Influenza A_NS1_flexible-tail_204")
                               Comment.append("Binds PDZ-domain containing proteins. The C-terminal motif is found primarily in avian isolates (common variants include av=ESEV/EPEV/KSEV; hu=RSKV/RSEV). Note: RSKV/RSEV is not thought to be a PDZ-ligand motif./Introduction of the PL motif at the C terminal in the virus A/WSN/33 conferred significant weight loss compared to WT. The virus variant showed severe alveolitis and hemorrhage in lung tissue of mice./Refers to poly(A)-binding protein II (PABPII)- binding region of NS1. May be involved in inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                               Reference.append("16439620/20702615/18334632/20410267/20686040/18334632/10205180/11421366/18585749/Uniprot:P03496")
                     
                     elif amino_site == 226:
                          Functional_relevant.append("Yes") 
                          Feature.append("Influenza A_NS1_Tissue-tropism and Clinical-symptoms-of-disease_226/Influenza A_NS1_PABPII-binding-site_223/Influenza A_NS1_flexible-tail_204")
                          Comment.append("Introduction of the PL motif at the C terminal in the virus A/WSN/33 conferred significant weight loss compared to WT. The virus variant showed severe alveolitis and hemorrhage in lung tissue of mice./Refers to poly(A)-binding protein II (PABPII)- binding region of NS1. May be involved in inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                          Reference.append("18334632/10205180/11421366/18585749/Uniprot:P03496")
                                             
                     elif amino_site == 224:
                          Functional_relevant.append("Yes") 
                          Feature.append("Influenza A_NS1_nuclear-localization-signal 2_219/Influenza A_NS1_PABPII-binding-site_223/Influenza A_NS1_flexible-tail_204")
                          Comment.append("These basic residues are essential for NLS2 function and are required for importin-alpha binding. The same residues have been shown to form a nucleolar localization signal in some strains./Refers to poly(A)-binding protein II (PABPII)- binding region of NS1. May be involved in inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                          Reference.append("2969057/17376915/10205180/11421366/18585749/Uniprot:P03496")
                                                  
                     else:
                         if "*" in Amino_mutation :
                            Functional_relevant.append("Yes")
                            Feature.append("NA")
                            Comment.append("NA")
                            Reference.append("NA")

                         else:
                            Functional_relevant.append("Yes") 
                            Feature.append("Influenza A_NS1_PABPII-binding-site_223/Influenza A_NS1_flexible-tail_204")
                            Comment.append("Refers to poly(A)-binding protein II (PABPII)- binding region of NS1. May be involved in inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs./The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                            Reference.append("10205180/11421366/18585749/Uniprot:P03496")
                    
                else:
                    if "*" in Amino_mutation :
                       Functional_relevant.append("Yes")
                       Feature.append("NA")
                       Comment.append("NA")
                       Reference.append("NA")

                    else:
                       Functional_relevant.append("Yes") 
                       Feature.append("Influenza A_NS1_flexible-tail_204")
                       Comment.append("The flexible tail appears to be unstructured and variable in length. It contains a number of motifs; including CDK/ERK phosphorylation; Crk/CrkL SH3 binding; PDZ ligand and NoLS/NLS2.")
                       Reference.append("18585749/Uniprot:P03496")
               
               
           elif amino_site > 230 and amino_site < 238:
                if amino_site == 231 or amino_site == 232:
                   Functional_relevant.append("Yes") 
                   Feature.append("Influenza A_NS1_nuclear-localization-signal 2_219/Influenza A_NS1_C-terminal-extension_231/Influenza A_NS1_PABPII-binding-site_223")
                   Comment.append("These basic residues are essential for NLS2 function and are required for importin-alpha binding. The same residues have been shown to form a nucleolar localization signal in some strains./Function unknown; but may contribute to NoLS/NLS2. Reported to have become prevalent among human influenza A viruses isolated between 1950 and 1987./Refers to poly(A)-binding protein II (PABPII)- binding region of NS1. May be involved in inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs.")
                   Reference.append(":2969057/17376915/17376915/10205180/11421366") 
                    
                    
                else:
                    if "*" in Amino_mutation :
                       Functional_relevant.append("Yes")
                       Feature.append("NA")
                       Comment.append("NA")
                       Reference.append("NA")

                    else:
                       Functional_relevant.append("Yes") 
                       Feature.append("Influenza A_NS1_C-terminal-extension_231/Influenza A_NS1_PABPII-binding-site_223")
                       Comment.append("Function unknown; but may contribute to NoLS/NLS2. Reported to have become prevalent among human influenza A viruses isolated between 1950 and 1987./Refers to poly(A)-binding protein II (PABPII)- binding region of NS1. May be involved in inhibiting the posttranscriptional 3'-end processing of cellular pre-mRNAs.")
                       Reference.append("17376915/10205180/11421366")
               
           else:
               if "*" in Amino_mutation :
                  Functional_relevant.append("Yes")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
               else:
                  Functional_relevant.append("No")
                  Feature.append("NA")
                  Comment.append("NA")
                  Reference.append("NA")
                  
        if str(Amino_mutation[0]) == str(Amino_mutation[-1]):
           Functional_relevant.append("No")
           Feature.append("NA")
           Comment.append("NA")
           Reference.append("NA") 

        if str(Amino_mutation_S[0]) != str(Amino_mutation_S[-1]):
           
           if amino_site_S > 11 and amino_site_S < 22:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_NS2_nuclear-export-signal-motif_12")
              Comment_S.append("NA")
              Reference_S.append("Uniprot:P03508")
              
           elif amino_site_S == 47 or amino_site_S == 51:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_NS2_determinants-of-virulence_47")
              Comment_S.append("Residues 47 and 51 of NS2 are associated with difference in virulence between high and low pathogenic H5N1 viruses in ferrets. Amino acid differences at residue 134 of HA; at 200 and 205 of NS1 also contribute to this phenotype.")
              Reference_S.append("20862325")
              
           elif amino_site_S > 84 and amino_site_S < 95:
              Functional_relevant_S.append("Yes")
              Feature_S.append("Influenza A_NS2_nuclear-export-signal motif_85")
              Comment_S.append("NA")
              Reference_S.append("Uniprot:P03508")
              
           else:
               if "*" in Amino_mutation_S :
                  Functional_relevant_S.append("Yes")
                  Feature_S.append("NA")
                  Comment_S.append("NA")
                  Reference_S.append("NA")
                  
               else:
                  Functional_relevant_S.append("No")
                  Feature_S.append("NA")
                  Comment_S.append("NA")
                  Reference_S.append("NA")
                  
        if str(Amino_mutation_S[0]) == str(Amino_mutation_S[-1]):
           Functional_relevant_S.append("No")
           Feature_S.append("NA")
           Comment_S.append("NA")
           Reference_S.append("NA")        
 
                        
for i in range(0,len(Annotate),1):
    print(str(Sample[i]),str(Position[i]),str(Ref[i]),str(Cons[i]),str(Reads1[i]),str(Reads2[i]),str(VarFreq[i]),str(Strands1[i]),str(Strands2[i]),str(Qual1[i]),str(Qual2[i]),str(Pvalue[i]),str(MapQual1[i]),str(MapQual2[i]),str(Reads1Plus[i]),str(Reads1Minus[i]),str(Reads2Plus[i]),str(Reads2Minus[i]),str(VarAllele[i]),str(Day[i]),str(Treatment[i]),str(Segment[i]),str(Annotate[i]),str(Annotate_S[i]),str(Mutation[i]),str(Mutation_S[i]),str(Functional_relevant[i]),str(Feature[i]),str(Comment[i]),str(Reference[i]),str(Functional_relevant_S[i]),str(Feature_S[i]),str(Comment_S[i]),str(Reference_S[i]),sep=',', end='\n', file= out_file)

out_file.close()


