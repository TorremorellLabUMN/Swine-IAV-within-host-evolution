# Evolutionary rate

This script is for calculating the evolutionary rates of synonymous, non-synonymous and stop-gained mutations for each BALF sample at indivdual protein and whole genome level based on the method described before (Xue KS, Bloom JD. Linking influenza virus evolution within and between human hosts. Virus Evol. 2020;6: veaa010). Briefly, the evolutionary rates are calculated by summing the frequencies of corresponding variants and divided by the number of available sites and infection days (seven days).

The number of available sites for synonymous, non-synonymous and stop-gained mutations are equal to the 72%, 25% and 3% of the full nucleotide length of corresponding gene segments.

To execute the script, the information of identified SNVs are imported from the "H1N1_SNV.csv" and "H3N2_SNV.csv". The whole genome evolutionary rates are calculated based on obtained sequences from the Illumina Nextseq sequencing, the available sites of missing gene sequences are not included in the calculations. The information of failed sequenced genes are imported from "failed_sequencing.csv". 

Note: As only H1N1 NS gene amplified for sample 4931 and no SNV identified in this gene segment, the evolutionary rates of H1N1 NS1 and NS2 were mannually added in the exported files (BALF_H1_evolution_rate.csv) and we didn't calculate and use the H1N1 whole genome evolutionary rates of sample 4931 for statistically analysis. 
