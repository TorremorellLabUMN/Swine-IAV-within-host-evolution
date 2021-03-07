# Evolutionary rate

This script is for calculating the evolutionary rates of synonymous, non-synonymous and stop-gained mutations for each BALF sample at indivdual protein and whole genome level by summing the frequencies of corresponding variants and divided by the number of available sites and infection days (seven days).

The number of available sites for synonymous, non-synonymous and stop-gained mutations are equal to the 72%, 25% and 3% of the full nucleotide length of corresponding gene segments.

To execute the script, the information of identified SNVs are imported from the "BALF_H1.csv" and "BALF_H3.csv". The whole genome evolutionary rates are calculated based on obtained sequences from the Illumina Nextseq sequencing, the available sites of missing gene sequences are not included in the calculations. The information of failed sequenced genes are imported from "failed_sequencing.csv". 
