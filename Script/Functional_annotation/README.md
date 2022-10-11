# Functional annotation
This script is used for annotating the functional related non-synonymous variants based on the Sequence Feature Variant Types data from the Influenza Research Database (Noronha JM, Liu M, Squires RB, Pickett BE, Hale BG, Air GM, et al. Influenza virus sequence feature variant type analysis: evidence of a role for NS1 in influenza virus host range restriction. J Virol. 2012;86: 5857–5866). 
# Command line
python Functional_annotation_V2.py input_file.csv output_file.csv

The input files used for this script should be the identifed variant files which exported from the script "SNV_annotation.py". 

Each input file should only contain the variants from the same coding region.
