# SNV annotation

This python script is used to annotate the affection of identified variants on amino acid changes for reference sequences from H1N1 or H3N2 challenge viruses. 
Each variant will be annotated as synonymous, stop-gained or non-synonymous based on whether they will change the coding amino acid region of challenge viruses.   

# Command line
The python script can be executed by command line:python SNV_annotation.py input_file_name.vcf output_file_name.vcf

The input file is the corrected vcf file generated by VarScan which listed the information of all the report variants. 
