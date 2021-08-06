# Permutation test

The permutation test is performed for shared amino acid changes between pigs from different rooms within each treatment group for H1N1 and H3N2 influenza viruses. The scripts are adapted from https://github.com/blab/h5n1-cambodia/blob/master/figures/figure-5b-shared-sites-permutation-test.ipynb which was in previously published research article (Moncla LH, Bedford T, Dussart P, Horm SV, Rith S, Buchy P, et al. Quantifying within-host diversity of H5N1 influenza viruses in humans and poultry in Cambodia. PLoS Pathog. 2020;16: e1008191). 

# Script execution requirement
To execute the script, the information of number of amino acid variants identified in each gene segment and treatment groups need to be imported from the documents in import_files folder

The permutation test are performed separately based on 100% and 70% (assume 30% mutations are lethal) available coding region for amino acids changes.
