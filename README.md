# Ipyrad to BPP in R

Script to edit .loci files from ipyrad pipeline to fit it as a input in BP&P (Flouri et al. 2018);

##1. INPUTS FOR THIS TUTORIAL ----
#A. A .LOCI FILE, A OUTPUT FROM ipyrad OR stacks PIPELINE.

#B. THE FILE ".CSV" WITH INFORMATION ON SAMPLES ID, ON WHICH CLADE/POP/LOCATION EACH INDIVIDUAL IS LOCATED IN THE TREE, AND A PROVISIONAL NAME INCLUDING ID AND CLADE, SPLIT BY UNDERSCORES(_)


##2. GOALS FOR THIS STEP ----
#A. EDIT .LOCI FILES FROM GENOMICS DATA TO USE IN BP&P ANALYSES;
#### OPTIONAL OPTIONS:
#B. SELECT LOCI THAT HAVE LESS MISSING DATA;
#C. SELECT A SPECIFIC NUMBER OF LOCI;
#D. REDUCE THE NUMBER OF INDIVIDUALS BY POPULATION;


Eaton D.A.R., Overcast I. (2020). ipyrad: Interactive assembly and analysis of RADseq datasets. Bioinformatics, 36(8): 2592–2594. doi:10.1093/bioinformatics/btz966

Flouri T., Jiao X., Rannala B., Yang Z. (2018) Species Tree Inference with BPP using Genomic Sequences and the Multispecies Coalescent. Molecular Biology and Evolution, 35(10):2585-2593. doi:10.1093/molbev/msy147
