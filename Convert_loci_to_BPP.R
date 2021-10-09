####################### VALE INSTITUTE OF TECHNOLOGY ##########################

########################## CONVERT .LOCI TO BPP ############################



### Script prepared by Jeronymo Dalapicolla ###


#### PRE-ANALYSIS #### 

##1. INPUTS FOR THIS TUTORIAL ----
#A. A .LOCI FILE, A OUTPUT FROM ipyrad OR stacks PIPELINE.

#B. THE FILE ".CSV" WITH INFORMATION ON SAMPLES ID, ON WHICH CLADE/POP/LOCATION EACH INDIVIDUAL IS LOCATED IN THE TREE, AND A PROVISIONAL NAME INCLUDING ID AND CLADE, SPLIT BY UNDERSCORES(_)


##2. GOALS FOR THIS STEP ----
#A. EDIT .LOCI FILES FROM GENOMICS DATA TO USE IN BP&P ANALYSES;
#### OPTIONAL OPTIONS:
#B. SELECT LOCI THAT HAVE LESS MISSING DATA;
#C. SELECT A SPECIFIC NUMBER OF LOCI;
#D. REDUCE THE NUMBER OF INDIVIDUALS BY POPULATION;



##3. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! ----
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H

#B. MOVE TO YOUR WORKING DIRECTORY ALL NECESSARY FILES.


##4. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT ----
rm(list=ls())


##5. INSTALL AND LOAD THE PACKAGES ----
#From CRAN R:
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('tidyr'))        install.packages("tidyr");             library('tidyr')



#### ANALYSIS ---- 

###1. Read the file as matrix. Regular alignment functions take long time do to this ----
matrix_alig = scan('proechimys_MPE1.loci', what = 'character', sep = '\n')
class(matrix_alig)



###2. Define the number of characters to names ----
matrix_alig[1]
strsplit(matrix_alig[1], "")
names_char = 26 #26 in my case



###3. Define break lines for loci and other objects necessary to other steps ----
break.lines = grep('//', matrix_alig)
index.lines = c(0, break.lines)
loci_number = length(break.lines)
count_loci = c(0:loci_number)
count_ind = c(1, break.lines)



###4. Create a dataframe with information about the loci ----
len_loci = as.data.frame(matrix(NA, loci_number, 6))
colnames(len_loci) = c("Loci", "BreakLine", "ReadLength", "SampleSize", "InitialLine", "FinalLine")
loop = 0

for(i in 1:length(break.lines)){
  loop = loop+1
  num.1 = index.lines[i] + 1
  num.2 = index.lines[i + 1] - 1
  temp.text = matrix_alig[c(num.1:num.2)]
  len_loci[loop,] = c(count_loci[i], break.lines[i], nchar(as.character(matrix_alig[i])) - names_char, length(temp.text), num.1, num.2)
}

head(len_loci)
tail(len_loci)

#In case you wish use all loci regardless missing data, run these lines:
loci_subset = length(break.lines)
index_loci = len_loci

loci = list()
loop = 0
for (j in 1:loci_subset){
  loop = loop+1
  initial_line = index_loci$InitialLine[j]
  final_line = index_loci$FinalLine[j]
  loci[[loop]] = matrix_alig[c(initial_line:final_line)]
}
#verify first and last locus
loci[[1]]
loci[[loci_subset]]



###5. OPTIONAL: Subset a number of loci with more individuals (less missing data) ----
loci_subset = 200 #choose number of loci here!
index_loci = len_loci[order(len_loci$SampleSize, decreasing = T),]
index_loci = index_loci[1:loci_subset,1:length(names(len_loci))]
head(index_loci)

loci = list()
loop = 0
for (j in 1:loci_subset){
  loop = loop+1
  initial_line = index_loci$InitialLine[j]
  final_line = index_loci$FinalLine[j]
  loci[[loop]] = matrix_alig[c(initial_line:final_line)]
  }
#verify first and last locus
loci[[1]]
loci[[loci_subset]]



###6. Load file with information on populations or clades, putting a code for this in front of IDs ----
index_replace = read.csv("replace_names.csv")
head(index_replace)

#ID: original ID from .loci file
#ORDER: Individual order from alignment that generate .loci 
#CLADES: Clades where the sample is found in the phylogenetic trees. Could be a code for locations or pops
#NAME: A name including the code for pop or clades + _ + number of individual by population. For example, population A01 has 5 individuals, so we have "A01_01" "A01_02" "A01_03" "A01_04" "A01_05". This columns is important if you wanna remove individuals in populations with large sample size.
#REPLACE: Similar to previous however you use the original ID after the underscores "_".



###7. OPTIONAL: Reduce to three individuals to population. You can edit this number in the loop ---- 
ind_renamed = list()

for (k in 1:loci_subset){
  renamed_ind = loci[[k]]
  for (l in 1:length(index_replace$ID)){
    renamed_ind = stringr::str_replace(renamed_ind, index_replace$ID[l], index_replace$NAME[l])
  }
  ind_renamed[[k]] = renamed_ind
} 
#verify first and last locus
ind_renamed[[1]]
ind_renamed[[loci_subset]]

#7. Keep only three individuals by population to reduced sample size
ind_reduced = list()
for (m in 1:loci_subset){
  input = ind_renamed[[m]]
  ref01 = input[input %in% grep(paste0("_01", collapse = "|"), input, value = T)]
  ref02 = input[input %in% grep(paste0("_02", collapse = "|"), input, value = T)]
  ref03 = input[input %in% grep(paste0("_03", collapse = "|"), input, value = T)] # you can copy this line to put 4, 5... individuals
  input = c(ref01, ref02, ref03) #change here if you keep more then 3 individuals
  string.ed = input[1]
  find.p = '[[:alnum:]_]+[[:space:]]+'
  sub.pa = ''
  temp.dna.seq = sub(find.p,sub.pa,string.ed)
  dna.seq = strsplit(temp.dna.seq, '')
  loci.length = length(dna.seq[[1]])
  header = paste(length(input), loci.length , collapse = " ")
  ind_reduced[[m]] = c(header, input)
}
#verify first and last locus
ind_reduced[[1]]
ind_reduced[[loci_subset]]



###8. Rename individuals by order in the tree files ----
head(index_replace)


#If you run step 7 to remove individuals:
final_renamed = list()
for (k in 1:loci_subset){
  temp = ind_reduced[[k]]
  for (l in 1:length(index_replace$ID)){
  temp = gsub(index_replace$NAME[l], paste0(index_replace$CLADES[l],"_",index_replace$ORDER[l]), temp)
  }
  final_renamed[[k]] = temp
}  
#verify first and last locus
final_renamed[[1]]
final_renamed [[loci_subset]]


#If you DIDN'T run step 7 to remove individuals:
ind_reduced = loci

final_renamed = list()
for (k in 1:loci_subset){
  temp = ind_reduced[[k]]
  for (l in 1:length(index_replace$ID)){
    temp = gsub(index_replace$ID[l], paste0(index_replace$CLADES[l],"_",index_replace$ORDER[l]), temp)
  }
  final_renamed[[k]] = temp
} 

#verify first and last locus
final_renamed[[1]]
final_renamed [[loci_subset]]



####9. Merge all loci into a single files replaces "_" and save it.
bpp_file = c()
for (p in 1:loci_subset){
  add_file = as.character(final_renamed[[p]])
  without_first = add_file[-1]
  without_first = paste('', without_first, sep="       ")
  without_first = stringr::str_replace(without_first, "_", "^")
  add_file = c(add_file[1],without_first)
  bpp_file = c(bpp_file, add_file)
}


write.table(bpp_file, "input_genomics_bpp.loci", append = FALSE, sep = "\n", quote = F, row.names = F, col.names = F)



#END;
