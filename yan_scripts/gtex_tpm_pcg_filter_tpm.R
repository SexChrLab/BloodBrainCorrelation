## Continuation of gtex_tpm_filter_pcg.R

library(tidyverse)
library(dplyr)

## Read files from all_samples
working_dir <- "/Users/Sloth/Downloads/SexChrLab/BloodBrain/OutputFile/matched_sad_tissue/protein_coding"
setwd(working_dir)

pcg_genes_all_samples_folder <- paste(working_dir, "/all_samples", sep = "")
pcg_genes_male_samples_folder <- paste(working_dir, "/male_samples", sep ="")
pcg_genes_female_samples_folder <- paste(working_dir, "/female_samples", sep= "")

## In all_samples folder
## Identifies all files with the pattern _tpm.tsv
pcg_genes_all_samples_tsv <- list.files(pcg_genes_all_samples_folder, pattern = ".tsv")

## Reads my_files into R environment
## reset working directory, because read_tsv only works in working directory
setwd(pcg_genes_all_samples_folder)
pcg_genes_all_matched_samples <- lapply(pcg_genes_all_samples_tsv, read_tsv)

## Name all_tsv with their respective blood/brain tissue type
names(pcg_genes_all_matched_samples) <- gsub(".tsv", "", pcg_genes_all_samples_tsv, fixed=TRUE)

filenames <- names(pcg_genes_all_matched_samples)

######
## Create directory
pcg_genes_all_samples_tpm_1 <- paste(pcg_genes_all_samples_folder,"/gene_list_1", sep="")
dir.create(pcg_genes_all_samples_tpm_1)

pcg_genes_all_samples_tpm_5 <- paste(pcg_genes_all_samples_folder,"/gene_list_5", sep="")
dir.create(pcg_genes_all_samples_tpm_5)

pcg_genes_all_samples_tpm_10 <- paste(pcg_genes_all_samples_folder,"/gene_list_10", sep="")
dir.create(pcg_genes_all_samples_tpm_10)

#####
## Filter for TPM 1
# in_brain_1 - shows logical output FALSE/TRUE, tells you which line shows the gene that match criteria
pcg_genes_total_in_brain_1 <- list()
pcg_genes_total_in_blood_1 <- list()

## genes with TPM >= 1 in individual brain tissue type across all samples
## continuation of in_tissue_1, this one shows the genes with their TPM values
pcg_genes_total_brain_1 <- list()

## genes that do not have TPM >= 1 across all samples in individual brain tissue type
pcg_genes_total_brain_1_not <- list()

## genes expressed in blood across all samples with TPM >= 1
pcg_genes_total_blood_1 <- list()

## genes that do not have TPM >= 1 across all samples in blood
pcg_genes_total_blood_1_not <- list()

pcg_genes_total_both_brain_1 <- list()
pcg_genes_total_both_blood_1 <- list()

## genes that meet the TPM threshold in brain only
## gene in brain meets TPM threshold, gene in blood does not meet TPM threshold
pcg_genes_total_onlybrain_brain_1 <- list()
pcg_genes_total_onlybrain_blood_1 <- list()

## genes that meet the TPM threshold in blood only
pcg_genes_total_onlyblood_brain_1 <- list()
pcg_genes_total_onlyblood_blood_1 <- list()

## genes that meet the TPM threshold in neither blood or brain
pcg_genes_total_neither_brain_1 <- list()
pcg_genes_total_neither_blood_1 <- list()

num_sample_types <- 14
tpm_values  = c(1,5,10) # does this work? how to use if conditions to change tpm values
for (i in 1:num_sample_types){
  ## blood is 1:14, brain is 15:28
  ## this happens because the output files are blood first, then brain
  pcg_genes_total_in_blood_1[[i]] <- (rowSums(pcg_genes_all_matched_samples[[i]]>=1) == length(pcg_genes_all_matched_samples[[i]]))
  pcg_genes_total_in_brain_1[[i]] <- (rowSums(pcg_genes_all_matched_samples[[i+14]]>=1) == length(pcg_genes_all_matched_samples[[i+14]]))
  
  pcg_genes_total_blood_1[[i]] <- pcg_genes_all_matched_samples[[i]] %>% filter(pcg_genes_total_in_blood_1[[i]])
  pcg_genes_total_brain_1[[i]] <- pcg_genes_all_matched_samples[[i+14]] %>% filter(pcg_genes_total_in_brain_1[[i]])
  
  pcg_genes_total_blood_1_not[[i]] <- pcg_genes_all_matched_samples[[i]] %>% filter(!(pcg_genes_total_in_blood_1[[i]]))
  pcg_genes_total_brain_1_not[[i]] <- pcg_genes_all_matched_samples[[i+14]] %>% filter(!(pcg_genes_total_in_brain_1[[i]]))
  
  ## Both blood and brain
  pcg_genes_total_both_brain_1[[i]] <- pcg_genes_total_brain_1[[i]] %>% filter(pcg_genes_total_brain_1[[i]][[1]] %in% pcg_genes_total_blood_1[[i]][[1]])
  pcg_genes_total_both_blood_1[[i]] <- pcg_genes_total_blood_1[[i]] %>% filter(pcg_genes_total_blood_1[[i]][[1]] %in% pcg_genes_total_brain_1[[i]][[1]])
  
  ## Brain only 
  #pcg_genes_total_onlybrain_brain_1[[i]] <- pcg_genes_total_brain_1[[i]] %>% filter(pcg_genes_total_brain_1[[i]][[1]] %in% pcg_genes_total_blood_1_not[[i]][[1]])
  #pcg_genes_total_onlybrain_blood_1[[i]] <- pcg_genes_total_blood_1_not[[i]] %>% filter(pcg_genes_total_blood_1_not[[i]][[1]] %in% pcg_genes_total_brain_1[[i]][[1]])
  
  ## Blood only
  #pcg_genes_total_onlyblood_brain_1[[i]] <- pcg_genes_total_brain_1_not[[i]] %>% filter(pcg_genes_total_brain_1_not[[i]][[1]] %in% pcg_genes_total_blood_1[[i]][[1]])
  #pcg_genes_total_onlyblood_blood_1[[i]] <- pcg_genes_total_blood_1[[i]] %>% filter(pcg_genes_total_blood_1[[i]][[1]] %in% pcg_genes_total_brain_1_not[[i]][[1]])
  
  ## Neither blood or brain
  #pcg_genes_total_neither_brain_1[[i]] <- pcg_genes_total_brain_1_not[[i]] %>% filter(pcg_genes_total_brain_1_not[[i]][[1]] %in% pcg_genes_total_blood_1_not[[i]][[1]])
  #pcg_genes_total_neither_blood_1[[i]] <- pcg_genes_total_blood_1_not[[i]] %>% filter(pcg_genes_total_blood_1_not[[i]][[1]] %in% pcg_genes_total_brain_1_not[[i]][[1]])
  
}

#####
## File write out for TPM >= 1 only
## Directory and file name can be modified to include other TPMs
pcg_both_tpm_1 <- paste(pcg_genes_all_samples_tpm_1,"/both_tpm_1",sep="")
dir.create(pcg_both_tpm_1)

## only brain
pcg_onlybrain_tpm_1 <- paste(pcg_genes_all_samples_tpm_1, "/onlybrain_tpm_1", sep="")
dir.create(pcg_onlybrain_tpm_1)

## only blood
pcg_onlyblood_tpm_1 <- paste(pcg_genes_all_samples_tpm_1, "/onlyblood_tpm_1", sep="")
dir.create(pcg_onlyblood_tpm_1)

## neither
pcg_neither_tpm_1 <- paste(pcg_genes_all_samples_tpm_1,"/neither_tpm_1", sep="")
dir.create(pcg_neither_tpm_1)


for (i in 1:num_sample_types){
  ## both
  pcg_both_blood_tsv_1 <- paste(pcg_both_tpm_1, "/both_blood_",filenames[[i]],"_1.tsv", sep="")
  write.table(pcg_genes_total_both_blood_1[[i]], pcg_both_blood_tsv_1, row.names = FALSE, sep ="\t")
  
  pcg_both_brain_tsv_1 <- paste(pcg_both_tpm_1, "/both_brain_",filenames[[i+14]],"_1.tsv", sep="")
  write.table(pcg_genes_total_both_brain_1[[i]], pcg_both_brain_tsv_1, row.names = FALSE, sep ="\t")
  
  ## brain only
  pcg_onlybrain_blood_tsv_1 <- paste(pcg_onlybrain_tpm_1, "/onlybrain_blood_", filenames[[i]], "_1.tsv", sep="")
  write.table(pcg_genes_total_onlybrain_blood_1[[i]], pcg_onlybrain_blood_tsv_1, row.names =FALSE, sep ="\t")
  
  pcg_onlybrain_brain_tsv_1 <- paste(pcg_onlybrain_tpm_1, "/onlybrain_brain_", filenames[[i+14]], "_1.tsv", sep="")
  write.table(pcg_genes_total_onlybrain_brain_1[[i]], pcg_onlybrain_brain_tsv_1, row.names =FALSE, sep ="\t")
  
  ## blood only
  pcg_onlyblood_blood_tsv_1 <- paste(pcg_onlyblood_tpm_1, "/onlyblood_blood_", filenames[[i]], "_1.tsv", sep="")
  write.table(pcg_genes_total_onlyblood_blood_1[[i]], pcg_onlyblood_blood_tsv_1, row.names =FALSE, sep ="\t")
  
  pcg_onlyblood_brain_tsv_1 <- paste(pcg_onlyblood_tpm_1, "/onlyblood_brain_", filenames[[i+14]], "_1.tsv", sep="")
  write.table(pcg_genes_total_onlyblood_brain_1[[i]], pcg_onlyblood_brain_tsv_1, row.names =FALSE, sep ="\t")
  
  ## neither
  pcg_neither_blood_tsv_1 <- paste(pcg_neither_tpm_1,"/neither_blood_", filenames[[i]],"_1.tsv",sep="")
  write.table(pcg_genes_total_neither_blood_1[[i]], pcg_neither_blood_tsv_1, row.names=FALSE, sep="\t")
  
  pcg_neither_brain_tsv_1 <- paste(pcg_neither_tpm_1,"/neither_brain_",filenames[[i+14]],"_1.tsv",sep="")
  write.table(pcg_genes_total_neither_brain_1[[i]], pcg_neither_brain_tsv_1, row.names=FALSE, sep="\t")
  
}


#####
## Filter for TPM 5
# in_brain_1 - shows logical output FALSE/TRUE, tells you which line shows the gene that match criteria
pcg_genes_total_in_brain_5 <- list()
pcg_genes_total_in_blood_5 <- list()

## genes with TPM >= 5 in individual brain tissue type across all samples
## continuation of in_tissue_5, this one shows the genes with their TPM values
pcg_genes_total_brain_5 <- list()

## genes that do not have TPM >= 5 across all samples in individual brain tissue type
pcg_genes_total_brain_5_not <- list()

## genes expressed in blood across all samples with TPM >= 5
pcg_genes_total_blood_5 <- list()

## genes that do not have TPM >= 5 across all samples in blood
pcg_genes_total_blood_5_not <- list()

pcg_genes_total_both_brain_5 <- list()
pcg_genes_total_both_blood_5 <- list()

## genes that meet the TPM threshold in brain only
## gene in brain meets TPM threshold, gene in blood does not meet TPM threshold
pcg_genes_total_onlybrain_brain_5 <- list()
pcg_genes_total_onlybrain_blood_5 <- list()

## genes that meet the TPM threshold in blood only
pcg_genes_total_onlyblood_brain_5 <- list()
pcg_genes_total_onlyblood_blood_5 <- list()

## genes that meet the TPM threshold in neither blood or brain
pcg_genes_total_neither_brain_5 <- list()
pcg_genes_total_neither_blood_5 <- list()

num_sample_types <- 14
tpm_values  = c(1,5,10) # does this work? how to use if conditions to change tpm values
for (i in 1:num_sample_types){
  ## blood is 1:14, brain is 15:28
  ## this happens because the output files are blood first, then brain
  pcg_genes_total_in_blood_5[[i]] <- (rowSums(pcg_genes_all_matched_samples[[i]]>=5) == length(pcg_genes_all_matched_samples[[i]]))
  pcg_genes_total_in_brain_5[[i]] <- (rowSums(pcg_genes_all_matched_samples[[i+14]]>=5) == length(pcg_genes_all_matched_samples[[i+14]]))
  
  pcg_genes_total_blood_5[[i]] <- pcg_genes_all_matched_samples[[i]] %>% filter(pcg_genes_total_in_blood_5[[i]])
  pcg_genes_total_brain_5[[i]] <- pcg_genes_all_matched_samples[[i+14]] %>% filter(pcg_genes_total_in_brain_5[[i]])
  
  pcg_genes_total_blood_5_not[[i]] <- pcg_genes_all_matched_samples[[i]] %>% filter(!(pcg_genes_total_in_blood_5[[i]]))
  pcg_genes_total_brain_5_not[[i]] <- pcg_genes_all_matched_samples[[i+14]] %>% filter(!(pcg_genes_total_in_brain_5[[i]]))
  
  ## Both blood and brain
  pcg_genes_total_both_brain_5[[i]] <- pcg_genes_total_brain_5[[i]] %>% filter(pcg_genes_total_brain_5[[i]][[1]] %in% pcg_genes_total_blood_5[[i]][[1]])
  pcg_genes_total_both_blood_5[[i]] <- pcg_genes_total_blood_5[[i]] %>% filter(pcg_genes_total_blood_5[[i]][[1]] %in% pcg_genes_total_brain_5[[i]][[1]])
  
  ## Brain only 
  #pcg_genes_total_onlybrain_brain_5[[i]] <- pcg_genes_total_brain_5[[i]] %>% filter(pcg_genes_total_brain_5[[i]][[1]] %in% pcg_genes_total_blood_5_not[[i]][[1]])
  #pcg_genes_total_onlybrain_blood_5[[i]] <- pcg_genes_total_blood_5_not[[i]] %>% filter(pcg_genes_total_blood_5_not[[i]][[1]] %in% pcg_genes_total_brain_5[[i]][[1]])
  
  ## Blood only
  #pcg_genes_total_onlyblood_brain_5[[i]] <- pcg_genes_total_brain_5_not[[i]] %>% filter(pcg_genes_total_brain_5_not[[i]][[1]] %in% pcg_genes_total_blood_5[[i]][[1]])
  #pcg_genes_total_onlyblood_blood_5[[i]] <- pcg_genes_total_blood_5[[i]] %>% filter(pcg_genes_total_blood_5[[i]][[1]] %in% pcg_genes_total_brain_5_not[[i]][[1]])
  
  ## Neither blood or brain
  #pcg_genes_total_neither_brain_5[[i]] <- pcg_genes_total_brain_5_not[[i]] %>% filter(pcg_genes_total_brain_5_not[[i]][[1]] %in% pcg_genes_total_blood_5_not[[i]][[1]])
  #pcg_genes_total_neither_blood_5[[i]] <- pcg_genes_total_blood_5_not[[i]] %>% filter(pcg_genes_total_blood_5_not[[i]][[1]] %in% pcg_genes_total_brain_5_not[[i]][[1]])
  
}

#####
## File write out for TPM >= 5 only
## Directory and file name can be modified to include other TPMs
pcg_both_tpm_5 <- paste(pcg_genes_all_samples_tpm_5,"/both_tpm_5",sep="")
dir.create(pcg_both_tpm_5)

## only brain
pcg_onlybrain_tpm_5 <- paste(pcg_genes_all_samples_tpm_5, "/onlybrain_tpm_5", sep="")
dir.create(pcg_onlybrain_tpm_5)

## only blood
pcg_onlyblood_tpm_5 <- paste(pcg_genes_all_samples_tpm_5, "/onlyblood_tpm_5", sep="")
dir.create(pcg_onlyblood_tpm_5)

## neither
pcg_neither_tpm_5 <- paste(pcg_genes_all_samples_tpm_5,"/neither_tpm_5", sep="")
dir.create(pcg_neither_tpm_5)


for (i in 1:num_sample_types){
  ## both
  pcg_both_blood_tsv_5 <- paste(pcg_both_tpm_5, "/both_blood_",filenames[[i]],"_5.tsv", sep="")
  write.table(pcg_genes_total_both_blood_5[[i]], pcg_both_blood_tsv_5, row.names = FALSE, sep ="\t")
  
  pcg_both_brain_tsv_5 <- paste(pcg_both_tpm_5, "/both_brain_",filenames[[i+14]],"_5.tsv", sep="")
  write.table(pcg_genes_total_both_brain_5[[i]], pcg_both_brain_tsv_5, row.names = FALSE, sep ="\t")
  
  ## brain only
  pcg_onlybrain_blood_tsv_5 <- paste(pcg_onlybrain_tpm_5, "/onlybrain_blood_", filenames[[i]], "_5.tsv", sep="")
  write.table(pcg_genes_total_onlybrain_blood_5[[i]], pcg_onlybrain_blood_tsv_5, row.names =FALSE, sep ="\t")
  
  pcg_onlybrain_brain_tsv_5 <- paste(pcg_onlybrain_tpm_5, "/onlybrain_brain_", filenames[[i+14]], "_5.tsv", sep="")
  write.table(pcg_genes_total_onlybrain_brain_5[[i]], pcg_onlybrain_brain_tsv_5, row.names =FALSE, sep ="\t")
  
  ## blood only
  pcg_onlyblood_blood_tsv_5 <- paste(pcg_onlyblood_tpm_5, "/onlyblood_blood_", filenames[[i]], "_5.tsv", sep="")
  write.table(pcg_genes_total_onlyblood_blood_5[[i]], pcg_onlyblood_blood_tsv_5, row.names =FALSE, sep ="\t")
  
  pcg_onlyblood_brain_tsv_5 <- paste(pcg_onlyblood_tpm_5, "/onlyblood_brain_", filenames[[i+14]], "_5.tsv", sep="")
  write.table(pcg_genes_total_onlyblood_brain_5[[i]], pcg_onlyblood_brain_tsv_5, row.names =FALSE, sep ="\t")
  
  ## neither
  pcg_neither_blood_tsv_5 <- paste(pcg_neither_tpm_5,"/neither_blood_", filenames[[i]],"_5.tsv",sep="")
  write.table(pcg_genes_total_neither_blood_5[[i]], pcg_neither_blood_tsv_5, row.names=FALSE, sep="\t")
  
  pcg_neither_brain_tsv_5 <- paste(pcg_neither_tpm_5,"/neither_brain_",filenames[[i+14]],"_5.tsv",sep="")
  write.table(pcg_genes_total_neither_brain_5[[i]], pcg_neither_brain_tsv_5, row.names=FALSE, sep="\t")
  
}



#####
## Filter for TPM 10
# in_brain_1 - shows logical output FALSE/TRUE, tells you which line shows the gene that match criteria
pcg_genes_total_in_brain_10 <- list()
pcg_genes_total_in_blood_10 <- list()

## genes with TPM >= 1 in individual brain tissue type across all samples
## continuation of in_tissue_1, this one shows the genes with their TPM values
pcg_genes_total_brain_10 <- list()

## genes that do not have TPM >= 1 across all samples in individual brain tissue type
pcg_genes_total_brain_10_not <- list()

## genes expressed in blood across all samples with TPM >= 1
pcg_genes_total_blood_10 <- list()

## genes that do not have TPM >= 1 across all samples in blood
pcg_genes_total_blood_10_not <- list()

pcg_genes_total_both_brain_10 <- list()
pcg_genes_total_both_blood_10 <- list()

## genes that meet the TPM threshold in brain only
## gene in brain meets TPM threshold, gene in blood does not meet TPM threshold
pcg_genes_total_onlybrain_brain_10 <- list()
pcg_genes_total_onlybrain_blood_10 <- list()

## genes that meet the TPM threshold in blood only
pcg_genes_total_onlyblood_brain_10 <- list()
pcg_genes_total_onlyblood_blood_10 <- list()

## genes that meet the TPM threshold in neither blood or brain
pcg_genes_total_neither_brain_10 <- list()
pcg_genes_total_neither_blood_10 <- list()

num_sample_types <- 14
tpm_values  = c(1,5,10) # does this work? how to use if conditions to change tpm values
for (i in 1:num_sample_types){
  ## blood is 1:14, brain is 15:28
  ## this happens because the output files are blood first, then brain
  pcg_genes_total_in_blood_10[[i]] <- (rowSums(pcg_genes_all_matched_samples[[i]]>=10) == length(pcg_genes_all_matched_samples[[i]]))
  pcg_genes_total_in_brain_10[[i]] <- (rowSums(pcg_genes_all_matched_samples[[i+14]]>=10) == length(pcg_genes_all_matched_samples[[i+14]]))
  
  pcg_genes_total_blood_10[[i]] <- pcg_genes_all_matched_samples[[i]] %>% filter(pcg_genes_total_in_blood_10[[i]])
  pcg_genes_total_brain_10[[i]] <- pcg_genes_all_matched_samples[[i+14]] %>% filter(pcg_genes_total_in_brain_10[[i]])
  
  pcg_genes_total_blood_10_not[[i]] <- pcg_genes_all_matched_samples[[i]] %>% filter(!(pcg_genes_total_in_blood_10[[i]]))
  pcg_genes_total_brain_10_not[[i]] <- pcg_genes_all_matched_samples[[i+14]] %>% filter(!(pcg_genes_total_in_brain_10[[i]]))
  
  ## Both blood and brain
  pcg_genes_total_both_brain_10[[i]] <- pcg_genes_total_brain_10[[i]] %>% filter(pcg_genes_total_brain_10[[i]][[1]] %in% pcg_genes_total_blood_10[[i]][[1]])
  pcg_genes_total_both_blood_10[[i]] <- pcg_genes_total_blood_10[[i]] %>% filter(pcg_genes_total_blood_10[[i]][[1]] %in% pcg_genes_total_brain_10[[i]][[1]])
  
  ## Brain only 
  #pcg_genes_total_onlybrain_brain_10[[i]] <- pcg_genes_total_brain_10[[i]] %>% filter(pcg_genes_total_brain_10[[i]][[1]] %in% pcg_genes_total_blood_10_not[[i]][[1]])
  #pcg_genes_total_onlybrain_blood_10[[i]] <- pcg_genes_total_blood_10_not[[i]] %>% filter(pcg_genes_total_blood_10_not[[i]][[1]] %in% pcg_genes_total_brain_10[[i]][[1]])
  
  ## Blood only
  #pcg_genes_total_onlyblood_brain_10[[i]] <- pcg_genes_total_brain_10_not[[i]] %>% filter(pcg_genes_total_brain_10_not[[i]][[1]] %in% pcg_genes_total_blood_10[[i]][[1]])
  #pcg_genes_total_onlyblood_blood_10[[i]] <- pcg_genes_total_blood_10[[i]] %>% filter(pcg_genes_total_blood_10[[i]][[1]] %in% pcg_genes_total_brain_10_not[[i]][[1]])
  
  ## Neither blood or brain
  #pcg_genes_total_neither_brain_10[[i]] <- pcg_genes_total_brain_10_not[[i]] %>% filter(pcg_genes_total_brain_10_not[[i]][[1]] %in% pcg_genes_total_blood_10_not[[i]][[1]])
  #pcg_genes_total_neither_blood_10[[i]] <- pcg_genes_total_blood_10_not[[i]] %>% filter(pcg_genes_total_blood_10_not[[i]][[1]] %in% pcg_genes_total_brain_10_not[[i]][[1]])
  
}

#####
## File write out for TPM >= 10 only
## Directory and file name can be modified to include other TPMs
pcg_both_tpm_10 <- paste(pcg_genes_all_samples_tpm_10,"/both_tpm_10",sep="")
dir.create(pcg_both_tpm_10)

## only brain
pcg_onlybrain_tpm_10 <- paste(pcg_genes_all_samples_tpm_10, "/onlybrain_tpm_10", sep="")
dir.create(pcg_onlybrain_tpm_10)

## only blood
pcg_onlyblood_tpm_10 <- paste(pcg_genes_all_samples_tpm_10, "/onlyblood_tpm_10", sep="")
dir.create(pcg_onlyblood_tpm_10)

## neither
pcg_neither_tpm_10 <- paste(pcg_genes_all_samples_tpm_10,"/neither_tpm_10", sep="")
dir.create(pcg_neither_tpm_10)


for (i in 1:num_sample_types){
  ## both
  pcg_both_blood_tsv_10 <- paste(pcg_both_tpm_10, "/both_blood_",filenames[[i]],"_10.tsv", sep="")
  write.table(pcg_genes_total_both_blood_10[[i]], pcg_both_blood_tsv_10, row.names = FALSE, sep ="\t")
  
  pcg_both_brain_tsv_10 <- paste(pcg_both_tpm_10, "/both_brain_",filenames[[i+14]],"_10.tsv", sep="")
  write.table(pcg_genes_total_both_brain_10[[i]], pcg_both_brain_tsv_10, row.names = FALSE, sep ="\t")
  
  ## brain only
  pcg_onlybrain_blood_tsv_10 <- paste(pcg_onlybrain_tpm_10, "/onlybrain_blood_", filenames[[i]], "_10.tsv", sep="")
  write.table(pcg_genes_total_onlybrain_blood_10[[i]], pcg_onlybrain_blood_tsv_10, row.names =FALSE, sep ="\t")
  
  pcg_onlybrain_brain_tsv_10 <- paste(pcg_onlybrain_tpm_10, "/onlybrain_brain_", filenames[[i+14]], "_10.tsv", sep="")
  write.table(pcg_genes_total_onlybrain_brain_10[[i]], pcg_onlybrain_brain_tsv_10, row.names =FALSE, sep ="\t")
  
  ## blood only
  pcg_onlyblood_blood_tsv_10 <- paste(pcg_onlyblood_tpm_10, "/onlyblood_blood_", filenames[[i]], "_10.tsv", sep="")
  write.table(pcg_genes_total_onlyblood_blood_10[[i]], pcg_onlyblood_blood_tsv_10, row.names =FALSE, sep ="\t")
  
  pcg_onlyblood_brain_tsv_10 <- paste(pcg_onlyblood_tpm_10, "/onlyblood_brain_", filenames[[i+14]], "_10.tsv", sep="")
  write.table(pcg_genes_total_onlyblood_brain_10[[i]], pcg_onlyblood_brain_tsv_10, row.names =FALSE, sep ="\t")
  
  ## neither
  pcg_neither_blood_tsv_10 <- paste(pcg_neither_tpm_10,"/neither_blood_", filenames[[i]],"_10.tsv",sep="")
  write.table(pcg_genes_total_neither_blood_10[[i]], pcg_neither_blood_tsv_10, row.names=FALSE, sep="\t")
  
  pcg_neither_brain_tsv_10 <- paste(pcg_neither_tpm_10,"/neither_brain_",filenames[[i+14]],"_10.tsv",sep="")
  write.table(pcg_genes_total_neither_brain_10[[i]], pcg_neither_brain_tsv_10, row.names=FALSE, sep="\t")
  
}

num_pcg_filter_tpm_1 <- data.frame(matrix(ncol=num_sample_types,nrow=1))
num_pcg_filter_tpm_5 <- data.frame(matrix(ncol=num_sample_types,nrow=1))
num_pcg_filter_tpm_10 <- data.frame(matrix(ncol=num_sample_types,nrow=1))


for (i in 1:num_sample_types){
  num_pcg_filter_tpm_1[[i]] <- length(pcg_genes_total_both_brain_1[[i]][[1]])
  num_pcg_filter_tpm_5[[i]] <- length(pcg_genes_total_both_brain_5[[i]][[1]])
  num_pcg_filter_tpm_10[[i]] <- length(pcg_genes_total_both_brain_10[[i]][[1]])
}
num_pcg_filter_tpm <-as.data.frame(t(rbind(num_pcg_filter_tpm_1, num_pcg_filter_tpm_5, num_pcg_filter_tpm_10)), row.names = filenames)
colnames(num_pcg_filter_tpm)=c("pcg_tpm_1", "pcg_tpm_5", "pcg_tpm_10")
write.table(num_pcg_filter_tpm, paste0(pcg_genes_all_samples_folder, "/num_pcg_filter_tpm.tsv"), sep="\t", col.names=NA)
