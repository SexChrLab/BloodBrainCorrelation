## Following outputs genes with TPM >= 1
## To allow outputs with other TPMs, change >=1 to >=5 to >=10
## The following includes all genes
## A different section includes protein coding genes

library(tidyverse)
library(dplyr)

working_dir <- "/Users/Sloth/Downloads/SexChrLab/BloodBrain/OutputFile/matched_sad_tissue/all_genes"
setwd(working_dir)

all_genes_all_samples_folder <- paste(working_dir, "/all_samples", sep = "")
all_genes_male_samples_folder <- paste(working_dir, "/male_samples", sep ="")
all_genes_female_samples_folder <- paste(working_dir, "/female_samples", sep= "")

#####
## In all_samples folder
## Identifies all files with the pattern _tpm.tsv
all_genes_all_samples_tsv <- list.files(all_genes_all_samples_folder, pattern = ".tsv")

## Reads my_files into R environment
## reset working directory, because read_tsv only works in working directory
setwd(all_genes_all_samples_folder)
all_genes_all_matched_samples <- lapply(all_genes_all_samples_tsv, read_tsv)

## Name all_tsv with their respective blood/brain tissue type
names(all_genes_all_matched_samples) <- gsub(".tsv", "", all_genes_all_samples_tsv, fixed=TRUE)

filenames <- names(all_genes_all_matched_samples)

#####
## Create directory for TPM 1,5,10
all_genes_all_samples_tpm_1 <- paste(all_genes_all_samples_folder,"/gene_list_1", sep="")
dir.create(all_genes_all_samples_tpm_1)

all_genes_all_samples_tpm_5<- paste(all_genes_all_samples_folder,"/gene_list_5", sep="")
dir.create(all_genes_all_samples_tpm_5)

all_genes_all_samples_tpm_10 <- paste(all_genes_all_samples_folder,"/gene_list_10", sep="")
dir.create(all_genes_all_samples_tpm_10)


#####
# in_tissue_1 - shows logical output FALSE/TRUE, tells you which line shows the gene that match criteria
all_genes_total_in_brain_1 <- list()
all_genes_total_in_blood_1 <- list()

## genes with TPM >= 1 in individual brain tissue type across all samples
## continuation of in_tissue_1, this one shows the genes with their TPM values
all_genes_total_brain_1 <- list()

## genes that do not have TPM >= 1 across all samples in individual brain tissue type
all_genes_total_brain_1_not <- list()

## genes expressed in blood across all samples with TPM >= 1
all_genes_total_blood_1 <- list()

## genes that do not have TPM >= 1 across all samples in blood
all_genes_total_blood_1_not <- list()

## genes that meet the TPM threshold in BOTH brain and blood
## both - present in both brain and blood,
## brain - brain TPM value
## blood - blood TPM value
all_genes_total_both_brain_1 <- list()
all_genes_total_both_blood_1 <- list()

## genes that meet the TPM threshold in brain only
## gene in brain meets TPM threshold, gene in blood does not meet TPM threshold
all_genes_total_onlybrain_brain_1 <- list()
all_genes_total_onlybrain_blood_1 <- list()

## genes that meet the TPM threshold in blood only
all_genes_total_onlyblood_brain_1 <- list()
all_genes_total_onlyblood_blood_1 <- list()

## genes that meet the TPM threshold in neither blood or brain
all_genes_total_neither_brain_1 <- list()
all_genes_total_neither_blood_1 <- list()

num_sample_types <- 14
for (i in 1:num_sample_types){
  #all_genes_total_in_blood_1[[i]] <- (rowSums(all_genes_all_matched_samples[[i]]>=1) == length(all_genes_all_matched_samples[[i]]))
  #all_genes_total_in_brain_1[[i]] <- (rowSums(all_genes_all_matched_samples[[i+14]]>=1) == length(all_genes_all_matched_samples[[i+14]]))

  #all_genes_total_blood_1[[i]] <- all_genes_all_matched_samples[[i]] %>% filter(all_genes_total_in_blood_1[[i]])
  #all_genes_total_brain_1[[i]] <- all_genes_all_matched_samples[[i+14]] %>% filter(all_genes_total_in_brain_1[[i]])
  
  #all_genes_total_blood_1_not[[i]] <- all_genes_all_matched_samples[[i]] %>% filter(!(all_genes_total_in_blood_1[[i]]))
  #all_genes_total_brain_1_not[[i]] <- all_genes_all_matched_samples[[i+14]] %>% filter(!(all_genes_total_in_brain_1[[i]]))
  
  ## Both blood and brain
  all_genes_total_both_brain_1[[i]] <- all_genes_total_brain_1[[i]] %>% filter(all_genes_total_brain_1[[i]][[1]] %in% all_genes_total_blood_1[[i]][[1]])
  all_genes_total_both_blood_1[[i]] <- all_genes_total_blood_1[[i]] %>% filter(all_genes_total_blood_1[[i]][[1]] %in% all_genes_total_brain_1[[i]][[1]])
  
  ## Brain only 
  all_genes_total_onlybrain_brain_1[[i]] <- all_genes_total_brain_1[[i]] %>% filter(all_genes_total_brain_1[[i]][[1]] %in% all_genes_total_blood_1_not[[i]][[1]])
  all_genes_total_onlybrain_blood_1[[i]] <- all_genes_total_blood_1_not[[i]] %>% filter(all_genes_total_blood_1_not[[i]][[1]] %in% all_genes_total_brain_1[[i]][[1]])
  
  ## Blood only
  all_genes_total_onlyblood_brain_1[[i]] <- all_genes_total_brain_1_not[[i]] %>% filter(all_genes_total_brain_1_not[[i]][[1]] %in% all_genes_total_blood_1[[i]][[1]])
  all_genes_total_onlyblood_blood_1[[i]] <- all_genes_total_blood_1[[i]] %>% filter(all_genes_total_blood_1[[i]][[1]] %in% all_genes_total_brain_1_not[[i]][[1]])
  
  ## Neither blood or brain
  all_genes_total_neither_brain_1[[i]] <- all_genes_total_brain_1_not[[i]] %>% filter(all_genes_total_brain_1_not[[i]][[1]] %in% all_genes_total_blood_1_not[[i]][[1]])
  all_genes_total_neither_blood_1[[i]] <- all_genes_total_blood_1_not[[i]] %>% filter(all_genes_total_blood_1_not[[i]][[1]] %in% all_genes_total_brain_1_not[[i]][[1]])
  
}

#####
## File write out for TPM >= 1 only
## Directory and file name can be modified to include other TPMs
both_tpm_1 <- paste(all_genes_all_samples_tpm_1,"/both_tpm_1",sep="")
dir.create(both_tpm_1)

## only brain
onlybrain_tpm_1 <- paste(all_genes_all_samples_tpm_1, "/onlybrain_tpm_1", sep="")
dir.create(onlybrain_tpm_1)

## only blood
onlyblood_tpm_1 <- paste(all_genes_all_samples_tpm_1, "/onlyblood_tpm_1", sep="")
dir.create(onlyblood_tpm_1)

## neither
neither_tpm_1 <- paste(all_genes_all_samples_tpm_1,"/neither_tpm_1", sep="")
dir.create(neither_tpm_1)

for (i in 1:num_sample_types){
  ## both
  both_blood_tsv_1 <- paste(both_tpm_1, "/both_blood_",filenames[[i]],"_1.tsv", sep="")
  write.table(all_genes_total_both_blood_1[[i]], both_blood_tsv_1, row.names = FALSE, sep ="\t")
  
  both_brain_tsv_1 <- paste(both_tpm_1, "/both_brain_",filenames[[i+14]],"_1.tsv", sep="")
  write.table(all_genes_total_both_brain_1[[i]], both_brain_tsv_1, row.names = FALSE, sep ="\t")
 
  ## brain only
  onlybrain_blood_tsv_1 <- paste(onlybrain_tpm_1, "/onlybrain_blood_", filenames[[i]], "_1.tsv", sep="")
  write.table(all_genes_total_onlybrain_blood_1[[i]], onlybrain_blood_tsv_1, row.names =FALSE, sep ="\t")
  
  onlybrain_brain_tsv_1 <- paste(onlybrain_tpm_1, "/onlybrain_brain_", filenames[[i+14]], "_1.tsv", sep="")
  write.table(all_genes_total_onlybrain_brain_1[[i]], onlybrain_brain_tsv_1, row.names =FALSE, sep ="\t")
  
  ## blood only
  onlyblood_blood_tsv_1 <- paste(onlyblood_tpm_1, "/onlyblood_blood_", filenames[[i]], "_1.tsv", sep="")
  write.table(all_genes_total_onlyblood_blood_1[[i]], onlyblood_blood_tsv_1, row.names =FALSE, sep ="\t")
  
  onlyblood_brain_tsv_1 <- paste(onlyblood_tpm_1, "/onlyblood_brain_", filenames[[i+14]], "_1.tsv", sep="")
  write.table(all_genes_total_onlyblood_brain_1[[i]], onlyblood_brain_tsv_1, row.names =FALSE, sep ="\t")

  ## neither
  neither_blood_tsv_1 <- paste(neither_tpm_1,"/neither_blood_", filenames[[i]],"_1.tsv",sep="")
  write.table(all_genes_total_neither_blood_1[[i]], neither_blood_tsv_1, row.names=FALSE, sep="\t")
  
  neither_brain_tsv_1 <- paste(neither_tpm_1,"/neither_brain_",filenames[[i+14]],"_1.tsv",sep="")
  write.table(all_genes_total_neither_brain_1[[i]], neither_brain_tsv_1, row.names=FALSE, sep="\t")
  
}


#####
## Filter for TPM 5

# in_tissue_1 - shows logical output FALSE/TRUE, tells you which line shows the gene that match criteria
all_genes_total_in_brain_5 <- list()
all_genes_total_in_blood_5 <- list()

## genes with TPM >= 5 in individual brain tissue type across all samples
## continuation of in_tissue_5, this one shows the genes with their TPM values
all_genes_total_brain_5 <- list()

## genes that do not have TPM >= 5 across all samples in individual brain tissue type
all_genes_total_brain_5_not <- list()

## genes expressed in blood across all samples with TPM >= 5
all_genes_total_blood_5 <- list()

## genes that do not have TPM >= 5 across all samples in blood
all_genes_total_blood_5_not <- list()

## genes that meet the TPM threshold in BOTH brain and blood
## both - present in both brain and blood,
## brain - brain TPM value
## blood - blood TPM value
all_genes_total_both_brain_5 <- list()
all_genes_total_both_blood_5 <- list()

## genes that meet the TPM threshold in brain only
## gene in brain meets TPM threshold, gene in blood does not meet TPM threshold
all_genes_total_onlybrain_brain_5 <- list()
all_genes_total_onlybrain_blood_5 <- list()

## genes that meet the TPM threshold in blood only
all_genes_total_onlyblood_brain_5 <- list()
all_genes_total_onlyblood_blood_5 <- list()

## genes that meet the TPM threshold in neither blood or brain
all_genes_total_neither_brain_5 <- list()
all_genes_total_neither_blood_5 <- list()

num_sample_types <- 14
for (i in 1:num_sample_types){
  all_genes_total_in_blood_5[[i]] <- (rowSums(all_genes_all_matched_samples[[i]]>=5) == length(all_genes_all_matched_samples[[i]]))
  all_genes_total_in_brain_5[[i]] <- (rowSums(all_genes_all_matched_samples[[i+14]]>=5) == length(all_genes_all_matched_samples[[i+14]]))
  
  all_genes_total_blood_5[[i]] <- all_genes_all_matched_samples[[i]] %>% filter(all_genes_total_in_blood_5[[i]])
  all_genes_total_brain_5[[i]] <- all_genes_all_matched_samples[[i+14]] %>% filter(all_genes_total_in_brain_5[[i]])
  
  all_genes_total_blood_5_not[[i]] <- all_genes_all_matched_samples[[i]] %>% filter(!(all_genes_total_in_blood_5[[i]]))
  all_genes_total_brain_5_not[[i]] <- all_genes_all_matched_samples[[i+14]] %>% filter(!(all_genes_total_in_brain_5[[i]]))
  
  ## Both blood and brain
  all_genes_total_both_brain_5[[i]] <- all_genes_total_brain_5[[i]] %>% filter(all_genes_total_brain_5[[i]][[1]] %in% all_genes_total_blood_5[[i]][[1]])
  all_genes_total_both_blood_5[[i]] <- all_genes_total_blood_5[[i]] %>% filter(all_genes_total_blood_5[[i]][[1]] %in% all_genes_total_brain_5[[i]][[1]])
  
  ## Brain only 
  all_genes_total_onlybrain_brain_5[[i]] <- all_genes_total_brain_5[[i]] %>% filter(all_genes_total_brain_5[[i]][[1]] %in% all_genes_total_blood_5_not[[i]][[1]])
  all_genes_total_onlybrain_blood_5[[i]] <- all_genes_total_blood_5_not[[i]] %>% filter(all_genes_total_blood_5_not[[i]][[1]] %in% all_genes_total_brain_5[[i]][[1]])
  
  ## Blood only
  all_genes_total_onlyblood_brain_5[[i]] <- all_genes_total_brain_5_not[[i]] %>% filter(all_genes_total_brain_5_not[[i]][[1]] %in% all_genes_total_blood_5[[i]][[1]])
  all_genes_total_onlyblood_blood_5[[i]] <- all_genes_total_blood_5[[i]] %>% filter(all_genes_total_blood_5[[i]][[1]] %in% all_genes_total_brain_5_not[[i]][[1]])
  
  ## Neither blood or brain
  all_genes_total_neither_brain_5[[i]] <- all_genes_total_brain_5_not[[i]] %>% filter(all_genes_total_brain_5_not[[i]][[1]] %in% all_genes_total_blood_5_not[[i]][[1]])
  all_genes_total_neither_blood_5[[i]] <- all_genes_total_blood_5_not[[i]] %>% filter(all_genes_total_blood_5_not[[i]][[1]] %in% all_genes_total_brain_5_not[[i]][[1]])
  
}

#####
## File write out for TPM >= 5 only
## Directory and file name can be modified to include other TPMs
both_tpm_5 <- paste(all_genes_all_samples_tpm_5,"/both_tpm_5",sep="")
dir.create(both_tpm_5)

## only brain
onlybrain_tpm_5 <- paste(all_genes_all_samples_tpm_5, "/onlybrain_tpm_5", sep="")
dir.create(onlybrain_tpm_5)

## only blood
onlyblood_tpm_5 <- paste(all_genes_all_samples_tpm_5, "/onlyblood_tpm_5", sep="")
dir.create(onlyblood_tpm_5)

## neither
neither_tpm_5 <- paste(all_genes_all_samples_tpm_5,"/neither_tpm_5", sep="")
dir.create(neither_tpm_5)

for (i in 1:num_sample_types){
  ## both
  both_blood_tsv_5 <- paste(both_tpm_5, "/both_blood_",filenames[[i]],"_5.tsv", sep="")
  write.table(all_genes_total_both_blood_5[[i]], both_blood_tsv_5, row.names = FALSE, sep ="\t")
  
  both_brain_tsv_5 <- paste(both_tpm_5, "/both_brain_",filenames[[i+14]],"_5.tsv", sep="")
  write.table(all_genes_total_both_brain_5[[i]], both_brain_tsv_5, row.names = FALSE, sep ="\t")
  
  ## brain only
  onlybrain_blood_tsv_5 <- paste(onlybrain_tpm_5, "/onlybrain_blood_", filenames[[i]], "_5.tsv", sep="")
  write.table(all_genes_total_onlybrain_blood_5[[i]], onlybrain_blood_tsv_5, row.names =FALSE, sep ="\t")
  
  onlybrain_brain_tsv_5 <- paste(onlybrain_tpm_5, "/onlybrain_brain_", filenames[[i+14]], "_5.tsv", sep="")
  write.table(all_genes_total_onlybrain_brain_5[[i]], onlybrain_brain_tsv_5, row.names =FALSE, sep ="\t")
  
  ## blood only
  onlyblood_blood_tsv_5 <- paste(onlyblood_tpm_5, "/onlyblood_blood_", filenames[[i]], "_5.tsv", sep="")
  write.table(all_genes_total_onlyblood_blood_5[[i]], onlyblood_blood_tsv_5, row.names =FALSE, sep ="\t")
  
  onlyblood_brain_tsv_5 <- paste(onlyblood_tpm_5, "/onlyblood_brain_", filenames[[i+14]], "_5.tsv", sep="")
  write.table(all_genes_total_onlyblood_brain_5[[i]], onlyblood_brain_tsv_5, row.names =FALSE, sep ="\t")
  
  ## neither
  neither_blood_tsv_5 <- paste(neither_tpm_5,"/neither_blood_", filenames[[i]],"_5.tsv",sep="")
  write.table(all_genes_total_neither_blood_5[[i]], neither_blood_tsv_5, row.names=FALSE, sep="\t")
  
  neither_brain_tsv_5 <- paste(neither_tpm_5,"/neither_brain_",filenames[[i+14]],"_5.tsv",sep="")
  write.table(all_genes_total_neither_brain_5[[i]], neither_brain_tsv_5, row.names=FALSE, sep="\t")
  
}

#####
## Filter for TPM 10

# in_tissue_1 - shows logical output FALSE/TRUE, tells you which line shows the gene that match criteria
all_genes_total_in_brain_10 <- list()
all_genes_total_in_blood_10 <- list()

## genes with TPM >= 10 in individual brain tissue type across all samples
## continuation of in_tissue_10, this one shows the genes with their TPM values
all_genes_total_brain_10 <- list()

## genes that do not have TPM >= 10 across all samples in individual brain tissue type
all_genes_total_brain_10_not <- list()

## genes expressed in blood across all samples with TPM >= 10
all_genes_total_blood_10 <- list()

## genes that do not have TPM >= 10 across all samples in blood
all_genes_total_blood_10_not <- list()

## genes that meet the TPM threshold in BOTH brain and blood
## both - present in both brain and blood,
## brain - brain TPM value
## blood - blood TPM value
all_genes_total_both_brain_10 <- list()
all_genes_total_both_blood_10 <- list()

## genes that meet the TPM threshold in brain only
## gene in brain meets TPM threshold, gene in blood does not meet TPM threshold
all_genes_total_onlybrain_brain_10 <- list()
all_genes_total_onlybrain_blood_10 <- list()

## genes that meet the TPM threshold in blood only
all_genes_total_onlyblood_brain_10 <- list()
all_genes_total_onlyblood_blood_10 <- list()

## genes that meet the TPM threshold in neither blood or brain
all_genes_total_neither_brain_10 <- list()
all_genes_total_neither_blood_10 <- list()

num_sample_types <- 14
for (i in 1:num_sample_types){
  all_genes_total_in_blood_10[[i]] <- (rowSums(all_genes_all_matched_samples[[i]]>=10) == length(all_genes_all_matched_samples[[i]]))
  all_genes_total_in_brain_10[[i]] <- (rowSums(all_genes_all_matched_samples[[i+14]]>=10) == length(all_genes_all_matched_samples[[i+14]]))
  
  all_genes_total_blood_10[[i]] <- all_genes_all_matched_samples[[i]] %>% filter(all_genes_total_in_blood_10[[i]])
  all_genes_total_brain_10[[i]] <- all_genes_all_matched_samples[[i+14]] %>% filter(all_genes_total_in_brain_10[[i]])
  
  all_genes_total_blood_10_not[[i]] <- all_genes_all_matched_samples[[i]] %>% filter(!(all_genes_total_in_blood_10[[i]]))
  all_genes_total_brain_10_not[[i]] <- all_genes_all_matched_samples[[i+14]] %>% filter(!(all_genes_total_in_brain_10[[i]]))
  
  ## Both blood and brain
  all_genes_total_both_brain_10[[i]] <- all_genes_total_brain_10[[i]] %>% filter(all_genes_total_brain_10[[i]][[1]] %in% all_genes_total_blood_10[[i]][[1]])
  all_genes_total_both_blood_10[[i]] <- all_genes_total_blood_10[[i]] %>% filter(all_genes_total_blood_10[[i]][[1]] %in% all_genes_total_brain_10[[i]][[1]])
  
  ## Brain only 
  all_genes_total_onlybrain_brain_10[[i]] <- all_genes_total_brain_10[[i]] %>% filter(all_genes_total_brain_10[[i]][[1]] %in% all_genes_total_blood_10_not[[i]][[1]])
  all_genes_total_onlybrain_blood_10[[i]] <- all_genes_total_blood_10_not[[i]] %>% filter(all_genes_total_blood_10_not[[i]][[1]] %in% all_genes_total_brain_10[[i]][[1]])
  
  ## Blood only
  all_genes_total_onlyblood_brain_10[[i]] <- all_genes_total_brain_10_not[[i]] %>% filter(all_genes_total_brain_10_not[[i]][[1]] %in% all_genes_total_blood_10[[i]][[1]])
  all_genes_total_onlyblood_blood_10[[i]] <- all_genes_total_blood_10[[i]] %>% filter(all_genes_total_blood_10[[i]][[1]] %in% all_genes_total_brain_10_not[[i]][[1]])
  
  ## Neither blood or brain
  all_genes_total_neither_brain_10[[i]] <- all_genes_total_brain_10_not[[i]] %>% filter(all_genes_total_brain_10_not[[i]][[1]] %in% all_genes_total_blood_10_not[[i]][[1]])
  all_genes_total_neither_blood_10[[i]] <- all_genes_total_blood_10_not[[i]] %>% filter(all_genes_total_blood_10_not[[i]][[1]] %in% all_genes_total_brain_10_not[[i]][[1]])
  
}

#####
## File write out for TPM >= 10 only
## Directory and file name can be modified to include other TPMs
both_tpm_10 <- paste(all_genes_all_samples_tpm_10,"/both_tpm_10",sep="")
dir.create(both_tpm_10)

## only brain
onlybrain_tpm_10 <- paste(all_genes_all_samples_tpm_10, "/onlybrain_tpm_10", sep="")
dir.create(onlybrain_tpm_10)

## only blood
onlyblood_tpm_10 <- paste(all_genes_all_samples_tpm_10, "/onlyblood_tpm_10", sep="")
dir.create(onlyblood_tpm_10)

## neither
neither_tpm_10 <- paste(all_genes_all_samples_tpm_10,"/neither_tpm_10", sep="")
dir.create(neither_tpm_10)

for (i in 1:num_sample_types){
  ## both
  both_blood_tsv_10 <- paste(both_tpm_10, "/both_blood_",filenames[[i]],"_10.tsv", sep="")
  write.table(all_genes_total_both_blood_10[[i]], both_blood_tsv_10, row.names = FALSE, sep ="\t")
  
  both_brain_tsv_10 <- paste(both_tpm_10, "/both_brain_",filenames[[i+14]],"_10.tsv", sep="")
  write.table(all_genes_total_both_brain_10[[i]], both_brain_tsv_10, row.names = FALSE, sep ="\t")
  
  ## brain only
  onlybrain_blood_tsv_10 <- paste(onlybrain_tpm_10, "/onlybrain_blood_", filenames[[i]], "_10.tsv", sep="")
  write.table(all_genes_total_onlybrain_blood_10[[i]], onlybrain_blood_tsv_10, row.names =FALSE, sep ="\t")
  
  onlybrain_brain_tsv_10 <- paste(onlybrain_tpm_10, "/onlybrain_brain_", filenames[[i+14]], "_10.tsv", sep="")
  write.table(all_genes_total_onlybrain_brain_10[[i]], onlybrain_brain_tsv_10, row.names =FALSE, sep ="\t")
  
  ## blood only
  onlyblood_blood_tsv_10 <- paste(onlyblood_tpm_10, "/onlyblood_blood_", filenames[[i]], "_10.tsv", sep="")
  write.table(all_genes_total_onlyblood_blood_10[[i]], onlyblood_blood_tsv_10, row.names =FALSE, sep ="\t")
  
  onlyblood_brain_tsv_10 <- paste(onlyblood_tpm_10, "/onlyblood_brain_", filenames[[i+14]], "_10.tsv", sep="")
  write.table(all_genes_total_onlyblood_brain_10[[i]], onlyblood_brain_tsv_10, row.names =FALSE, sep ="\t")
  
  ## neither
  neither_blood_tsv_10 <- paste(neither_tpm_10,"/neither_blood_", filenames[[i]],"_10.tsv",sep="")
  write.table(all_genes_total_neither_blood_10[[i]], neither_blood_tsv_10, row.names=FALSE, sep="\t")
  
  neither_brain_tsv_10 <- paste(neither_tpm_10,"/neither_brain_",filenames[[i+14]],"_10.tsv",sep="")
  write.table(all_genes_total_neither_brain_10[[i]], neither_brain_tsv_10, row.names=FALSE, sep="\t")
  
}


