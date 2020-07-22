## -----------------
## 
## Script name: gtex_tpm_sad.R
## sad: acronym for sex, age, death
## refers to the document for filtering: sex_age_death
##
## Purpose of script:
## 1) Filter GTEx subset files by:
##  - TPM threshold
##  - sex
##  - age
## 
## Documents needed to run this script:
## 1) Individual tissue type TPM tsv files
##  - Blood tsv file should be named as WB or whole_blood
##  - This is to make sure the file is imported last for script purpose
##
## 2) Sample IDs with sex and age info, downloaded from GTEx portal
##
## Note:
## 1) Whole blood is included as the last document in all lists that include individual brain tissues
## 2) This script uses for loop instead of the apply family. 
## -----------------

#####
## Set working directory
## Working directory should contain the documents listed above

current_dir <- "/Users/Sloth/Downloads/SexChrLab/BloodBrain"
setwd(current_dir)

## Set output directory
out_folder <- paste(getwd(),"/OutputFile",sep="")
# Don't run dir.create if folder has been created
dir.create(out_folder)

## Load applicable library
library(tidyverse)

#####
## Identifies all files with the pattern _tpm.tsv
my_files <- list.files(current_dir,pattern="*_tpm.tsv")

## Reads my_files into R environment
all_tsv <- lapply(my_files, read_tsv)

## Name all_tsv with their respective blood/brain tissue type
names(all_tsv) <- gsub(".tsv", "",
                       list.files(my_files, full.names=FALSE),
                       fixed=TRUE)

## Create a list of blood/brain tissue type
filenames <- names(all_tsv)

#####
## Rename GTEx sample IDs to only patient IDs
## This is because other downloaded GTEx files (age and sex info) contain only patient IDs, no sample IDs
## Current sample IDs contain extra information such as tissue type
## Patient ID format: GTEX-AAAA(A)
## Sample ID format: GTEX-AAAA(A)-BBBB-...
## The first 9 or 10 letters of sample ID is the patient ID

## Substring to the first 10 letters, then remove the last hyphen if present
cnames <- colnames(all_tsv)

all_tsv_try <- lapply(colnames(all_tsv), substr, 1, 10)

for (i in 1:length(all_tsv)){
  colnames(all_tsv[[i]]) <- substr(colnames(all_tsv[[i]]),1,10)
  colnames(all_tsv[[i]]) <- gsub("\\-$","",colnames(all_tsv[[i]]))
}

#####
## Read file that contains sex and age info
## 1 is male, 2 is female
sex_age_death <- read.table("gtex_sex_age_death.txt", sep = "\t", header = TRUE)

sample_male_age <- filter(sex_age_death, sex_age_death[[2]] == "1")
sample_female_age <- filter(sex_age_death, sex_age_death[[2]] == "2")

## Current all_tsv file contain samples that do not have corresponding sex info
## Since sex info is important for our analysis, we remove any samples that do not contain this info
## Intersect with all_tsv file to get samples with sex info
total_samples <- list()

male_samples <- list()
female_samples <- list()

num_sample_total <-  data.frame(matrix(ncol=length(all_tsv),nrow=1))
num_sample_male <-  data.frame(matrix(ncol=length(all_tsv),nrow=1))
num_sample_female <-  data.frame(matrix(ncol=length(all_tsv),nrow=1))

Name <- all_tsv[[1]][[1]]
Description <- all_tsv[[1]][[2]]

for (i in 1:length(all_tsv)){
  total_samples[[i]] <- all_tsv[[i]][,intersect(colnames(all_tsv[[i]]), sex_age_death[[1]])]
  total_samples[[i]] <- add_column(total_samples[[i]], Name, Description, .before = 1)
  num_sample_total[[i]] <- length(total_samples[[i]]) - 2
  
  male_samples[[i]] <- all_tsv[[i]][,intersect(colnames(all_tsv[[i]]), sample_male_age[[1]])]
  male_samples[[i]] <- add_column(male_samples[[i]], Name, Description, .before = 1)
  num_sample_male[[i]] <- length(male_samples[[i]]) - 2
  
  female_samples[[i]] <- all_tsv[[i]][,intersect(colnames(all_tsv[[i]]), sample_female_age[[1]])]
  female_samples[[i]] <- add_column(female_samples[[i]], Name, Description, .before = 1)
  num_sample_female[[i]] <- length(female_samples[[i]]) - 2
}

matched_sad_folder <- paste(out_folder,"/matched_sad_tissue",sep="")
dir.create(matched_sad_folder)

matched_tissue_num_sample <-as.data.frame(t(rbind(num_sample_total, num_sample_male, num_sample_female)), row.names = filenames)
colnames(matched_tissue_num_sample)=c("Total samples","Male samples","Female samples")
write.table(matched_tissue_num_sample, paste0(matched_sad_folder, "/matched_tissue_num_sample.tsv"), sep="\t", col.names=NA)

#####
## We want to make sure each blood tissue sample also has a corresponding brain tissue sample
## The following compares samples in both blood and individual tissue and keeps only samples present in both
## Matched indicates all samples in brain tissue are also present in blood

matched_total_brain <- list()
matched_total_blood <- list()

matched_male_brain <- list()
matched_male_blood <- list()

matched_female_brain <- list()
matched_female_blood <- list()

num_matched_sample_total <-  data.frame(matrix(ncol=length(all_tsv),nrow=1))
num_matched_sample_male <-  data.frame(matrix(ncol=length(all_tsv),nrow=1))
num_matched_sample_female <-  data.frame(matrix(ncol=length(all_tsv),nrow=1))

## Note that there are 14 samples in total, and the 14th sample is whole blood

all_genes_folder <- paste(matched_sad_folder, "/all_genes", sep = "")
dir.create(all_genes_folder)

all_genes_all_samples_folder <- paste(all_genes_folder, "/all_samples", sep = "")
dir.create(all_genes_all_samples_folder)

all_genes_male_samples_folder <- paste(all_genes_folder, "/male_samples", sep = "")
dir.create(all_genes_male_samples_folder)

all_genes_female_samples_folder <- paste(all_genes_folder, "/female_samples", sep = "")
dir.create(all_genes_female_samples_folder)

for (i in 1:length(all_tsv)){
  ## total
  matched_total_brain[[i]] <- total_samples[[i]][,intersect(colnames(total_samples[[i]]),colnames(total_samples[[14]]))]
  matched_total_blood[[i]] <- total_samples[[14]][,intersect(colnames(total_samples[[i]]),colnames(total_samples[[14]]))]
  num_matched_sample_total[[i]] <- length(matched_total_brain[[i]]) - 2
  
  matched_total_brain_tsv <- paste0(all_genes_all_samples_folder,"/matched_all_genes_all_brain_", filenames[[i]],".tsv")
  matched_total_blood_tsv <- paste0(all_genes_all_samples_folder,"/matched_all_genes_all_blood_", filenames[[i]],".tsv")
  
  write.table(matched_total_brain[[i]], matched_total_brain_tsv, row.names = FALSE, sep = "\t")
  write.table(matched_total_blood[[i]], matched_total_blood_tsv, row.names = FALSE, sep = "\t")
  
  ## male samples
  matched_male_brain[[i]] <- male_samples[[i]][,intersect(colnames(male_samples[[i]]),colnames(male_samples[[14]]))]
  matched_male_blood[[i]] <- male_samples[[14]][,intersect(colnames(male_samples[[i]]),colnames(male_samples[[14]]))]
  num_matched_sample_male[[i]] <- length(matched_male_brain[[i]]) - 2
  
  matched_male_brain_tsv <- paste0(all_genes_male_samples_folder,"/matched_all_genes_male_brain_", filenames[[i]],".tsv")
  matched_male_blood_tsv <- paste0(all_genes_male_samples_folder,"/matched_all_genes_male_blood_", filenames[[i]],".tsv")
  
  write.table(matched_male_brain[[i]], matched_male_brain_tsv, row.names = FALSE, sep = "\t")
  write.table(matched_male_blood[[i]], matched_male_blood_tsv, row.names = FALSE, sep = "\t")
  
  ## female samples
  matched_female_brain[[i]] <- female_samples[[i]][,intersect(colnames(female_samples[[i]]),colnames(female_samples[[14]]))]
  matched_female_blood[[i]] <- female_samples[[14]][,intersect(colnames(female_samples[[i]]),colnames(female_samples[[14]]))]
  num_matched_sample_female[[i]] <- length(matched_female_brain[[i]]) - 2
  
  matched_female_brain_tsv <- paste0(all_genes_female_samples_folder,"/matched_all_genes_female_brain_", filenames[[i]],".tsv")
  matched_female_blood_tsv <- paste0(all_genes_female_samples_folder,"/matched_all_genes_female_blood_", filenames[[i]],".tsv")
  
  write.table(matched_female_brain[[i]], matched_female_brain_tsv, row.names = FALSE, sep = "\t")
  write.table(matched_female_blood[[i]], matched_female_blood_tsv, row.names = FALSE, sep = "\t")
  
}

num_before_after_matched <-as.data.frame(t(rbind(num_sample_total, num_matched_sample_total, num_sample_male, num_matched_sample_male, num_sample_female, num_matched_sample_female)), row.names = filenames)
colnames(num_before_after_matched)=c("Before Total", "After Total", "Before Male","After Male", "Before Female", "After Female")
write.table(num_before_after_matched, paste0(matched_sad_folder, "/num_sample_before_after_matched.tsv"), sep="\t", col.names=NA)

#####
## Filter by protein coding gene after filtering for samples with sex info

## Filter for protein coding genes only, no TPM threshold yet
## Read file with protein coding genes
## Filter by ID for now, it seems like there's more if filtered by names

protein_coding_id <- read.table("pcg_mod_id.txt", sep = "\t", header = FALSE)

pcg_id_matched_total_brain <- list()
pcg_id_matched_total_blood <- list()

pcg_id_matched_male_brain <- list()
pcg_id_matched_male_blood <- list()

pcg_id_matched_female_brain <- list()
pcg_id_matched_female_blood <- list()

pcg_folder <- paste(matched_sad_folder,"/protein_coding",sep="")
dir.create(pcg_folder)

pcg_all_samples <- paste(pcg_folder, "/all_samples", sep = "")
dir.create(pcg_all_samples)

pcg_male_samples <- paste(pcg_folder, "/male_samples", sep = "")
dir.create(pcg_male_samples)

pcg_female_samples <- paste(pcg_folder, "/female_samples", sep = "")
dir.create(pcg_female_samples)

for (i in 1:length(all_tsv)){
  ## total
  pcg_id_matched_total_brain[[i]] <- matched_total_brain[[i]] %>% filter(matched_total_brain[[i]][[1]] %in% protein_coding_id[[1]])
  pcg_id_matched_total_blood[[i]] <- matched_total_blood[[i]] %>% filter(matched_total_blood[[i]][[1]] %in% protein_coding_id[[1]])
  
  pcg_id_matched_total_brain_tsv <- paste0(pcg_all_samples,"/matched_pcg_all_brain_", filenames[[i]],".tsv")
  pcg_id_matched_total_blood_tsv <- paste0(pcg_all_samples,"/matched_pcg_all_blood_", filenames[[i]],".tsv")
  
  write.table(pcg_id_matched_total_brain[[i]], pcg_id_matched_total_brain_tsv, row.names = FALSE, sep = "\t")
  write.table(pcg_id_matched_total_blood[[i]], pcg_id_matched_total_blood_tsv, row.names = FALSE, sep = "\t")
  
  ## male samples
  pcg_id_matched_male_brain[[i]] <- matched_male_brain[[i]] %>% filter(matched_male_brain[[i]][[1]] %in% protein_coding_id[[1]])
  pcg_id_matched_male_blood[[i]] <- matched_male_blood[[i]] %>% filter(matched_male_blood[[i]][[1]] %in% protein_coding_id[[1]])
  
  pcg_id_matched_male_brain_tsv <- paste0(pcg_male_samples,"/matched_pcg_male_brain_", filenames[[i]],".tsv")
  pcg_id_matched_male_blood_tsv <- paste0(pcg_male_samples,"/matched_pcg_male_blood_", filenames[[i]],".tsv")
  
  write.table(pcg_id_matched_male_brain[[i]], pcg_id_matched_male_brain_tsv, row.names = FALSE, sep = "\t")
  write.table(pcg_id_matched_male_blood[[i]], pcg_id_matched_male_blood_tsv, row.names = FALSE, sep = "\t")
  
  ## female samples
  pcg_id_matched_female_brain[[i]] <- matched_female_brain[[i]] %>% filter(matched_female_brain[[i]][[1]] %in% protein_coding_id[[1]])
  pcg_id_matched_female_blood[[i]] <- matched_female_blood[[i]] %>% filter(matched_female_blood[[i]][[1]] %in% protein_coding_id[[1]])
  
  pcg_id_matched_female_brain_tsv <- paste0(pcg_female_samples,"/matched_pcg_female_brain_", filenames[[i]],".tsv")
  pcg_id_matched_female_blood_tsv <- paste0(pcg_female_samples,"/matched_pcg_female_blood_", filenames[[i]],".tsv")
  
  write.table(pcg_id_matched_female_brain[[i]], pcg_id_matched_female_brain_tsv, row.names = FALSE, sep = "\t")
  write.table(pcg_id_matched_female_blood[[i]], pcg_id_matched_female_blood_tsv, row.names = FALSE, sep = "\t")
}