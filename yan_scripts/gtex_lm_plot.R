## -----------------
## 
## Script name: gtex_lm_plot.R
##
## Purpose of script:
## 1) Perform linear regression analysis after filtering for TPM
##  - All genes
##  - Candidate genes
## 
## 2) Perform multiple linear regression
##  - Sex info located at 3rd last column of sad_merged_genes_TPM
##  - Age info located at 2nd last column of sad_merged_genes_TPM
##
## Documents needed to run this script:
## 1) Gene lists from gtex_tpm_sad.R
## - sample filename forms: 
##    a) both_blood_brainregion_tpm_#.tsv
##    b) both_brain_brainregion_tpm_#.tsv
## - 28 files for each TPM, 84 files in total
##
## Note: 
## 1) Currently, working directory needs to be changed manually to include multiple TPMs
## 2) This script uses for loops
## -----------------

#####
## The following outputs genes with TPM >= 10 only
## Note that TPM >= 10 is used here instead of TPM >= 1
## This is because there are less genes with TPM >= 10, 
## To allow outputs with other TPMs, change directory and TPM values where necessary

## Read files from existing matched tissue

## Set working directory
## Working directory should include the documents listed above
working_dir <- "/Users/Sloth/Downloads/SexChrLab/BloodBrain/OutputFile/matched_sad_tissue/all_genes/all_samples/gene_list_1/both_tpm_1"
setwd(working_dir)

## Set output directory
out_folder <- paste(getwd(),"/plots",sep="")
# Don't run dir.create if folder has been created
dir.create(out_folder)

## Load applicable library
library(tidyverse)

## sex info
sex_age_death <- read.table("gtex_sex_age_death.txt", sep = "\t", header = TRUE)

sample_male_age <- filter(sex_age_death, sex_age_death[[2]] == "1")
sample_female_age <- filter(sex_age_death, sex_age_death[[2]] == "2")

#####
## TPM >= 1
tsv_genes_in_both_1 <- list.files(working_dir, pattern = ".tsv")
## Reads my_files into R environment
genes_in_both_1 <- lapply(tsv_genes_in_both_1, read_tsv)

## Name all_tsv with their respective blood/brain tissue type
names(genes_in_both_1) <- gsub(".tsv", "",
                               list.files(working_dir, pattern = ".tsv", full.names=FALSE),
                               fixed=TRUE)

## Create a list of blood/brain tissue type
filenames_genes_both_1 <- names(genes_in_both_1)

filenames <- sapply(strsplit(filenames_genes_both_1, "_"), function(x) x[8])

#####
## plotting graphs
blood_genes_1 <- list()
mod_blood_genes_1 <- list()

brain_genes_1 <- list()
mod_brain_genes_1 <- list()

merged_genes_1 <- list()

sad_column_1 <- list()
sad_merged_genes_1 <- list()

num_tissue_type <- 14

for (i in 1:num_tissue_type){
  blood_genes_1[[i]] <- genes_in_both_1[[i]][-c(1,2)]
  mod_blood_genes_1[[i]] <- as.data.frame(t(blood_genes_1[[i]]))
  names(mod_blood_genes_1[[i]]) <- paste0(genes_in_both_1[[i]][[2]],"_blood")
  
  brain_genes_1[[i]] <- genes_in_both_1[[i+ num_tissue_type]][-c(1,2)]
  mod_brain_genes_1[[i]] <- as.data.frame(t(brain_genes_1[[i]]))
  names(mod_brain_genes_1[[i]]) <- paste0(genes_in_both_1[[i+num_tissue_type]][[2]],"_brain")
  
  merged_genes_1[[i]] <- bind_cols(mod_blood_genes_1[[i]], mod_brain_genes_1[[i]])
  
  sad_column_1[[i]] <- sex_age_death %>% filter(sex_age_death[[1]] %in% rownames(merged_genes_1[[i]]))
  sad_merged_genes_1[[i]] <- bind_cols(merged_genes_1[[i]], sad_column_1[[i]])
  
}

gene_names <- sapply(strsplit(colnames(sad_merged_genes_1[[1]]), "_"), function(x) x[1])

#####
## Linear regression
gene_plot_1 <- list()
lin_reg_1 <- list()
graph_lin_reg_1 <- list()

for (i in 1:num_tissue_type){
  lin_reg_1[[i]] <- list()
  gene_plot_1[[i]] <- list()
  graph_lin_reg_1[[i]] <- list()
  
  pdfnames <- paste0(filenames[[i]],"_tpm_1.pdf")
  
  pdf(file=pdfnames)
  for (j in 1:((length(sad_merged_genes_1[[i]])-4)/2)){
    
    layout(matrix(c(1,1,2,3,4,1), nrow = 3, byrow = TRUE))
    par(mar=c(4,4.5,1,2), mgp = c(2,1,0), cex = 1)
    graph_lin_reg_1[[i]][[j]] <- plot(sad_merged_genes_1[[i]][[j]], sad_merged_genes_1[[i]][[j +length(merged_genes_1[[i]])/2]],  xlab = filenames[14], ylab = filenames[[i]], main = gene_names[[j]])
    lin_reg_1[[i]][[j]] <- lm(sad_merged_genes_1[[i]][[j +length(merged_genes_1[[i]])/2]] ~ sad_merged_genes_1[[i]][[j]])
    ## fix abline, abline does not look correct in graph
    abline(lin_reg_1[[i]][[j]])
    plot(lin_reg_1[[i]][[j]])
    
  }
  dev.off()
  
}

#####
## TPM >= 5
tsv_genes_in_both_5 <- list.files(working_dir, pattern = ".tsv")
## Reads my_files into R environment
genes_in_both_5 <- lapply(tsv_genes_in_both_5, read_tsv)

## Name all_tsv with their respective blood/brain tissue type
names(genes_in_both_5) <- gsub(".tsv", "",
                                list.files(working_dir, pattern = ".tsv", full.names=FALSE),
                                fixed=TRUE)

## Create a list of blood/brain tissue type
filenames_genes_both_5 <- names(genes_in_both_5)

filenames <- sapply(strsplit(filenames_genes_both_5, "_"), function(x) x[8])

## sex info
sex_age_death <- read.table("gtex_sex_age_death.txt", sep = "\t", header = TRUE)

sample_male_age <- filter(sex_age_death, sex_age_death[[2]] == "1")
sample_female_age <- filter(sex_age_death, sex_age_death[[2]] == "2")

#####
## plotting graphs
blood_genes_5 <- list()
mod_blood_genes_5 <- list()

brain_genes_5 <- list()
mod_brain_genes_5 <- list()

merged_genes_5 <- list()

sad_column_5 <- list()
sad_merged_genes_5 <- list()

num_tissue_type <- 14

for (i in 1:num_tissue_type){
  blood_genes_5[[i]] <- genes_in_both_5[[i]][-c(1,2)]
  mod_blood_genes_5[[i]] <- as.data.frame(t(blood_genes_5[[i]]))
  names(mod_blood_genes_5[[i]]) <- paste0(genes_in_both_5[[i]][[2]],"_blood")
  
  brain_genes_5[[i]] <- genes_in_both_5[[i+ num_tissue_type]][-c(1,2)]
  mod_brain_genes_5[[i]] <- as.data.frame(t(brain_genes_5[[i]]))
  names(mod_brain_genes_5[[i]]) <- paste0(genes_in_both_5[[i+num_tissue_type]][[2]],"_brain")
  
  merged_genes_5[[i]] <- bind_cols(mod_blood_genes_5[[i]], mod_brain_genes_5[[i]])
  
  sad_column_5[[i]] <- sex_age_death %>% filter(sex_age_death[[1]] %in% rownames(merged_genes_5[[i]]))
  sad_merged_genes_5[[i]] <- bind_cols(merged_genes_5[[i]], sad_column_5[[i]])
  
}

gene_names <- sapply(strsplit(colnames(sad_merged_genes_5[[1]]), "_"), function(x) x[1])

#####
## Linear regression
gene_plot_5 <- list()
lin_reg_5 <- list()
graph_lin_reg_5 <- list()

for (i in 1:num_tissue_type){
  lin_reg_5[[i]] <- list()
  gene_plot_5[[i]] <- list()
  graph_lin_reg_5[[i]] <- list()
  
  pdfnames <- paste0(filenames[[i]],"_tpm_5.pdf")
  
  pdf(file=pdfnames)
  for (j in 1:((length(sad_merged_genes_5[[i]])-4)/2)){
    
    layout(matrix(c(1,1,2,3,4,5), nrow = 3, byrow = TRUE))
    par(mar=c(4,4.5,1,2), mgp = c(2,1,0), cex = 1)
    graph_lin_reg_5[[i]][[j]] <- plot(sad_merged_genes_5[[i]][[j]], sad_merged_genes_5[[i]][[j +length(merged_genes_5[[i]])/2]],  xlab = filenames[14], ylab = filenames[[i]], main = gene_names[[j]])
    lin_reg_5[[i]][[j]] <- lm(sad_merged_genes_5[[i]][[j +length(merged_genes_5[[i]])/2]] ~ sad_merged_genes_5[[i]][[j]])
    ## fix abline, abline does not look correct in graph
    abline(lin_reg_5[[i]][[j]])
    plot(lin_reg_5[[i]][[j]])
    
  }
  dev.off()
  
}

#####
## TPM >= 10
## Identifies all files with the pattern _tpm.tsv
tsv_genes_in_both_10 <- list.files(working_dir, pattern = ".tsv")
## Reads my_files into R environment
genes_in_both_10 <- lapply(tsv_genes_in_both_10, read_tsv)

## Name all_tsv with their respective blood/brain tissue type
names(genes_in_both_10) <- gsub(".tsv", "",
                       list.files(working_dir, pattern = ".tsv", full.names=FALSE),
                       fixed=TRUE)

## Create a list of blood/brain tissue type
filenames_genes_both_10 <- names(genes_in_both_10)

filenames <- sapply(strsplit(filenames_genes_both_10, "_"), function(x) x[8])

#####
## plotting graphs - ggplots2
blood_genes_10 <- list()
mod_blood_genes_10 <- list()

brain_genes_10 <- list()
mod_brain_genes_10 <- list()

merged_genes_10 <- list()

sad_column_10 <- list()
sad_merged_genes_10 <- list()

num_tissue_type <- 14

for (i in 1:num_tissue_type){
  blood_genes_10[[i]] <- genes_in_both_10[[i]][-c(1,2)]
  mod_blood_genes_10[[i]] <- as.data.frame(t(blood_genes_10[[i]]))
  names(mod_blood_genes_10[[i]]) <- paste0(genes_in_both_10[[i]][[2]],"_blood")
  
  brain_genes_10[[i]] <- genes_in_both_10[[i+ num_tissue_type]][-c(1,2)]
  mod_brain_genes_10[[i]] <- as.data.frame(t(brain_genes_10[[i]]))
  names(mod_brain_genes_10[[i]]) <- paste0(genes_in_both_10[[i+num_tissue_type]][[2]],"_brain")
  
  merged_genes_10[[i]] <- bind_cols(mod_blood_genes_10[[i]], mod_brain_genes_10[[i]])
  
  sad_column_10[[i]] <- sex_age_death %>% filter(sex_age_death[[1]] %in% rownames(merged_genes_10[[i]]))
  sad_merged_genes_10[[i]] <- bind_cols(merged_genes_10[[i]], sad_column_10[[i]])

}

gene_names <- sapply(strsplit(colnames(sad_merged_genes_10[[1]]), "_"), function(x) x[1])

#####
## Linear regression
gene_plot_10 <- list()
lin_reg_10 <- list()
graph_lin_reg_10 <- list()

for (i in 1:num_tissue_type){
  lin_reg_10[[i]] <- list()
  gene_plot_10[[i]] <- list()
  graph_lin_reg_10[[i]] <- list()
  
  pdfnames <- paste0(filenames[[i]],"_tpm_10.pdf")

  pdf(file=pdfnames)
  for (j in 1:((length(sad_merged_genes_10[[i]])-4)/2)){

    layout(matrix(c(1,1,2,3,4,5), nrow = 3, byrow = TRUE))
    par(mar=c(4,4.5,1,2), mgp = c(2,1,0), cex = 1)
    graph_lin_reg_10[[i]][[j]] <- plot(sad_merged_genes_10[[i]][[j]], sad_merged_genes_10[[i]][[j +length(merged_genes_10[[i]])/2]],  xlab = filenames[14], ylab = filenames[[i]], main = gene_names[[j]])
    lin_reg_10[[i]][[j]] <- lm(sad_merged_genes_10[[i]][[j +length(merged_genes_10[[i]])/2]] ~ sad_merged_genes_10[[i]][[j]])
    ## fix abline, abline does not look correct in graph
    abline(lin_reg_10[[i]][[j]])
    plot(lin_reg_10[[i]][[j]])

  }
  dev.off()
  
}