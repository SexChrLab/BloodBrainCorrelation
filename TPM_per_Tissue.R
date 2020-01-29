#This script identifies the number of protein-coding genes within each brain region of interest and whole blood within the GTEx data set. It then goes to identify number of genes at certain expression levels (TPM >=1, >=5, >=10). Finally, it identifies the number of missed genes in the blood based on the brain patients that have expression in the blood (ie Blood1 = Cortex1). This also examines the missed genes at TPM >1, 5, and 10.

#set working directory
setwd("/Users/austinevanovich/Documents/GTEx_Blood_Brain")

#Read in data
cortex <- read.table("Cortex_Counts.tsv", header = T)

substantia <- read.table("SubstantiaNigra_Counts.tsv", header = T)

spinalcord <- read.table("SpinalCord_Counts.tsv", header = T)

putamen <- read.table("PutamenBasalGanglia_Counts.tsv", header = T)

nucleus <- read.table("NucleusAccumbensBasalGanglia_Counts.tsv", header = T)

hypothalamus <- read.table("Hypothalamus_Counts.tsv", header = T)

hippocampus <- read.table("Hippocampus_Counts.tsv", header = T)

frontalcortex <- read.table("FrontalCortex_Counts.tsv", header = T)

cerebellum <- read.table("Cerebellum_Counts.tsv", header = T)

cerebellarhemisphere <- read.table("CerebellarHemisphere_Counts.tsv", header = T)

caudate <- read.table("CaudateBasalGanglia_Counts.tsv", header = T)

anterior <- read.table("AnteriorCingulateCortex_Counts.tsv", header = T)

amygdala <- read.table("Amygdala_Counts.tsv", header = T)

wholeblood <- read.table("WholeBlood_Counts.tsv", header = T)

#read in data for the blood samples that correspond to their respective brain sample
blood_cortex <- read.table("Blood_Cortex_counts.tsv", header = T)

blood_substantia <- read.table("Blood_Substantia_counts.tsv", header = T)

blood_spinalcord <- read.table("Blood_SpinalCord_counts.tsv", header = T)

blood_putamen <- read.table("Blood_Putamen_counts.tsv", header = T)

blood_nucleus <- read.table("Blood_Nucleus_counts.tsv", header = T)

blood_hypothalamus <- read.table("Blood_Hypothalamus_counts.tsv", header = T)

blood_hippocampus <- read.table("Blood_Hippocampus_counts.tsv", header = T)

blood_frontalcortex <- read.table("Blood_FrontalCortex_counts.tsv", header = T)

blood_cerebellum <- read.table("Blood_Cerebellum_counts.tsv", header = T)

blood_cerebellarhemisphere <- read.table("Blood_CerebellarHemisphere_counts.tsv", header = T)

blood_caudate <- read.table("Blood_Caudate_counts.tsv", header = T)

blood_anterior <- read.table("Blood_Anterior_counts.tsv", header = T)

blood_amygdala <- read.table("Blood_Amygdala_counts.tsv", header = T)

######Limit each sample to only protein coding genes

#BiocManager::install("dplyr")
library(dplyr)

protein_ids <- read.csv("mart_export.txt", header = T)

colnames(protein_ids)
########################
#Subsetting each count object to just the protein-coding genes by their ensembl gene ID: http://www.ensembl.org/biomart/martview/58e679a255c6e23f443459a2779f080e
########################
cortex_protein <- subset(cortex, cortex$Name %in% protein_ids$Gene.stable.ID.version)

substantia_protein <- subset(substantia, substantia$Name %in% protein_ids$Gene.stable.ID.version)

spinalcord_protein <- subset(spinalcord, spinalcord$Name %in% protein_ids$Gene.stable.ID.version)

putamen_protein <- subset(putamen, putamen$Name %in% protein_ids$Gene.stable.ID.version)

nucleus_protein <- subset(nucleus, nucleus$Name %in% protein_ids$Gene.stable.ID.version)

hypothalamus_protein <- subset(hypothalamus, hypothalamus$Name %in% protein_ids$Gene.stable.ID.version)

hippocampus_protein <- subset(hippocampus, hippocampus$Name %in% protein_ids$Gene.stable.ID.version)

frontalcortex_protein <- subset(frontalcortex, frontalcortex$Name %in% protein_ids$Gene.stable.ID.version)

cerebellum_protein <- subset(cerebellum, cerebellum$Name %in% protein_ids$Gene.stable.ID.version)

cerebellarhemisphere_protein <- subset(cerebellarhemisphere, cerebellarhemisphere$Name %in% protein_ids$Gene.stable.ID.version)

caudate_protein <- subset(caudate, caudate$Name %in% protein_ids$Gene.stable.ID.version)

anterior_protein <- subset(anterior, anterior$Name %in% protein_ids$Gene.stable.ID.version)

amygdala_protein <- subset(amygdala, amygdala$Name %in% protein_ids$Gene.stable.ID.version)

wholeblood_protein <- subset(wholeblood, wholeblood$Name %in% protein_ids$Gene.stable.ID.version)

#do the above again, but for the blood samples that correspond to the brain samples
blood_cortex_protein <- subset(blood_cortex, blood_cortex$Name %in% protein_ids$Gene.stable.ID.version)

blood_substantia_protein <- subset(blood_substantia, blood_substantia$Name %in% protein_ids$Gene.stable.ID.version)

blood_spinalcord_protein <- subset(blood_spinalcord, blood_spinalcord$Name %in% protein_ids$Gene.stable.ID.version)

blood_putamen_protein <- subset(blood_putamen, blood_putamen$Name %in% protein_ids$Gene.stable.ID.version)

blood_nucleus_protein <- subset(blood_nucleus, blood_nucleus$Name %in% protein_ids$Gene.stable.ID.version)

blood_hypothalamus_protein <- subset(blood_hypothalamus, blood_hypothalamus$Name %in% protein_ids$Gene.stable.ID.version)

blood_hippocampus_protein <- subset(blood_hippocampus, blood_hippocampus$Name %in% protein_ids$Gene.stable.ID.version)

blood_frontalcortex_protein <- subset(blood_frontalcortex, blood_frontalcortex$Name %in% protein_ids$Gene.stable.ID.version)

blood_cerebellum_protein <- subset(blood_cerebellum, blood_cerebellum$Name %in% protein_ids$Gene.stable.ID.version)

blood_cerebellarhemisphere_protein <- subset(blood_cerebellarhemisphere, blood_cerebellarhemisphere$Name %in% protein_ids$Gene.stable.ID.version)

blood_caudate_protein <- subset(blood_caudate, blood_caudate$Name %in% protein_ids$Gene.stable.ID.version)

blood_anterior_protein <- subset(blood_anterior, blood_anterior$Name %in% protein_ids$Gene.stable.ID.version)

blood_amygdala_protein <- subset(blood_amygdala, blood_amygdala$Name %in% protein_ids$Gene.stable.ID.version)

#now with those, get the intersection of genes greater than 1TPM, 5TPM, 10 TPM within each tissue - this is to identify any missed genes
#first, subset to just the GTEx patient IDs
#then only sum the rows based on the condition you want
#merge the two objects together and select for the count columns you want to see
#then, based on the selected object, only keep rows with specific values

####TPM >1
brain2 <- cortex_protein[,-1]
brain3 <- brain2[,-1]  #Use this value to determine number of patients in the brain tissue

blood2 <- blood_cortex_protein[,-1]
blood3 <- blood2[,-1] #and this to determine number of patients in the brain-blood samples

blood_cortex_protein$bl_count_1 <- rowSums(blood3 >= 1)
blood_cortex_1 <- which(blood_cortex_protein$bl_count_1 >= 210)

cortex_protein$su_count_1 <- rowSums(brain3 >= 1)
cortex_1 <- which(cortex_protein$su_count_1 >= 215)
#View(blood_cortex_1)

blood_cortex_merged <- merge(blood_cortex_protein, cortex_protein, by = "Name")
blsu_count1 <- select(blood_cortex_merged, bl_count_1, su_count_1)
#View(blco_count1)

tpm1 <- blsu_count1[blsu_count1$bl_count_1 == 210 & blsu_count1$su_count_1 == 215,][,1:2]
tpm1 <- blsu_count1[blsu_count1$bl_count_1 == 210 & blsu_count1$su_count_1 < 215,][,1:2]
tpm1 <- blsu_count1[blsu_count1$bl_count_1 < 210 & blsu_count1$su_count_1 == 215,][,1:2]
#View(bcor2)

####For TPM >5
blood_cortex_protein$bl_count_5 <- rowSums(blood3 >= 5)
blood_cortex_5 <- which(blood_cortex_protein$bl_count_5 >= 210)

cortex_protein$su_count_5 <- rowSums(brain3 >= 5)
cortex_5 <- which(cortex_protein$su_count_5 >= 215)

blood_cortex_merged <- merge(blood_cortex_protein, cortex_protein, by = "Name")
blsu_count5 <- select(blood_cortex_merged, bl_count_5, su_count_5)
#View(blco_count1)

tpm5 <- blsu_count5[blsu_count5$bl_count_5 == 210 & blsu_count5$su_count_5 == 215,][,1:2]
tpm5 <- blsu_count5[blsu_count5$bl_count_5 == 210 & blsu_count5$su_count_5 < 215,][,1:2]
tpm5 <- blsu_count5[blsu_count5$bl_count_5 < 210 & blsu_count5$su_count_5 == 215,][,1:2]

###For TPM >10
blood_cortex_protein$bl_count_10 <- rowSums(blood3 >= 10)
blood_cortex_10 <- which(blood_cortex_protein$bl_count_10 >= 210)

cortex_protein$su_count_10 <- rowSums(brain3 >= 10)
cortex_10 <- which(cortex_protein$su_count_10 >= 215)
#View(blood_cortex_1)

blood_cortex_merged <- merge(blood_cortex_protein, cortex_protein, by = "Name")
blsu_count10 <- select(blood_cortex_merged, bl_count_10, su_count_10)
#View(blco_count1)

tpm10 <- blsu_count10[blsu_count10$bl_count_10 == 210 & blsu_count10$su_count_10 == 215,][,1:2]
tpm10 <- blsu_count10[blsu_count10$bl_count_10 == 210 & blsu_count10$su_count_10 < 215,][,1:2]
tpm10 <- blsu_count10[blsu_count10$bl_count_10 < 210 & blsu_count10$su_count_10 == 215,][,1:2]
