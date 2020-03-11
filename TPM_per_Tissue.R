#This script identifies the number of protein-coding genes within each brain region of interest and whole blood within the GTEx data set. It then goes to identify number of genes at certain expression levels (TPM >=1, >=5, >=10). Finally, it identifies the number of missed genes in the blood based on the brain patients that have expression in the blood (ie Blood1 = Cortex1). This also examines the missed genes at TPM >1, 5, and 10.

#set working directory
setwd("/Users/austinevanovich/Documents/GTEx_Blood_Brain")

#Read in data
cortex <- read.table("Cortex_counts.tsv", header = T)

substantia <- read.table("Substantia_counts.tsv", header = T)

spinalcord <- read.table("SpinalCord_counts.tsv", header = T)

putamen <- read.table("Putamen_counts.tsv", header = T)

nucleus <- read.table("Nucleus_counts.tsv", header = T)

hypothalamus <- read.table("Hypothalamus_counts.tsv", header = T)

hippocampus <- read.table("Hippocampus_counts.tsv", header = T)

frontalcortex <- read.table("FrontalCortex_counts.tsv", header = T)

cerebellum <- read.table("Cerebellum_counts.tsv", header = T)

cerebellarhemisphere <- read.table("CerebellarHemisphere_counts.tsv", header = T)

caudate <- read.table("Caudate_counts.tsv", header = T)

anterior <- read.table("Anterior_counts.tsv", header = T)

amygdala <- read.table("Amygdala_counts.tsv", header = T)

wholeblood <- read.table("WholeBlood_counts.tsv", header = T)

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
#Subsetting each count object to just the protein-coding genes by their ensembl gene ID: http://www.ensembl.org/biomart/martview/58e679a255c6e23f443473a2779f080e
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
brain2 <- wholeblood[,-1]
brain3 <- brain2[,-1]  #Use this value to determine number of patients in the brain tissue

blood2 <- blood_hippocampus[,-1]
blood3 <- blood2[,-1] #and this to determine number of patients in the brain-blood samples

#blood3 <- blood3[, -69]
#colnames(blood3)

blood_hippocampus$bl_count_1 <- rowSums(blood3 >= 1)
blood_hippocampus_1 <- which(blood_hippocampus$bl_count_1 >= 73)

wholeblood$su_count_1 <- rowSums(brain3 >= 1)
wholeblood <- which(wholeblood$su_count_1 >= 73)
#View(blood_hippocampus_1)

blood_hippocampus_merged <- merge(blood_hippocampus, hippocampus, by = "Name")
blsu_count1 <- select(blood_hippocampus_merged, bl_count_1, su_count_1)
#View(blco_count1)

tpm1 <- blsu_count1[blsu_count1$bl_count_1 == 73 & blsu_count1$su_count_1 == 73,][,1:2]
tpm1 <- blsu_count1[blsu_count1$bl_count_1 == 73 & blsu_count1$su_count_1 < 73,][,1:2]
tpm1 <- blsu_count1[blsu_count1$bl_count_1 < 73 & blsu_count1$su_count_1 == 73,][,1:2]
tpm1 <- blsu_count1[blsu_count1$bl_count_1 < 73 & blsu_count1$su_count_1 < 73,][,1:2]
#View(bcor2)

####For TPM >5
blood_hippocampus$bl_count_5 <- rowSums(blood3 >= 5)
blood_hippocampus_5 <- which(blood_hippocampus$bl_count_5 >= 73)

wholeblood$su_count_5 <- rowSums(brain3 >= 5)
wholeblood_5 <- which(wholeblood$su_count_5 >= 73)

blood_hippocampus_merged <- merge(blood_hippocampus, hippocampus, by = "Name")
blsu_count5 <- select(blood_hippocampus_merged, bl_count_5, su_count_5)
#View(blco_count1)

tpm5 <- blsu_count5[blsu_count5$bl_count_5 == 73 & blsu_count5$su_count_5 == 73,][,1:2]
tpm5 <- blsu_count5[blsu_count5$bl_count_5 == 73 & blsu_count5$su_count_5 < 73,][,1:2]
tpm5 <- blsu_count5[blsu_count5$bl_count_5 < 73 & blsu_count5$su_count_5 == 73,][,1:2]
tpm5 <- blsu_count5[blsu_count5$bl_count_5 < 93 & blsu_count5$su_count_5 < 93,][,1:2]

###For TPM >10
blood_hippocampus$bl_count_10 <- rowSums(blood3 >= 10)
blood_hippocampus_10 <- which(blood_hippocampus$bl_count_10 >= 73)

wholeblood$su_count_10 <- rowSums(brain3 >= 10)
wholeblood_10 <- which(wholeblood$su_count_10 >= 73)
#View(blood_hippocampus_1)

blood_hippocampus_merged <- merge(blood_hippocampus, hippocampus, by = "Name")
blsu_count10 <- select(blood_hippocampus_merged, bl_count_10, su_count_10)
#View(blood_cortex_merged)
#View(cortex_10)

tpm10 <- blsu_count10[blsu_count10$bl_count_10 == 73 & blsu_count10$su_count_10 == 73,][,1:2]
tpm10 <- blsu_count10[blsu_count10$bl_count_10 == 73 & blsu_count10$su_count_10 < 73,][,1:2]
tpm10 <- blsu_count10[blsu_count10$bl_count_10 < 73 & blsu_count10$su_count_10 == 73,][,1:2]
tpm10 <- blsu_count10[blsu_count10$bl_count_10 < 73 & blsu_count10$su_count_10 < 73,][,1:2]

#get the names of genes that are present in both, just brain, just blood, or neither
amygdala_genes_1 <- select(blood_amygdala_merged, Description.x, su_count_1, bl_count_1)
amygdala_genes_1_both <- subset(amygdala_genes_1, su_count_1 >= 59 & bl_count_1 >= 59)
amygdala_genes_1_blood <- subset(amygdala_genes_1, su_count_1 < 59 & bl_count_1 >= 59)
amygdala_genes_1_brain <- subset(amygdala_genes_1, su_count_1 >= 59 & bl_count_1 < 59)
amygdala_genes_1_neither <- subset(amygdala_genes_1, su_count_1 < 59 & bl_count_1 < 59)

View(amygdala_genes_1_both)

CD4 <- subset(blood_amygdala_merged, blood_amygdala_merged$Description.x == "CD4")

APOE_amyg <- subset(amygdala, amygdala$Description == "APOE")
APOE_blood <- subset(blood_amygdala, blood_amygdala$Description == "APOE")

CD4_amyg <- subset(amygdala, amygdala$Description == "CD4")
CD4_blood <- subset(blood_amygdala, blood_amygdala$Description == "CD4")

RPS20_amyg <- subset(amygdala, amygdala$Description == "RPS")
RPS20_blood <- subset(blood_amygdala, blood_amygdala$Description == "CD4")

#apoe2 <- subset(cortex_protein, cortex_protein$Description == "APOE")
geneselect <- function(DF, COL, GENE){
  results <- subset(DF, DF[[COL]] == GENE)
  return(results)
}

chooseRow <- function(DF, COL, VAL){
  results <- DF[DF$COL == VAL, ]
  return(results)
}

test2 <- amygdala[amygdala$su_count_1 == 59, ]
test2_b <- blood_amygdala[blood_amygdala$bl_count_1 == 59, ]
test3 <- chooseRow(DF = amygdala, COL = su_count_1, VAL = 59)

amygdala$su_count_1

#test <- geneselect(DF = blood_amygdala_merged, COL = "Description.x", GENE = "CNBD1")

APOE_anterior <- geneselect(DF = anterior, COL = "Description", GENE = "APOE")
APOE_blood <- geneselect(DF = blood_anterior, COL = "Description", GENE = "APOE")

CD4_anterior <- geneselect(DF = anterior, COL = "Description", GENE = "CD4")
CD4_blood <- geneselect(DF = blood_anterior, COL = "Description", GENE = "CD4")

RPS20_anterior <- geneselect(DF = anterior, COL = "Description", GENE = "RPS20")
RPS20_blood <- geneselect(DF = blood_anterior, COL = "Description", GENE = "RPS20")


#View(CD4)
#View(test)
#View(apoe2)

#col.names.remove <- c("su_count_1", "su_count_5","su_count_10", "Name", #"Description")
#
#col.names.remove.bl <- c("bl_count_1", "bl_count_5","bl_count_10", "Name", #"Description")
#
#RPS20_blood <- RPS20_blood[,!(colnames(RPS20_blood) %in% col.names.remove.bl)]


#start to fit a linear model for TPM levels

Drop_Cols <- function(DF, COL){
  res <- DF[, -COL]
  return(res)
}

t1 <- colnames(amygdala)
test <- amygdala[, -which(t1 %in% c("Name", "Description", "su_count_1", "su_count_5", "su_count_10"))]

colnames(test)

drop_brain_APOE <- Drop_Cols(DF = APOE_anterior, COL = c(1,2,71,72,73))
drop_blood_APOE <- Drop_Cols(DF = APOE_blood, COL = c(1,2,71,72,73))

drop_brain_CD4 <- Drop_Cols(DF = CD4_anterior, COL = c(1,2,71,72,73))
drop_blood_CD4 <- Drop_Cols(DF = CD4_blood, COL = c(1,2,71,72,73))

drop_brain_RPS20 <- Drop_Cols(DF = RPS20_anterior, COL = c(1,2,71,72,73))
drop_blood_RPS20 <- Drop_Cols(DF = RPS20_blood, COL = c(1,2,71,72,73))


transpose <- function(DF, COL){
  res <- as.data.frame(t(DF))
  colnames(res) = c(COL)
  return(res)
}

APOE_brain_turn <- transpose(DF = drop_brain_APOE, COL = "TPM")
APOE_blood_turn <- transpose(DF = drop_blood_APOE, COL = "TPM")

CD4_brain_turn <- transpose(DF = drop_brain_CD4, COL = "TPM")
CD4_blood_turn <- transpose(DF = drop_blood_CD4, COL = "TPM")

RPS20_brain_turn <- transpose(DF = drop_brain_RPS20, COL = "TPM")
RPS20_blood_turn <- transpose(DF = drop_blood_RPS20, COL = "TPM")

#View(CD4_blood_turn)


#APOE_model <- lm(APOE_at$`48985` ~ APOE_bt$`23729`)
APOE_model <- lm(APOE_brain_turn$TPM ~ APOE_blood_turn$TPM)
summary(APOE_model)


plot(APOE_blood_turn$TPM, APOE_brain_turn$TPM,  xlab = "TPM in Blood", ylab = "TPM in anterior", pch = c(16), col = c("blue"), main = "APOE")
abline(APOE_model)

#================================

CD4_model <- lm(CD4_brain_turn$TPM ~ CD4_blood_turn$TPM)
summary(CD4_model)


plot(CD4_blood_turn$TPM, CD4_brain_turn$TPM,  xlab = "TPM in Blood", ylab = "TPM in anterior", pch = c(16), col = c("blue"), main = "CD4")
abline(CD4_model)

#================================

RPS20_model <- lm(RPS20_brain_turn$TPM ~ RPS20_blood_turn$TPM)
summary(RPS20_model)


plot(RPS20_blood_turn$TPM, RPS20_brain_turn$TPM, xlab = "TPM in Blood", ylab = "TPM in anterior", pch = c(16), col = c("blue"), main = "RPS20")
abline(RPS20_model)

#legend("topright", legend = c("Amygdala", "Blood"), col = c("blue", "red"),
#pch = c(16,18))
plot(amygdala_model)

amygdala_model_CD4 <- lm(CD4_at$`32653` ~ CD4_bt$`32653`)
summary(amygdala_model_rps20)

plot(RPS20_bt$`23729`, RPS20_at$`23729`, xlab = "TPM in Blood", ylab = "TPM in Amygdala", pch = c(16), col = c("blue"), main = "RPS20")
#legend("topright", legend = c("Amygdala"), col = c("blue"),
      # pch = c(16))
abline(amygdala_model_rps20)
#test <- amygdala_protein %>% filter(
#  Description == "RPS20")

APOE_amyg <- APOE_amyg[, -1]
View(APOE_amyg)
APOE_at <- as.data.frame(t(APOE_amyg), colnames = c("ID", "TPM"))

APOE_blood <- APOE_blood[, -1]
APOE_bt <- as.data.frame(t(APOE_blood), colnames = c("ID", "TPM"))
View(APOE_at)
write.table(RPS20bl, file = "blood_amygdala_rps20.csv", sep = ",")

RPS20f <- cbind(RPS20bl, RPS20t)

RPS20t <- as.matrix(RPS20t)
colnames(RPS20t) <- c( "TPM")
RPS20bl <- RPS20bl[-1, ]
View(RPS20t)

colnames(RPS20f) <- c("TPM", "ID")

library(stringr)
library(broom)

RPS20blood <- RPS20bl[-c(1,2),]
RPS20amyg <- RPS20t[-c(1,2),]
colnames(RPS20amyg) <- c("test","second")
View(RPS20amyg)
RPS20_model <- merge(RPS20blood, RPS20amyg)

row.names.remove <- c("su_count_1", "su_count_5","su_count_10")

RPS20amyg <- RPS20amyg[!(row.names(RPS20amyg) %in% row.names.remove), ]

library(tidyverse)

RPS20amyg <- RPS20amyg[grepl('GTEX', RPS20amyg),]

RPS20amyg <- data.frame(stringsAsFactors = FALSE,
                            EA = c("Los Angeles, CA", "Other text")
)

large_dataset %>% 
  filter(str_detect(EA, pattern = "CA"))