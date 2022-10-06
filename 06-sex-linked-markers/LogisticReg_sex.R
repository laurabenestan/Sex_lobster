setwd("~/Documents/CEFE Langouste")

### Dowload librairies
library(adegenet)
library(ade4)
library(ggplot2)
library(reshape)
library(stringr)
library(radiator)
library(mapdata)
library(viridis)
library(snpStats)
library(dplyr)
library(VariantAnnotation)

### Load PCadapt outlier names
outliers<-unlist(read.table("02-Adaptive_SNPs/01-PCA_adapt/833pcadapt_id_palinurus.txt", sep='\t', header=F))

### Download sex data
sex <- read.table("02-Sex/256ind_palinurus_sex.txt", sep="\t")
colnames(sex) <- c("IND", "Sex")

### Load genetic data for the 833 outlier loci
vcf <- read.vcfR("00-Datasets/03-Adaptive/833snps_243ind_pcadapt.recode.vcf")
### Convert the vcf data in SNPs count matrix (0,1,2)
genind <- vcfR2genind(vcf)@tab
odd <- seq(1,1666,2)
snps <- cbind(rownames(genind), genind[,odd])
colnames(snps)[1] <- "IND"

write.table(snps, "02-Sex/833snps_243ind_pcadapt.AlleleCount.txt", sep="\t")

### Match sex with SNP data
snps_sex <- merge(snps, sex, by="IND")

### Do the logistic regression

## Prepare result matrices
results <- matrix(0, 3, ncol(snps))
rownames(results) <- c("Loci","t_val", "Pr(>|t|)")

### Run the regressions
for (i in 2:ncol(snps)) { # For each locus
  y_i <- snps_sex[,i]
  mod_i <- lm(y_i ~ ., data=cbind.data.frame(y_i, snps_sex$Sex))
  
  #Fill result matrix
  results[1,i] <- colnames(snps)[i]
  results[2:3,i] <- summary(mod_i)$coefficients[2,c("t value", "Pr(>|t|)")]
  
  rm(x_i)
  rm(mod_i)
}

write.table(results, file="Logistic_Regression_SNPs_outliers_PCAadapt_results.txt", sep="\t")

### Loci correlated to Sex for each set of outliers (pval<0.01 for t test)
#### Threshold 0.01 
candidates <- results[,which(as.numeric(results[3,])<0.01)] # 278 loci
write.table(candidates, file="SNP_correlated_to_sex_278.txt", sep="\t")
