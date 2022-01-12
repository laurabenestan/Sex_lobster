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
library(vcfR)
library(dplyr)

### Load PCadapt outlier names
outliers<-unlist(read.table("02-outlier_detection/833snps_outliers.txt", sep='\t', header=F))

### Download sex data
sex <- read.table("00-data/256ind_palinurus_sex.txt", sep="\t")
colnames(sex) <- c("IND", "Sex")

### Load genetic data for the 833 outlier loci
vcf <- read.vcfR("02-outlier_detection/833snps_243ind_pcadapt.recode.vcf")
### Convert the vcf data in SNPs count matrix (0,1,2)
genind <- vcfR2genind(vcf)@tab
odd <- seq(1,1666,2)
snps <- cbind(rownames(genind), genind[,odd])
colnames(snps)[1] <- "IND"

write.table(snps, "07-Sex Environment/833snps_243ind_pcadapt.AlleleCount.txt", sep="\t")

### Match sex with SNP data
snps_sex <- merge(snps, sex, by="IND") %>%
  # remove individual with unknown sex
  filter(Sex %in% c("Female", "Male")) %>% 
  # make all columns numerical (except Sex and IND)
  mutate_at(2:834, as.numeric) %>% 
  # make sex factor
  mutate_at(835, as.factor)
  

### Do the logistic regression

## Prepare result matrices
results <- as.data.frame(matrix(0, 3, ncol(snps)))
rownames(results) <- c("Loci","z_val", "Pr(>|z|)")

### Run the regressions
for (i in 2:834) { # For each locus
  dat_i <- snps_sex[,c(i,835)]
  locus_i <- colnames(snps_sex)[i]
  colnames(dat_i) <- c("Y", "Sex")
  mod_i <- glm(Y ~ Sex, data=dat_i, family="poisson")
  
  #Fill result matrix
  results[1,(i-1)] <- locus_i
  results[2:3,(i-1)] <- summary(mod_i)$coefficients[2,c("z value", "Pr(>|z|)")]
  
  rm(dat_i)
  rm(mod_i)
}

write.table(results, file="07-Sex Environment/Logistic_Regression_SNPs_outliers_PCAadapt_results.txt", sep="\t")

### Loci correlated to Sex for each set of outliers (pval<0.01 for t test)
#### Threshold 0.01 
candidates <- results[,which(as.numeric(results[3,])<0.01)] %>%
  t(.) %>%
  as.data.frame(.) # 289 loci
write.table(candidates, file="05-sex-linked-markers/SNP_correlated_to_sex_289.txt", sep="\t", row.names=F)

cand_old <- read.table("05-sex-linked-markers/SNP_correlated_to_sex_278.txt", sep="\t")
setdiff(cand_old$V1, candidates$Loci) # 3
setdiff(candidates$Loci, cand_old$V1) # 15
intersect(cand_old$V1, candidates$Loci) # 274
