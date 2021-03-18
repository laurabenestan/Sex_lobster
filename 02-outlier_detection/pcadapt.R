' --------------------------------------------------------------------------   @Header
#'
#' @title Outlier detection
#'
#' @description
#' This R script
#'
#' @author   Laura Benestan, \email{lmbenestan@gmail.com}
#'
#' @date 2021/03/15
#'
#' --------------------------------------------------------------------------   @libraries

### Download libraries
library(pcadapt)
library(vcfR)
library(ggplot2)
library(dplyr)
library(cowplot)
library(colorBlindness)

### Pick up color-blind
safeColors
col2 <- c("#D55E00", "#009E73")

### Download vcf data
data <- read.pcadapt("../00-data/83372snps_243ind.recode.vcf", type = "vcf")

### Download individual names 
ind <- read.table("../00-data/83372snps_243ind.tfam")

### Download marine reserve information
ind_pop_mpa <- read.table("../00-data/population_map_palinurus_243ind_mpa.txt", header=TRUE, sep="\t",dec=".")

### Import gender information
sex <- read.table("../00-data/256ind_palinurus_sex.txt", sep="\t")

### Merge individuals and pop information
ind_mpa_sex <- merge(ind_pop_mpa, sex, by.x="INDIVIDUALS",by.y="V1")
colnames(ind_mpa_sex) <- c("IND","MPA","SEX")

### Run pcadapt with large K, for instance K = 20
data_pcadapt_trial <- pcadapt(data, K = 20, min.maf = 0.01) 

### Make screenplot and score plot to determine ideal K:
plot(data_pcadapt_trial, option = "screeplot", col, snp.info = NULL, plt.pkg = "ggplot")

### Reduce the number of PCs to see in order to better define the threshold
plot(data_pcadapt_trial, option = "screeplot", K = 10)
plot(data_pcadapt_trial, option = "scores", col, snp.info = NULL, plt.pkg = "ggplot")
#ideal K = 2 as we can see two clusters

### Run pcadapt again with K = 2
data_pcadapt <- pcadapt(data, K = 2, min.maf = 0.01) 

### Create a dataframe gathering PCAadapt info and MPA info
pca_adapt_mpa <-cbind(data_pcadapt$scores, ind_mpa_info) 

### Merge sex info for palinurus
sex_pcadapt <- merge(pca_adapt_mpa, sex, by=("V1"))
colnames(sex_pcadapt) <- c("IND","PC1","PC2","MPA","SEX")

### Check the percent of vbariation explained by each axis
pca_axis <- data_pcadapt$singular.values

### Visualizing the results according to marine reserve
mpa <- ggplot(sex_pcadapt, aes(x=PC1, y=PC2,fill= factor(MPA)))+
         geom_point(size=1.5, shape=21)+
  scale_fill_brewer(palette="RdYlBu", name="Marine reserves")+
  theme_classic()+
  xlab(paste("PC1 (",round(pca_axis[1]*100,2), "%)"))+
  ylab(paste("PC2 (",round(pca_axis[2]*100,2), "%)"))+ 
  theme_bw()
mpa

### Find the two uncorrect genotypes
remove_ind <- subset(sex_pcadapt,subset=sex_pcadapt$PC2<= -0.5)

### Save the graph
ggsave("FigureS1.pdf", width = 5, height=5)

### Remove NA
sex_pcadapt2 <- na.omit(sex_pcadapt)

### Remove ind
sex_pcadapt3 <- subset(sex_pcadapt2, subset=sex_pcadapt2$SEX!="Ind")

### Visualizing the results according to gender info
sex <- ggplot(sex_pcadapt3, aes(x=PC1, y=PC2,fill= factor(SEX)))+
  geom_point(pch=21, size=1.5)+
  scale_fill_manual(values=col2,name="Gender information")+
  theme_bw()+
  xlab(paste("PC1 (",round(pca_axis[1]*100,2), "%)"))+
  ylab(paste("PC2 (",round(pca_axis[2]*100,2), "%)")) 
sex

### Save the graph
ggsave("FigureS2.pdf", width = 4, height=5)

### Combine figures
plot_grid(mpa, sex, nrow=2,ncol=1,labels=c("a","b"))

### Save the graph
ggsave("FigureS12.pdf", width = 4, height=8)

### To get e.g. list of p-values:
vcf <- read.vcfR("../00-data/83372snps_243ind.recode.vcf")
loci <- vcf@fix[,3]
snps_pvalues <- cbind(loci, data_pcadapt$pvalues)
snps_pvalues_no_na <- na.omit(as.data.frame(snps_pvalues))
colnames(snps_pvalues_no_na) <- c("SNP","P-value")
write.table(snps_pvalues, "All_Pvalues.txt", sep="\t", quote=FALSE)

### Graphical tools. The excess of small P-values indicate the presence of outliers.
plot(data_pcadapt, option = "qqplot", col, snp.info = NULL, plt.pkg = "ggplot")
plot(data_pcadapt, option = "manhattan", col, snp.info = NULL, plt.pkg = "ggplot")
hist(data_pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(data_pcadapt, option = "stat.distribution", col, snp.info = NULL, plt.pkg = "ggplot")

### Change class of P-value
snps_pvalues_no_na$Pvalue <- as.numeric(as.character(snps_pvalues_no_na$`P-value`))

### Visualizing the distribution of p-values
quantile(snps_pvalues_no_na$Pvalue, probs = c(0.01, 0.99))

### Get only the top 1% of the markers
top_1percent <- subset(snps_pvalues_no_na, snps_pvalues_no_na$Pvalue <= 4.823761e-39)
write.table(top_1percent, "Outliers.txt", sep="\t", quote=FALSE, row.names = FALSE)
