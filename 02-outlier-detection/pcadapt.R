' --------------------------------------------------------------------------   @Header
#'
#' @title Outlier detection with pcadapt
#'
#' @description
#' This R script
#'
#' @author   Laura Benestan, \email{lmbenestan@gmail.com}
#'
#' @date 2021/01/28
#'
#' --------------------------------------------------------------------------   @libraries

### Download libraries
library(pcadapt)
library(vcfR)
library(ggplot2)
library(dplyr)
library(viridis)

### Download data 
data <- read.pcadapt("../00-data/83372snps_243ind.recode.vcf", type = "vcf")
ind <- read.table("../../../05-Landscape/00-Data/00-All/palinurus/00-all/83372snps_243ind.tfam")
ind_pop_mpa <- read.table("../../../05-Landscape/00-Data/04-Geo/population_map_palinurus_243ind_mpa.txt", header=TRUE, sep="\t")
#snp <- read.table("../../../05-Landscape/00-Data/00-All/palinurus/00-all/")
sex <- read.table("256ind_palinurus_sex.txt", sep="\t")

### Merge individuals and pop information
ind_mpa <- merge(ind, ind_pop_mpa, by.x=c("V1"), by.y=c("IND"))
ind_mpa_info <- select(ind_mpa, V1,LAT,LON,DISTANCE,MPA,CATEGORY)
#ind_sex <- merge(ind, sex, by=c("V1"))

### Run pcadapt with large K, for instance K = 20
data_pcadapt_trial <- pcadapt(data, K = 20, min.maf = 0.01) 

### Make screeplot and score plot to determine ideal K:
plot(data_pcadapt_trial, option = "screeplot", col, snp.info = NULL, plt.pkg = "ggplot")

### Reduce the number of PCs to see in order to better define the threshold
plot(data_pcadapt_trial, option = "screeplot", K = 10)
plot(data_pcadapt_trial, option = "scores", col, snp.info = NULL, plt.pkg = "ggplot")
#ideal K = 2

### Run pcadapt again with K = 2
data_pcadapt <- pcadapt(data, K = 2, min.maf = 0.01) 

### Create a dataframe gathering PCAadapt info and MPA info
pca_adapt_mpa <-cbind(data_pcadapt$scores, ind_pop_mpa) 

### Merge sex info for palinurus
sex_pcadapt <- merge(pca_adapt_mpa, sex, by.x="INDIVIDUALS",by.y=c("V1"))

### Visualising the results according to MPA and INSIDE OUTSIDE
ggplot(pca_adapt_mpa, aes(x=pca_adapt_mpa$`1`, y=pca_adapt_mpa$`2`, shape = CATEGORY, fill= factor(MPA)))+
         geom_point(size=1.5)+
  scale_shape_manual(values=c(21, 24))+
  scale_fill_brewer(palette="Accent", guide=FALSE)+
  facet_wrap(~MPA)+
  theme_classic()+
  xlab("PC1")+
  ylab("PC2")+ 
  theme_bw()

### Save the graph
ggsave("PCAdapt_MPA.pdf", width = 10, height=5)

### Visualising the results according to LON
ggplot(pca_adapt_mpa, aes(x=pca_adapt_mpa$`1`, y=pca_adapt_mpa$`2`, fill = pca_adapt_mpa$LON))+
  geom_point(pch=21, size=1.5)+
  scale_fill_viridis()+
  theme_classic()+
  xlab("PC1")+
  ylab("PC2")+ 
  theme(legend.title = element_blank()) 

### Visualising the results for sex info
library(viridis)
ggplot(sex_pcadapt, aes(x=sex_pcadapt$`1`, y=sex_pcadapt$`2`, fill = sex_pcadapt$V2))+
  geom_point(pch=21, size=1.5)+
  scale_fill_viridis(discrete=TRUE)+
  theme_classic()+
  xlab("PC1")+
  ylab("PC2")+ 
  theme(legend.title = element_blank()) 

### Save the graph
ggsave("PCAdapt_sex.pdf", width = 5, height=5)

### To get e.g. list of p-values:
snps_pvalues <- cbind(snp, data_pcadapt$pvalues)
snps_pvalues_no_na <- na.omit(snps_pvalues)
write.table(snps_pvalues, "All_Pvalues.txt", sep="\t", quote=FALSE)

### Graphical tools. The excess of small P-values indicate the presence of outliers.
plot(data_pcadapt, option = "qqplot", col, snp.info = NULL, plt.pkg = "ggplot")
plot(data_pcadapt, option = "manhattan", col, snp.info = NULL, plt.pkg = "ggplot")
hist(data_pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(data_pcadapt, option = "stat.distribution", col, snp.info = NULL, plt.pkg = "ggplot")

### Visualising the distribution of p-values
quantile(snps_pvalues_no_na$`data_pcadapt$pvalues`, probs = c(0.01, 0.99))

### Get onlt the top 1% of the markers
top_1percent <- subset(snps_pvalues_no_na, snps_pvalues_no_na$`data_pcadapt$pvalues` <= 4.823761e-39)
write.table(top_1percent, "Outliers.txt", sep="\t", quote=FALSE, row.names = FALSE)
