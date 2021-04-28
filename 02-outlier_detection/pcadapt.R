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
data <- read.pcadapt("../00-data-dryad/83372snps_243ind.recode.vcf", type = "vcf")
ind <- read.table("../00-data-dryad/83372snps_243ind.tfam")
ind_pop_mpa <- read.table("../00-data/population_map_palinurus_243ind_mpa.txt", header=TRUE, sep="\t")
sex <- read.table("../00-data/256ind_palinurus_sex.txt", sep="\t")
snp <- as.data.frame(vcfR::getID(read.vcfR("../00-data-dryad/83372snps_243ind.recode.vcf")))
  
### Merge individuals and pop information
ind_mpa <- merge(ind, ind_pop_mpa, by.x=c("V1"), by.y=c("INDIVIDUALS"))

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
colnames(sex_pcadapt) <- c("INDIVIDUALS","PC1","PC2","MPA","SEX")

### Visualizing the results according to MPA and INSIDE OUTSIDE
g1 <- sex_pcadapt %>%
ggplot(aes(x=PC1, y=PC2, fill= MPA))+
         geom_point(size=2, shape=21)+
  scale_fill_brewer(palette="Accent")+
  theme_classic()+
  xlab("PC1")+
  ylab("PC2")

### Save the graph
ggsave("PCAdapt_MPA.pdf", width = 10, height=5)

### Pick up nice colorblind palette
col2 <- c("#D55E00", "#009E73")

### Subset only individuals with known sex
subset_sex <- subset(sex_pcadapt,(sex_pcadapt$SEX=="Male"|sex_pcadapt$SEX=="Female"))

### Visualizing the results according to SEX
g2 <- subset_sex %>%
  ggplot(aes(x=PC1, y=PC2, fill= SEX))+
  geom_point(size=2, shape=21)+
  scale_fill_manual(values=col2)+
  theme_classic()+
  xlab("PC1")+
  ylab("PC2")

### Combine the graphs
cowplot::plot_grid(g1,g2)

### Save the graph
ggsave("FigureS1.pdf", width = 5, height=5)

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
top_1percent <- subset(snps_pvalues_no_na, snps_pvalues_no_na$`data_pcadapt$pvalues` <= 6.184358e-29)
write.table(top_1percent, "833snps_outliers.txt", sep="\t", quote=FALSE, row.names = FALSE)
