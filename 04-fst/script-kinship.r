' --------------------------------------------------------------------------   @Header
#'
#' @title Relatedness values for Palinurus elephas according to gender
#'
#' @description
#' This R script
#'
#' @author   Laura Benestan, \email{lmbenestan@gmail.com}
#'
#' @date 2021/01/28
#'
#' --------------------------------------------------------------------------   @libraries

### Dowload librairies
library(adegenet)
library(ade4)
library(ggplot2)
library(reshape)
library(stringr)
library(radiator)
library(mapdata)
library(viridis)
library(pegas)
library(dplyr)
library(hierfstat)
library(colorBlindness)

### Pick up nice colorblind palette
col2 <- c("#D55E00", "#009E73")

### Download genepop and sex info
vcf_pal <- vcfR::read.vcfR("../00-data-dryad/25230snps_243ind.vcf")
genind_pal <- vcfR::vcfR2genind(vcf_pal)
sex <- read.table("../00-data/256ind_palinurus_sex.txt", sep="\t", header=FALSE)

### Select only the individuals with sex info
ind <- indNames(genind_pal)
sex_ind <- sex[which(sex$V1 %in% ind),]
sex_genind <- genind_pal[i=sex_ind$V1]
sex_genind@pop <- as.factor(sex_ind$V2)
genind.imp <- poppr::missingno(sex_genind, type="mean")

## Convert genind object into hierfstat data
data_fstat <- genind2hierfstat(genind.imp)
data_fstat<-data_fstat[-230,] # remove the empty row at the end

# Sex biased test
test_res <- sexbias.test(data_fstat.imp,sex_ind$V2, nperm=10)
test_res$p.value

############ KINSHIP ####
### Download libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)

### Download kinship matrix
kinship <- read.table("kinship_palinurus.txt",header=TRUE)

### Transform matrix to long format
kinship_df <- melt(kinship)
colnames(kinship_df) <- c("INDV1","INDV2","KINSHIP")
sex <- read.table("../00-data/256ind_palinurus_sex.txt",header=FALSE)
colnames(sex) <- c("IND","SEX")

### Merge sex and kinship
sex_indv <- merge(x=kinship_df, y=sex, by.x=c("INDV1"),by.y=c("IND"))
sex_indv <- merge(x=sex_indv, y=sex, by.x=c("INDV2"),by.y=c("IND"))

### Make categories
sex_indv$SEX_CAT <- ifelse(sex_indv$SEX.x=="Female" & sex_indv$SEX.y=="Female", "Female/Female", ifelse(sex_indv$SEX.x=="Male" &sex_indv$SEX.y=="Male","Male/Male","Female/Male"))

### Find highly related individuals
sex_indv_no_self <- subset(sex_indv, subset=sex_indv$KINSHIP <=0.3)
sex_indv_no_self_no_na <-na.omit(sex_indv_no_self)

sex_indv_no_self_no_na %>%
ggplot( aes(x=SEX_CAT, y=KINSHIP, fill=SEX_CAT)) +
  geom_boxplot()

### Remove NA
sex_indv_no_na <- na.omit(sex_indv)
kruskal.test(KINSHIP ~ SEX_CAT, data = sex_indv_no_self_no_na)

### Remove category female/male
sex_indv_no_self_no_na_subset <- subset(sex_indv_no_self_no_na,subset=sex_indv_no_self_no_na$SEX_CAT!="Female/Male")
montest=kruskal.test(KINSHIP~SEX_CAT, data=sex_indv_no_self_no_na)
montest

# sample size
sample_size = sex_indv_no_self_no_na %>% group_by(SEX_CAT) %>% summarize(num=n())

# Plot
sex_indv_no_self_no_na_subset %>%
  left_join(sample_size) %>%
   ggplot( aes(x=SEX_CAT, y=KINSHIP, fill=SEX_CAT)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_manual(values=col2,name="Gender comparison") +
  theme_bw() +
  ylab("Kinship of Loiselle")+
  xlab("Gender comparison")

# avreage kinship per sex category
sex_indv_no_self_no_na %>% group_by(SEX_CAT) %>% summarize(mean=mean(KINSHIP),sd=sd(KINSHIP))


### Make boxplot on categories
sex_indv_no_self_no_na_subset %>%
  ggline(x = "SEX_CAT", y = "KINSHIP", 
         add = c("mean_se"), 
         ylab = "Kinship of Loiselle", xlab = "Comparison between genders")

### Save the graph
ggsave("Kinship_lobster_sex.pdf")
