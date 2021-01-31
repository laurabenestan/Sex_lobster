library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)

### Download kinship matrix
kinship <- read.table("kinship_palinurus.txt",header=TRUE)

### Transform matrix to long format
kinship_df <- melt(kinship)
colnames(kinship_df) <- c("INDV1","INDV2","KINSHIP")
sex <- read.table("256ind_palinurus_sex.txt",header=FALSE)
colnames(sex) <- c("IND","SEX")

### Merge sex and kinship
sex_indv <- merge(x=kinship_df, y=sex, by.x=c("INDV1"),by.y=c("IND"))
sex_indv <- merge(x=sex_indv, y=sex, by.x=c("INDV2"),by.y=c("IND"))

### Make categories
sex_indv$SEX_CAT <- ifelse(sex_indv$SEX.x=="Female" & sex_indv$SEX.y=="Female", "Female/Female", ifelse(sex_indv$SEX.x=="Male" &sex_indv$SEX.y=="Male","Male/Male","Female/Male"))
wilcox.test(NEUTRAL_HET ~ BUFFER5, data = serranus_morpho_fish_het_north)

### Find highly related individuals
sex_indv_no_self <- subset(sex_indv, subset=sex_indv$KINSHIP <=0.4)
sex_indv_no_self_no_na <-na.omit(sex_indv_no_self)

### Remove category female/male
sex_indv_no_self_no_na_subset <- subset(sex_indv_no_self_no_na,subset=sex_indv_no_self_no_na$SEX_CAT!="Female/Male")
montest=wilcox.test(KINSHIP~SEX_CAT, data=sex_indv_no_self_no_na_subset)

### Make boxplot on categories
sex_indv_no_self_no_na_subset %>%
  ggline(x = "SEX_CAT", y = "KINSHIP", 
         add = c("mean_se"), 
         ylab = "Kinship of Loiselle", xlab = "Comparison between genders")

### Save the graph
ggsave("Kinship_lobster_sex.pdf")

g2 <- het_all_adaptative %>%
  ggline(x = "SEX", y = "HET", 
         add = c("mean_se"), 
         ylab = "Heterozygosity observed", xlab = "Comparison between genders")
