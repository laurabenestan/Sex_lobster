' --------------------------------------------------------------------------   @Header
#'
#' @title Sex-linked markers in Palinurus elephas according 
#'
#' @description
#' This R script
#'
#' @author   Alicia Dalongeville 
#'
#' @date 2021/04/21
#'
#' --------------------------------------------------------------------------   @libraries
### Download libraries
library(adegenet)
library(vcfR)
library(dplyr)
library(hierfstat)

### Download sex-linked markers info
sex_linked <- read.table("05-sex-linked-markers/SNP_correlated_to_sex_288.txt", header=T)

### Download sex info
sex <- read.table("00-data/243ind_sex.txt",header=TRUE,sep="\t")

### Create a vector of this list
vector_288snps <- as.vector(sex_linked$Loci)

### Keep only sex markers
ind_to_keep <- as.vector(sex$INDIVIDUALS)

### Read vcf
vcf <- read.vcfR("00-data/288snps_243ind_sexlinked.recode.vcf")

### Transform vcf to genind
genind_langouste <- vcfR2genind(vcf)

### Order the population strata
order <- as.data.frame(rownames(genind_langouste$tab))
colnames(order)="V1"
strata_order <- merge(order, sex, by.x="V1",by.y="INDIVIDUALS", sort=FALSE)

### Remove individuals with unknown genotype
genotyped <- sex %>%
  filter(STRATA %in% c("Male", "Female")) %>%
  pull(INDIVIDUALS)
genind2 <- genind_langouste[row.names(genind_langouste@tab) %in% genotyped]

sex_genotyped <- dplyr::filter(strata_order,V1 %in% genotyped)

### Specifying the sex of each individual as population in the genind
genind2@pop <- factor(sex_genotyped$STRATA)

### Transform genind to hierfstat
hierfstat_sex <- hierfstat::genind2hierfstat(genind2) 

### To avoid the error is that your levels variable (i.e. pop) is a factor, and because you removed whole levels without refactoring, it creates an error message. 
hierfstat_sex$pop<-factor(hierfstat_sex$pop)

### Calculate observed heterozygosity
basicstat <- hierfstat::basic.stats(hierfstat_sex, diploid = TRUE, digits = 2)
het_obs <- as.data.frame(basicstat$Ho)

### Add marker info
het_SL <- het_obs %>%
  mutate(Marker = rep("Sex-linked", nrow(het_obs)))

###########################################################
## Same for outliers non sex-linked markers
#import vcf file
vcf <- read.vcfR("00-data/545snps_243ind_outliers.recode.vcf")
genind_langouste <- vcfR2genind(vcf)
### Order the population strata
order <- as.data.frame(rownames(genind_langouste$tab))
colnames(order)="V1"
strata_order <- merge(order, sex, by.x="V1",by.y="INDIVIDUALS", sort=FALSE)
### Remove individuals with unknown genotype
genotyped <- sex %>%
  filter(STRATA %in% c("Male", "Female")) %>%
  pull(INDIVIDUALS)
genind2 <- genind_langouste[row.names(genind_langouste@tab) %in% genotyped]

sex_genotyped <- dplyr::filter(strata_order,V1 %in% genotyped)
### Specifying the sex of each individual as population in the genind
genind2@pop <- factor(sex_genotyped$STRATA)
### Transform genind to hierfstat
hierfstat_sex <- hierfstat::genind2hierfstat(genind2) 
hierfstat_sex$pop<-factor(hierfstat_sex$pop)
### Calculate observed heterozygosity
basicstat <- hierfstat::basic.stats(hierfstat_sex, diploid = TRUE, digits = 2)
het_obs <- as.data.frame(basicstat$Ho)
### Add marker info
het_out <- het_obs %>%
  mutate(Marker = rep("Outliers", nrow(het_obs)))

## Same for neutral markers
#import vcf file
vcf <- read.vcfR("/Users/alicia/Documents/MARBEC/CEFE_Langouste/00-Datasets/02-Neutral/25230snps_243ind.vcf")
genind_langouste <- vcfR2genind(vcf)
### Order the population strata
order <- as.data.frame(rownames(genind_langouste$tab))
colnames(order)="V1"
strata_order <- merge(order, sex, by.x="V1",by.y="INDIVIDUALS", sort=FALSE)
### Remove individuals with unknown genotype
genotyped <- sex %>%
  filter(STRATA %in% c("Male", "Female")) %>%
  pull(INDIVIDUALS)
genind2 <- genind_langouste[row.names(genind_langouste@tab) %in% genotyped]

sex_genotyped <- dplyr::filter(strata_order,V1 %in% genotyped)
### Specifying the sex of each individual as population in the genind
genind2@pop <- factor(sex_genotyped$STRATA)
### Transform genind to hierfstat
hierfstat_sex <- hierfstat::genind2hierfstat(genind2) 
hierfstat_sex$pop<-factor(hierfstat_sex$pop)
### Calculate observed heterozygosity
basicstat <- hierfstat::basic.stats(hierfstat_sex, diploid = TRUE, digits = 2)
het_obs <- as.data.frame(basicstat$Ho)
### Add marker info
het_neut <- het_obs %>%
  mutate(Marker = rep("Neutral", nrow(het_obs)))

###################################################
## Combine the 3 datasets together
ho_all <- rbind.data.frame(het_neut, het_SL, het_out)
## summarize results per marker type and sex
het_obs_loci <- ho_all %>%
  group_by(Marker) %>%
  summarise(across(everything(), list(mean = mean, min = min, max = max, sd = sd))) %>%
  as.data.frame(.)

write.csv(het_obs_loci, file="06-genetic_diversity/Hobs_per_loci.csv", row.names=F)

##### Plot
### Put in long format
ho_all_long <- reshape2::melt(ho_all)

### Pick up color-blind
col2 <- c("#D55E00", "#009E73")

# Plot
ho_all_long %>%
  ggplot(aes(x=variable, y=value, fill=variable)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_manual(values=col2,name="Gender") +
  theme_bw() +
  facet_wrap(~Marker)+
  ylab("Genetic diversity (Hobs)")+
  xlab("")+
  stat_compare_means(method= "wilcox.test",
                     label = "p.signif", 
                     symnum.args=list(
                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                       symbols = c("****", "***", "**", "*", "ns")))

ggsave("06-genetic_diversity/Hobs_loci.pdf",width=8, height=5)

##### GENETIC DIVERSITY WITH VCFTOOLS ########
### Download libraries
library(ggpubr)
het_all <- read.table("het_all.txt",sep="\t", header=TRUE,dec=",")

### Transform to long format
het_obs <- dplyr::select(het_all, INDV,Sexlinked, Adaptive, Neutral)
het_all_long <- tidyr::gather(het_obs, "MARKERS", "Hobs", Sexlinked:Neutral)

### Add sex info
hobs_sex <- merge(het_all_long,sex, by.x="INDV",by.y="INDIVIDUALS")

### Select only male and female 
female <- subset(hobs_sex, subset= STRATA=="Female")
male <- subset(hobs_sex,subset= STRATA=="Male")
female_male <- rbind(female, male)
female_male %>% group_by(STRATA,MARKERS) %>% dplyr::summarize(mean_hobs = mean(Hobs, na.rm = TRUE) , min = min(Hobs),max = max(Hobs),sd=sd(Hobs))

# sample size
sample_size = female_male %>% group_by(STRATA) %>% summarize(num=n())

### Pick up color-blind
col2 <- c("#D55E00", "#009E73")

# Plot
female_male %>%
 dplyr::left_join(sample_size) %>%
  ggplot(aes(x=STRATA, y=Hobs, fill=STRATA)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_manual(values=col2,name="Gender comparison") +
  theme_bw() +
  facet_wrap(~MARKERS)+
  ylab("Genetic diversity (Hobs)")+
  xlab("Gender comparison")+
  stat_compare_means(method= "wilcox.test")

ggsave("06-genetic_diversity/Figure5.pdf",width=8, height=5)


