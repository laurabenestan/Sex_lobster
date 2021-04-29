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

### Download sex-linked markers info
sex_linked <- read.table("SNP_correlated_to_sex_278.txt")

### Download sex info
sex <- read.table("../00-data/243ind_sex.txt",header=TRUE,sep="\t")

### Create a vector of this list
vector_278snps <- as.vector(sex_linked$V1)

### Keep only sex markers
ind_to_keep <- as.vector(sex$INDIVIDUALS)

### Read vcf
vcf <- read.vcfR("277snps-181ind.recode.vcf")

### Transform vcf to genind
genind_langouste <- vcfR2genind(vcf)

### Order the population strata
order <- read.table("277snps-181ind.tfam",header=FALSE)
strata_order <- merge(order, sex, by.x="V1",by.y="INDIVIDUALS", sort=FALSE)

### Specifying the sex of each individual
genind_langouste@pop <- factor(strata_order$STRATA)

### Add population strata
genotyped <- indNames(genind_langouste)
sex_genotyped <- dplyr::filter(strata_order,V1 %in% genotyped)
genind_langouste@pop <- factor(sex_genotyped$STRATA)

### Transform genind to hierfstat
hierfstat_sex <- hierfstat::genind2hierfstat(genind_langouste) 

### To avoid the error is that your levels variable (i.e. pop) is a factor, and because you removed whole levels without refactoring, it creates an error message. 
hierfstat_sex$pop<-factor(hierfstat_sex$pop)

### Calculate observed heterozygosity
basicstat <- hierfstat::basic.stats(hierfstat_sex, diploid = TRUE, digits = 2)
het_obs <- as.data.frame(basicstat$Ho)

### Put in long format
het_obs_long <- reshape2::melt(het_obs)

### Add color
col2 <- c("#88CCEE", "#CC6677")

### Distribution of Heterozygosity
het_density <- het_obs_long %>%
  ggplot(aes(x=value, group=variable, fill=variable)) +
  geom_density(adjust=1.5)+
  xlab("Heterozygosity observed")+
  scale_fill_manual(values= col2,name="Gender")+
  ylab("Number of SNP markers")+
  theme_classic()
het_density
ggsave("Heterozygosity_female_male.pdf")

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

ggsave("Figure4.pdf",width=8, height=5)

##### GENETIC DIVERSITY WITH ADEGENET ########

### Read vcf
vcf <- read.vcfR("../00-data-dryad/556snps-243ind-palinurus.recode.vcf")

### Transform vcf to genind
genind_langouste <- vcfR2genind(vcf)

#### Keep only male
male <- subset(sex, subset=sex$STRATA=="Male")
genind_male <- genind_langouste[i=male$INDIVIDUALS]
summary(genind_male)

#### Keep only female
female <- subset(sex, subset=sex$STRATA=="Female")
genind_female <- genind_langouste[i=female$INDIVIDUALS]
summary(genind_female)
