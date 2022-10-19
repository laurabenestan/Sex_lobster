' --------------------------------------------------------------------------   @Header
#'
#' @title Sex-linked markers in Palinurus elephas according 
#'
#' @description
#' This R script
#'
#' @author  Laura Benestan
#'
#' @date 2022/10/18
#'
#' --------------------------------------------------------------------------   @libraries
### Download libraries
library(adegenet)
library(vcfR)
library(dplyr)
library(ggplot2)
library(ggpubr)

### Download sex-linked markers info
sex_linked <- read.table("SNP_correlated_to_sex_331.txt", header=TRUE)

### Download sex info
sex <- read.table("../00-data/256ind_palinurus_sex.txt",header=FALSE,sep="\t")
colnames(sex) <- c("INDV","STRATA")

### Download heterozygosity values
het_all <- read.table("het-all-241ind.txt",sep="\t", header=TRUE,dec=",")

### Transform to long format
het_obs <- dplyr::select(het_all, INDV,SEX, OUTLIERS, NEUTRAL)
het_all_long <- tidyr::gather(het_obs, "MARKERS", "Hobs", SEX:NEUTRAL)

### Add sex info
hobs_sex <- merge(het_all_long,sex, by="INDV")

### Select only male and female 
female <- subset(hobs_sex, subset= STRATA=="Female")
male <- subset(hobs_sex,subset= STRATA=="Male")
female_male <- rbind(female, male)

### Make a summary of the heterozygosity values
summary_het <- female_male %>% 
  group_by(STRATA,MARKERS) %>% 
  dplyr::summarize(mean_hobs = mean(Hobs, na.rm = TRUE) , min = min(Hobs),max = max(Hobs),sd=sd(Hobs))
summary_het 

### Check group sample size
sample_size <- female_male %>% 
  group_by(STRATA) %>% 
  summarize(num=n())
sample_size

### Pick up color-blind
col2 <- c("#D55E00", "#009E73")

### Make a nice graph
female_male %>%
# dplyr::left_join(sample_size) %>%
  ggplot(aes(x=STRATA, y=Hobs, fill=STRATA)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_manual(values=col2,name="Gender comparison") +
  theme_bw() +
  facet_wrap(~MARKERS)+
  ylab("Genetic diversity (Hobs)")+
  xlab("Gender comparison")+
  stat_compare_means(method= "wilcox.test")

ggsave("Figure4-241ind.pdf",width=8, height=5)
