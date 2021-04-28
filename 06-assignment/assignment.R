' --------------------------------------------------------------------------   @Header
#'
#' @title Sexing Palinurus elephas according to sex-linked markers
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
library(tidyr)
library(vcfR)
library(ggplot2)
library(dplyr)
library(forcats)
library(colorBlindness)
library(ggpubr)

################### FILTERING STEPS IN R ##################

### Download vcf
vcf <- vcfR::read.vcfR("../00-data/277snps-sex-243ind-palinurus.recode.vcf")

### Download pop info
pop <- read.table("../00-data/population_map_palinurus_243ind_mpa.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)

### Transform vcf to genind
genind <- vcfR::vcfR2genind(vcf)

### Define population map 
genind@pop <- pop$STRATA

################### POPULATION STRUCTURE IN R ##################

########### DPCA without a priori ######

## Use Bayesian Information Criterion on your filtered dataset.
grp_all <- find.clusters(genind, max.n.clust=9, n.pca=nInd(genind)/3) # select N/3
# Select K =2

### Observe the statistics
grp_all$Kstat

### Observe the size of the group
grp_all$size

### Perform the DPCA
dapc <- dapc(genind, grp_all$grp,n.pca=nInd(genind)/3, nda=1)
# Choose the 1 discriminant function to retain.

### Pick up color-blind
safeColors
col2 <- c("#D55E00","#009E73")

### Check the result 
dpca_result <- scatter(dapc, col=col2)
dpca_result

### Visualize your results with complot
compoplot(dapc,col=col2,cleg=.6, posi=list(x=0,y=1.2), lab=loci.individuals.maf.het@pop)

### Save the posterior
dapc_compoplot <- as.data.frame(dapc$posterior)
dapc_compoplot$IND <- row.names(dapc_compoplot)

### Add sex info
sex <- read.table("../00-data/256ind_palinurus_sex.txt")
colnames(sex) <- c("IND","SEX")

### Combine dapc and sex info
dapc_compoplot_sex <- merge(dapc_compoplot, sex, by="IND")
dapc_compoplot_sex_ggplot <- gather(dapc_compoplot_sex, CLUSTER,value, `1`:`2`)

### Make a graph
compoplot_sex <- dapc_compoplot_sex_ggplot %>%
  mutate(IND = fct_reorder(IND,value)) %>%
  ggplot() +
  geom_bar(aes(x = IND, y = value, fill = CLUSTER),stat="identity", width = 1)  +
  scale_fill_manual(values = col2) +
  labs(y = "Posterior membership probabilities") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank())+ facet_wrap(~SEX, scales = "free")
compoplot_sex

### Save the graph
ggsave("Figure5.pdf", width=10, height=5)

### Extract cluster info 
assignment_groups <- as.data.frame(dapc$grp)
assignment_groups$INDV <- row.names(assignment_groups)
colnames(assignment_groups) <- c("CLUSTER","INDV")

### Subset only ind and NA
dapc_assignment_final <- subset(dapc_assignment, subset=dapc_assignment$SEX=="Ind"|is.na(dapc_compoplot_sex_ggplot$SEX))

### Combine dapc assignment and sex info
dapc_assignment <- merge(dapc_assignment_final, sex, by.x="INDV",by.y="IND")
write.table(dapc_assignment,"Assignment_results.txt", quote=FALSE,row.names=FALSE, sep="\t")

### CALCULATE in VCFTOOLS your heterozygosity values vcftools --vcf ../00-data/277snps-sex-243ind-palinurus.recode.vcf --het

### Download vcftools output
het_assignment <- read.table("out.het",sep="\t", header=TRUE,dec=",")

### Estimate Hobs
het_assignment$HOBS <- het_assignment$O.HOM./het_assignment$N_SITES

### Merge Hobs and clusters
het_assignment_cluster <- merge(dapc_assignment, het_assignment, by="INDV")

### Transform to long format
het_assignment_cluster_long <- select(het_assignment_cluster, CLUSTER,INDV,HOBS)

### Select only male and female 
female <- subset(het_assignment_cluster_long, subset= CLUSTER==2)
male <- subset(het_assignment_cluster_long,subset= CLUSTER==1)
female_male <- rbind(female, male)
het_assignment_cluster_long %>% group_by(CLUSTER) %>% summarize(mean_hobs = mean(HOBS, na.rm = TRUE) , min = min(HOBS),max = max(HOBS),sd=sd(HOBS))

# sample size
sample_size = female_male %>% group_by(CLUSTER) %>% summarize(num=n())

### Pick up color-blind
safeColors
col2 <- c("#D55E00","#009E73")

# Plot
het_assignment_cluster_long %>%
  left_join(sample_size) %>%
  mutate(CLUSTER = fct_reorder(CLUSTER, HOBS))   %>%
  ggplot(aes(x=CLUSTER, y=HOBS, fill=CLUSTER)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_manual(values=col2,name="Gender comparison") +
  theme_bw() +
  ylab("Genetic diversity (Hobs)")+
  xlab("Gender comparison")+
  stat_compare_means(method= "wilcox.test")

ggsave("Figure5b.pdf",width=5, height=5)

#### Calculate statistics
subset_males <- subset(het_assignment_cluster_long, het_assignment_cluster_long$HOBS>=0.90)
subset_female <- subset(het_assignment_cluster_long, het_assignment_cluster_long$HOBS<0.90)
