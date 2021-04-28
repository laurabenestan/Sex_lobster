' --------------------------------------------------------------------------   @Header
#'
#' @title Adaptive population genomic structure for Palinurus elephas according to gender
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
library(hierfstat)
library(mapdata)
library(viridis)
library(tidyr)
library(adegenet)
library(vcfR)
library(ggplot2)
library(poppr)
library(gridExtra)
library(dplyr)
library(cowplot)
library(viridis)
library(forcats)
library(colorBlindness)

################### FILTERING STEPS IN R ##################

### Download vcf
vcf <- vcfR::read.vcfR("../00-data/833snps_243ind_pcadapt.recode.vcf")

### Download pop info
pop <- read.table("../00-data/population_map_palinurus_243ind_mpa.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)

### Transform vcf to genind
genind <- vcfR::vcfR2genind(vcf)

### Define population map 
genind@pop <- pop$STRATA

################### POPULATION STRUCTURE IN R ##################

########### DPCA without a priori ######

## Use Bayesian Information Criterion on your filtered dataset.
grp_all <- find.clusters(genind, max.n.clust=5, n.pca=nInd(genind)/3) # select N/3
# Select K =2

### Observe the statistics
grp_all$Kstat

### Observe the size of the group
grp_all$size

### Perform the DPCA
dapc <- dapc(genind, grp_all$grp,n.pca=nInd(genind)/3, nda=1)
# Choose the 1 discriminant function to retain.

### Pick up nice colorblind palette
col2 <- c("#D55E00", "#009E73")

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

### Make a graph with all individuals
compoplot_ind <- dapc_compoplot_sex_ggplot %>%
  mutate(IND = fct_reorder(IND,value)) %>%
  ggplot() +
  geom_bar(aes(x = IND, y = value, fill = CLUSTER),stat="identity", width = 1, color="lightgrey")  +
  scale_fill_manual(values = col2) +
  guides(fill = "none") +
  labs(y = "Posterior membership probabilities") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank())+ facet_wrap(~SEX, scales = "free")

### Subset only male and female
dapc_compoplot_sex_final <- subset(dapc_compoplot_sex_ggplot, subset=dapc_compoplot_sex_ggplot$SEX=="Female"|dapc_compoplot_sex_ggplot$SEX=="Male")

### Make a graph
compoplot_sex <- dapc_compoplot_sex_final %>%
  mutate(IND = fct_reorder(IND,value)) %>%
  ggplot() +
  geom_bar(aes(x = IND, y = value, fill = CLUSTER),stat="identity", width = 1)  +
  scale_fill_manual(values = col2) +
  guides(fill = "none") +
  labs(y = "Posterior membership probabilities") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank())+ facet_wrap(~SEX, scales = "free")
compoplot_sex

### Save the graph
ggsave("Figure3.pdf", width=10, height=5)

### Extract DPCA information 
dapc$IND <- row.names(dapc$ind.coord)
dapc_info <- as.data.frame(cbind(dapc$IND, dapc$ind.coord, grp_all$grp))
colnames(dapc_info) <- c("IND","DPC1", "K")

### Merge the results with sex info
dapc_sex <- merge(dapc_info, sex, by="IND")

### Remove NA and Ind individuals
dapc_sex_all_info <- filter(dapc_sex, SEX!="Ind")

### Create geom_density
ggplot(dapc_sex_all_info, aes(x="", y=K, fill=K)) +
  geom_bar(stat="identity", width=1) +
  facet_wrap(~SEX)+
  coord_polar("y", start=0)+
  theme_bw()+
  scale_fill_manual(values=col2)+
  xlab("")
ggsave("Male_female_piechart_dapc.pdf")

################### ADMIXTURE IN R ##################

### Download libraries
library(stringr)
library(ggplot2)
library(dplyr)

### Download the **cross-validation** results 
cv <- read.table("cross_validation.txt")

### Analyze the **cross-validation** results
cv$K <-gsub("[\\(\\)]", "", regmatches(cv$V3, gregexpr("\\(.*?\\)", cv$V3))) 
CV <- select(cv, V4,K)
CV2 <- separate(data = CV, col = K, into = c("K", "CLUSTER"), sep = "=")

### Rename your two columns CV and K-cluster
colnames(CV2) <- c("CV","K","CLUSTER")

### Allow R to sort as numeric and not alphabetic.
CV2$CLUSTER <- as.numeric(as.character(CV2$CLUSTER))

###Do a **graph showing the cross validation results**.
graph_title="Cross-Validation plot" 
x_title="K"
y_title="Cross-validation error" 
cross_validation <- ggplot(CV2, aes(x=CLUSTER,y=CV,group=1))+
  geom_point()+ geom_line(linetype = "dashed")+ 
  scale_x_continuous(breaks = seq(0, 25, by = 1))+
  labs(title=graph_title)+ labs(x=x_title)+ labs(y=y_title)+ 
  theme_classic()
cross_validation

### Save the graph
ggsave("Admixture_cross-validation_adaptive.pdf",width=7,height=5,dpi=600)

#### Observe information get from ADMIXTURE output ####

#### Download libraries
library(reshape2)
library(plyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

### Read the file.Q with the Q = the optimal K
admixture <- read.table("833snps_243ind_pcadapt.2.Q")

### Add one column with the individuals names and the population they belong to.
admixture <- cbind(pop,admixture)

### Rename columns.
colnames(admixture) <- c("IND","POP","K1","K2")

### Transform the admixture object into a long format.
admixture_long <- reshape::melt(admixture,id.vars=c("IND","POP"),variable.name="ANCESTRY",value.name="PERC")
colnames(admixture_long) <- c("IND","POP","ANCESTRY","PERC")
class(admixture_long$ANCESTRY)
levels(admixture_long$ANCESTRY)

### Make a **graph with ADMIXTURE results**.
graph_title="Stacked barplot of Admixture analysis in species"
x_title="Individuals"
y_title="Ancestry"
mpa_ancestry<-ggplot(admixture_long,aes(x=POP,y=PERC,fill=ANCESTRY))
mpa_ancestry+geom_bar(stat="identity")+
  scale_fill_manual(values=col2, name= "K", labels=c("K1","K2"))+ 
  labs(y=y_title)+
  labs(x=x_title)+
  theme_classic()

### Estimate at which each genetic group belong each individual regarding its maximum % of ancestry.
admixture[, "max"] <- apply(admixture[, 3:4], 1, max)
head(admixture)

### Check the individuals that could not be clearly attributed to one genetic cluster.
admixture_subset <- subset(admixture, subset=admixture$max >= 0.7)

### Create a pop map regarding the cluster found.
admixture$CLUSTER <- colnames(admixture)[apply(admixture,1,which.max)]
admixture_results <- select(admixture, IND, CLUSTER, max)
head(admixture_results)

### Merge the results with sex info
admixture_sex <- merge(admixture, sex, by="IND")

### Create geom_density
ggplot(admixture_sex, aes(x="", y=max, fill=CLUSTER)) +
  geom_bar(stat="identity", width=1) +
  facet_wrap(~SEX)+
  coord_polar("y", start=0)+
  theme_bw()+
  scale_fill_manual(values=col2)+
  xlab("")
ggsave("Male_female_piechart_admixture.pdf")

### Save the files.
write.table(admixture_results, 'Admixture_results_K2.txt',quote=FALSE, row.names=FALSE, sep="\t", dec=".")

#### Produce a nice graph with pophelper ####

#### Download libraries
library(pophelper)
packageDescription("pophelper", fields="Version")

### Import labels for sampling sites
labset <- read.table("../00-data/population_map_palinurus_243ind_mpa.txt", sep="\t",header=TRUE,stringsAsFactors=F)
labset_order = labset[,2,drop=FALSE] # very important step for setting labels
labset_order$STRATA <- as.character(labset_order$STRATA)

### Check if labels are a character data type.
sapply(labset, is.character)
class(labset)

### Load admixture file
slist <- readQ("833snps_243ind_pcadapt.2.Q",filetype="basic")

### Sort the labels
labset_order$STRATA <-factor(labset_order$STRATA, levels = c("Cabo  de  Palos", 
                                                             "Cabo  de  Gata  Nijar", 
                                                             "Illa  de  Tabarca",
                                                             "Illes  Columbretes",
                                                             "Llevant  de  Mallorca",
                                                             "Norte  de  Menorca",
                                                             "Cap  de  Creus",
                                                             "Cerbere  Banyuls"))

### Check if labels are a character data type.
sapply(labset_order, is.character)
class(labset)
labset_order$STRATA <- as.character(labset_order$STRATA)

### Create a qplot for K = 2 considering two species.
slist1 <- alignK(slist[1]) 
admixture_plot <- plotQ(slist1,  clustercol= col2,grplab = labset_order, grplabsize=1.5,
      showsp=FALSE,ordergrp=T,imgtype="pdf",
      showlegend=TRUE, legendpos="right", legendkeysize = 4, legendtextsize = 4,
      legendmargin=c(2,2,2,0), width=20, height=4,sortind="all")

################### COMPARE DPCA AND ADMIXTURE IN R ##################
dapc_admixture <- merge(dapc_info, admixture_results, by="IND")
dapc_admixture$COMPARE <- substr(dapc_admixture$CLUSTER,2,3)
dapc_admixture$CONVERGENCE <- ifelse(dapc_admixture$COMPARE==dapc_admixture$K, "yes","no")

### Compare both analyses
dapc_admixture %>% group_by(CONVERGENCE) %>%
  tally()

################### COMPARE DPCA AND ADMIXTURE WITH SEX IN R ##################

### Merge sex info
sex_dapc_admixture <- merge(dapc_admixture, sex, by="IND")

### Compare both analyses
sex_dapc_admixture %>% group_by(SEX,K) %>%
  tally()
sex_dapc_admixture %>% group_by(SEX,CLUSTER) %>%
  tally()
