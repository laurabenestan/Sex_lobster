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

################### FILTERING STEPS IN R ##################

### Download vcf
vcf <- vcfR::read.vcfR("../../../00-Datasets/03-Adaptive/833snps_243ind_pcadapt.recode.vcf")

### Download pop info
pop <- read.table("../../../00-Datasets/04-Environmental_data/population_map_palinurus_243ind_mpa.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)

### Transform vcf to genind
genind <- vcfR::vcfR2genind(vcf)

### Define population map 
genind@pop <- pop$STRATA

################### POPULATION STRUCTURE IN R ##################

########### DPCA without a priori ######

## Use Bayesian Information Criterion on your filtered dataset.
grp_all <- find.clusters(genind, max.n.clust=10, n.pca=nInd(genind)/3) # select N/3
# Select K =2

### Observe the statistics
grp_all$Kstat

### Observe the size of the group
grp_all$size

### Perform the DPCA
dapc <- dapc(genind, grp_all$grp,n.pca=nInd(genind)/3, nda=1)
# Choose the 1 discriminant function to retain.

### Give a colorblind palet.
col2 <- c("#88CCEE", "#CC6677")

### Check the result 
dpca_result <- scatter(dapc, col=col2)
dpca_result

### Visualize your results with complot
compoplot(dapc,col=col2,cleg=.6, posi=list(x=0,y=1.2), lab=loci.individuals.maf.het@pop)

### Save the posterior
dapc_compoplot <- as.data.frame(dapc$posterior)
dapc_compoplot$IND <- row.names(dapc_compoplot)

### Combine dapc and sex info
dapc_compoplot_sex <- merge(dapc_compoplot, dapc_sex_all_info, by="IND")
dapc_compoplot_sex_ggplot <- gather(dapc_compoplot_sex, CLUSTER,value, `1`:`2`)

### Make a graph
compoplot_ggplot <- ggplot(dapc_compoplot_sex_ggplot) +
  geom_bar(aes(x = IND, y = value, fill = CLUSTER),stat="identity", width = 1)  +
  scale_fill_manual(values = col2) +
  guides(fill = "none") +
  labs(y = "Posterior membership probabilities") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank())+ facet_wrap(~GENDER, scales = "free")

### Extract DPCA information 
dapc$IND <- row.names(dapc$ind.coord)
dapc_info <- as.data.frame(cbind(dapc$IND, dapc$ind.coord, grp_all$grp))
colnames(dapc_info) <- c("IND","DPC1", "K")

### Add a site info by keeping only the three letters of the individual names
dapc_info$SITE <- substr(compoplot_ggplot, 1,3)
dapc_info 

### Merge the results with sex info
dapc_sex <- merge(dapc_info, sex, by="IND")

### Remove NA and Ind individuals
dapc_sex_all_info <- filter(dapc_sex, GENDER!="Ind")

### Create geom_density
piechart_dapc <- ggplot(compoplot_ggplot, aes(x="", y=K, fill=K)) +
  geom_bar(stat="identity", width=1) +
  facet_wrap(~GENDER)+
  coord_polar("y", start=0)+
  theme_bw()+
  scale_fill_manual(values=col2)+
  xlab("")
ggsave("Male_female_piechart.pdf")

################### VISUALIZE RESULTS ##################
pdf("Fig2.pdf",width=10, height=10)
plot_grid(compoplot_ggplot, piechart_dapc, nrow=2,labels=c("A", "B"))
dev.off()

#################### Create a map #####
### Import geographic coordinates file
sites <-read.table("env_geo_palinurus_243ind.txt",header=TRUE, dec=".",sep="\t",na.strings="NA",strip.white=T)
summary(sites)

### Download sex info
sex <- read.table("../../../02-Sex/256ind_palinurus_sex.txt")
colnames(sex) <- c("IND","GENDER")

### Merge sites and sex
sites_sex <- merge(sites, sex, by="IND")
length(which(is.na(sites_sex$GENDER)))
sites_sex %>% group_by(GENDER) %>%
  tally()

### Remove Na and Ind
target <- c("Female", "Male")
sites_sex2 <- sites_sex %>% filter(GENDER %in% target)

### Count the number of samples per lattitude and longitude points
sites_number <- sites_sex2 %>% group_by(LAT,LON,GENDER) %>%
  tally()

# Download the map for the Mediterranean Sea
wH <- map_data("worldHires",  xlim=c(-8,37), ylim=c(29.5,47)) # subset polygons surrounding med sea

### Map
x_title="Longitude"
y_title="Latitude"
map_lobster <- ggplot() +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = NA) +
  coord_fixed(xlim=c(-8,8), ylim=c(35.5,46.5), ratio=1.2)+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())+
  geom_point(aes(x = LON, y = LAT,size=n, color=GENDER), data=sites_number,shape = 19,alpha = 0.5)+
  xlab("Longitude") + ylab("Latitude") +
  theme_bw()+ 
  scale_color_manual(values=col2, name="Gender")+
  facet_wrap(~GENDER)+
  theme(axis.text.x=element_text(colour="black",size=14))+
  theme(axis.text.y=element_text(colour="black", size=14))+
  theme(axis.title=element_text(colour="black",size=14))
map_lobster

### Save the map
ggsave("Map_lobster.pdf",width=10, height=10)

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
  scale_x_continuous(breaks = seq(0, 10, by = 1))+
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
admixture_long <- melt(admixture,id.vars=c("IND","POP"),variable.name="ANCESTRY",value.name="PERC")
names(admixture_long)
class(admixture_long$ANCESTRY)
levels(admixture_long$ANCESTRY)

### Subset only the individuals with about 50% of ancestry with one genetic cluster
admixture_50 <- subset(admixture, subset=admixture$K1>0.40&admixture$K1<0.60)
dim(admixture_50)

### Keep only individuals assigned to more than 10% to a genetic cluster
admixture_long_info <- subset(admixture_long, subset=admixture_long$PERC >=0.10)

### Make a **graph with ADMIXTURE results**.
graph_title="Stacked barplot of Admixture analysis in species"
x_title="Individuals"
y_title="Ancestry"
mpa_ancestry<-ggplot(admixture_long_info,aes(x=POP,y=PERC,fill=ANCESTRY))
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
ggplot(admixture_sex, aes(x="", y=CLUSTER, fill=CLUSTER)) +
  geom_bar(stat="identity", width=1) +
  facet_wrap(~GENDER)+
  coord_polar("y", start=0)+
  theme_bw()+
  scale_fill_manual(values=col2)+
  xlab("")
ggsave("Male_female_piechart.pdf")

### Save the files.
write.table(admixture_results, 'Admixture_results_K2.txt',quote=FALSE, row.names=FALSE, sep="\t", dec=".")

#### Produce a nice graph with pophelper ####

#### Download libraries
library(pophelper)
packageDescription("pophelper", fields="Version")

### Import labels for sampling sites
labset <- read.table("../../00-Datasets/04-Environmental_data/population_map_palinurus_243ind_mpa.txt", header=TRUE,stringsAsFactors=F)
labset$REGIONS <- sites$LAT
labset_order = labset[,2,drop=FALSE] # very important step for setting labels
labset_order$STRATA <- as.character(labset_order$STRATA)

### Check if labels are a character data type.
sapply(labset, is.character)
class(labset)

### Load admixture file
slist <- readQ("25230snps_243ind.2.Q",filetype="basic")

### Create a qplot for K = 2 considering two species.
slist1 <- alignK(slist[1]) 
admixture_plot <- plotQ(slist1,  clustercol= col2,grplab = labset_order, grplabsize=1.5,
      showsp=FALSE,ordergrp=T,imgtype="pdf",selgrp="V3",
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
sex_dapc_admixture %>% group_by(GENDER,K) %>%
  tally()
sex_dapc_admixture %>% group_by(GENDER,CLUSTER) %>%
  tally()

################### VISUALIZE RESULTS ##################
plot_grid(map_lobster, dpca_result)
