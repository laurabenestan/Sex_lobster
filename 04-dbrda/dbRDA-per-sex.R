###########################################################################################################
# Distance-based Redundancy Analysis (db-RDA)
#
# Author : Laura Benestan
# Date : 30-03-2020
# 
#Input:
#Euclidean distances or genepop file = response variable
#Environmental table = explanatory variables
#Order of the samples 

#Libraries that we will need
library(vcfR)
library(codep)
library(adespatial)
library(adegraphics)
library(vegan)
library(ape)
library(car)
library(adegenet)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(tidyr)
library(fastDummies)
library(ggplot2)
library(tidyverse)

##### Y = GENOMIC DATA ####

### Genomic data to be used

#import vcf file
vcf <- read.vcfR("00-data-dryad/500snps-241ind-nosex-outliers.recode.vcf")

#Transform vcf to genind file
data <- vcfR2genind(vcf)

# Environmental data
env <- read.table("00-data/env_geo_palinurus_241ind.txt", header=TRUE, sep="\t",dec=".")
hab <- read.table("00-data/sex_morpho_habitats_env_palinurus.txt", header=TRUE, sep="\t",dec=".")

# Sex data
sex <- read.table("00-data/241ind_sex.txt", sep="\t", header=T)
colnames(sex) <- c("IND", "Sex")

# create list of individuals (rows) to remove from @tab slot 
remove1 <- sex[which(sex$Sex %in% c("Male", "Female") == F), "IND"]
remove2 <- row.names(data@tab)[which(row.names(data@tab) %in% sex$IND == F)]
removeInd <- c(remove1,remove2)
# remove individuals from *genind object*
# note that in this step there's no longer a comma needed before the closing square bracket
genind <- data[!row.names(data@tab) %in% removeInd]
indiv <- row.names(genind@tab) # vector of the individuals in genetic dataset

## Combine sex and environment data frame
habitat <- env %>%
  inner_join(sex, by = c("labels"="IND"))

### Create one dataset for female and one for male
female <- habitat[which(habitat$Sex == "Female" ), "labels"]
male <- habitat[which(habitat$Sex == "Male" ), "labels"]
genind_M <- genind[!row.names(genind@tab) %in% female]
genind_F <- genind[!row.names(genind@tab) %in% male]

### Remove 20 males from the dataset at random to have the same number of individuals for both sexes
#removeM <- row.names(genind_M1@tab)[sample(seq(from = 1, to = 101, by = 1), size = 20, replace = FALSE)]
#genind_M <- genind_M1[!row.names(genind_M1@tab) %in% removeM]

#male <- row.names(genind_M@tab)

### Calculate Euclidean distances
distgenEUCL_M <- dist(genind_M, method = 
                        "euclidean", diag = FALSE, upper = FALSE, p = 2)

distgenEUCL_F <- dist(genind_F, method = 
                        "euclidean", diag = FALSE, upper = FALSE, p = 2)

##### X = PCA ENV ####

# Remove the variables that have a lot of NAs (all Salinity water columns)
EnvMatrix1 <- habitat[,c(1:4,8:10,14:25)]
EnvMatrix <- EnvMatrix1[complete.cases(EnvMatrix1),] %>%
  column_to_rownames(var="labels")

# Select the palinurus individuals
EnvMatrix_pal <- EnvMatrix[indiv,]

### Scale before PCA
sEnvMatrix = scale(EnvMatrix_pal, center = T, scale = T)

### We can use the prcomp() function to perform the principal component analysis.
ePCA = prcomp(sEnvMatrix)
pca_env_axis <- factoextra::fviz_eig(ePCA, addlabels = TRUE, ylim = c(0, 50), 
                                     barfill = "white",
                                     barcolor = "grey",
                                     linecolor = "black",
                                     main = " ")+ theme_classic()

### Represent correlation among explanatory variables
pdf("04-dbrda/FigS2.pdf")
corrplot::corrplot(ePCA$rotation[,1:3], is.corr=FALSE, 
                   method = "color", rect.col = "black", 
                   rect.lwd = 5,outline = F,cl.cex = 0.25, 
                   number.cex = 0.5,addCoef.col = "black",tl.col = "black",
                   col = colorRampPalette(c("darkred","white","midnightblue"))(100))
dev.off()

### Add color 
br_pal <- brewer.pal(8,"RdYlBu") 


### Add MPA info
mpa <- read.table("00-data/243ind/distance_to_mpa_palinurus_235ind.txt", sep="\t", header=TRUE)
mpa <- mpa[which(mpa$IND %in% indiv),] %>%
  arrange(match(IND, indiv))

### Represent PCA
pca_env_plot12 <- factoextra::fviz_pca_ind(ePCA,axes = c(1, 2),
                                           geom.ind = "point", # Montre les points seulement (mais pas le "text")
                                           col.ind = mpa$MPA, # colorer by groups
                                           palette = br_pal,
                                           addEllipses = TRUE, # Ellipses de concentration
                                           legend.title = "Marine reserve", title=" "
)
pca_env_plot23 <- factoextra::fviz_pca_ind(ePCA,axes = c(2, 3),
                                           geom.ind = "point", # Montre les points seulement (mais pas le "text")
                                           col.ind = mpa$MPA, # colorer by groups
                                           palette = br_pal,
                                           addEllipses = TRUE, # Ellipses de concentration
                                           legend.title = "Marine reserve",title=" "
)

cowplot::plot_grid(pca_env_axis, pca_env_plot12,pca_env_plot23, ncol=1)
ggsave("04-dbrda/FigS3_PCA.pdf", width=5, height=10)

# The principal components can be found in the $x matrix:
pca_axis <- as.data.frame(ePCA$x)

# Create one dataset for males and one for females
env_M <- pca_axis[rownames(genind_M@tab),]
env_F <- pca_axis[rownames(genind_F@tab),]

#### X = MPA DUMMIES ####

### Create dummy variables
results <- fastDummies::dummy_cols(mpa, select_columns = "MPA",remove_first_dummy = TRUE)
mpa_dummy <- results[7:ncol(results)]   
rownames(mpa_dummy) <- mpa$IND

# Create one dataset for males and one for females
mpa_dummy_M <- mpa_dummy[rownames(genind_M@tab),]
mpa_dummy_F <- mpa_dummy[rownames(genind_F@tab),]

#### X = HABITATS ####

### Select habitats per species
habitat_select <- hab[,c("Label","habitat_prediction","Depth", "Sex")]
colnames(habitat_select) <-c("IND","Habitats", "Depth", "Sex")

### Produce dummy variables
habitats_dummy <- fastDummies::dummy_cols(habitat_select, select_columns = "Habitats",remove_first_dummy = TRUE)
habitats_dummy <- habitats_dummy[,3:9]
row.names(habitats_dummy) <- habitat_select$IND
habitats_dummy <- habitats_dummy[which(rownames(habitats_dummy) %in% indiv),]

# Create one dataset for males and one for females
habitats_dummy_M <- habitats_dummy[rownames(genind_M@tab),]
habitats_dummy_F <- habitats_dummy[rownames(genind_F@tab),]

##### ORDINATION ON Y ####
Pcoa_M=pcoa(distgenEUCL_M)
Pcoa_M

Pcoa_F=pcoa(distgenEUCL_F)
Pcoa_F

### Choose the number of pcoa axis
Pcoa_M$values
Pcoa_F$values

### Extract Pcoa principal components, which will be the response variable in the db-RDA
X_M=as.data.frame(Pcoa_M$vectors)
X_F=as.data.frame(Pcoa_F$vectors)
X_M <- X_M[rownames(env_M),]
X_F <- X_F[rownames(env_F),]

#Look at genotypes distribution in relation to the first 2 Pcoa axes
plot(X_F[,1], X_F[,2])

###### Merge the datasets by row names
dataf <- merge(X_F, env_F,  by=0, all.x = T)
rownames(dataf) <- dataf$Row.names
X_F <- dataf[,2:65]
env_F <- dataf[,66:83]
datam <- merge(X_M, env_M,  by=0, all.x = T)
rownames(datam) <- datam$Row.names
X_M <- datam[,2:81]
env_M <- datam[,82:99]

############### SELECT ENV ####

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0envm<-rda(X_M ~ 1, env_M)
rda0envf<-rda(X_F ~ 1, env_F)

#OrdiR2step will move towards the global model with all explanatory variables
rdaGenvm<- rda(X_M ~ ., env_M)
rdaGenvf<- rda(X_F ~ ., env_F)

### Select the variables
Selenvm <- ordistep(rda0envm, scope=formula(rdaGenvm), direction="both") # no var selected
Selenvf <- ordistep(rda0envf, scope=formula(rdaGenvf), direction="both") # PC13

### Summary table with selected variables    
Selenvm$anova
Selenvf$anova

### Build a significant model
Ysel_envM=env_M[,1:3]
Ysel_envF=env_F[,1:3]
rdaSm<- rda(X_M ~ .,Ysel_envM)
rdaSf<- rda(X_F ~ .,Ysel_envF)

# Check the RDA summary
RsquareAdj(rdaSm)
RsquareAdj(rdaSf)
anova(rdaSm, perm=10000)
anova(rdaSf, perm=10000)

############### SELECT MEM ####
Coor=dplyr::select(mpa, LAT, LON)
Coorxy=dplyr::select(mpa, LON,LAT)

#look at sites spatial distribution 
plot(Coorxy, asp=1) 

###Compute spatial distances among sites accounting for the earth curvature
DistSpatial=gcd.hf(Coor) 

#Compute MEM
# Keeping default setting for truncation (length of the longest edge of the minimum spanning tree will be used as the threshold) and just positive MEM
dbmem = dbmem(DistSpatial) 

#Look at general output
summary(dbmem)

#Add the MEM to the dataset
row.names(dbmem) <- mpa$IND

#Create one dataset for males and one for females
dbmem_M <- as.data.frame(dbmem[which(rownames(dbmem) %in% male),])
dbmem_F <- as.data.frame(dbmem[which(rownames(dbmem) %in% female),])
dbmem_M <- dbmem_M[match(rownames(X_M), rownames(dbmem_M)), ]
dbmem_F <- dbmem_F[match(rownames(X_F), rownames(dbmem_F)), ]

#Visualizing the links longer than the threshold
adegraphics::s.label(Coor, nb = attr(dbmem, "listw"))

#Visualizing the mem. Can create more spatial structure type when MEM are selected together in the same analysis
#The 1rst MEMs are large spatial scales, the last MEMs are small spatial scales
ade4::s.value(Coorxy, dbmem[,1])

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0dbmem_f<-rda(X_F ~ 1, dbmem_F)
rda0dbmem_m<-rda(X_M ~ 1, dbmem_M)

#OrdiR2step will move towards the global model with all explanatory variables
rdaGdbmem_f<- rda(X_F ~ ., dbmem_F)
rdaGdbmem_m<- rda(X_M ~ ., dbmem_M)

### Select the variables
Seldbmemf <- ordistep(rda0dbmem_f, scope = formula(rdaGdbmem_f), direction="both") 
Seldbmemm <- ordistep(rda0dbmem_m, scope = formula(rdaGdbmem_m), direction="both") 

### Summary table with selected variables    
Seldbmemf$anova
Seldbmemm$anova

### Select the MEM to keep
mem_select_F <- dbmem_F %>% dplyr::select(MEM1)
mem_select_M <- dbmem_M %>% dplyr::select(MEM1)

### Check the model ith only dbmem
rdamem <- rda(X_F ~., mem_select_F)
summary(rdamem, scaling=1)  

# Check the RDA summary
RsquareAdj(rdamem)
anova(rdamem, perm=10000)

############### SELECT INSIDE/OUTSIDE ####

### Produce dummy variables
inside_dummy <- fastDummies::dummy_cols(mpa, select_columns = "CATEGORY",remove_first_dummy = TRUE)

# Add the inside variable to the dataset
row.names(inside_dummy) <- inside_dummy$IND
inside_dummy_F <- inside_dummy[which(rownames(inside_dummy) %in% rownames(env_F)),]
inside_dummy_M <- inside_dummy[which(rownames(inside_dummy) %in% rownames(env_M)),]

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda_category <-rda(X_F ~ inside_dummy_F$CATEGORY_OUTSIDE)

# Check the RDA summary
RsquareAdj(rda_category)
anova(rda_category, perm=10000)

############### CHECK SIGNIFICANCE OF THE MODEL ####

### Check the model with only dummy mpa
rda_mpa <- rda(X_F, inside_dummy_F$CATEGORY_OUTSIDE)
summary(rda_mpa, scaling=1)  

# Check the RDA summary
RsquareAdj(rda_mpa)
anova(rda_mpa, perm=10000)

##### FINAL RDA ####

### Build significant model
Yallsel <- cbind(pca_axis[,1:3], dbmem)
YallselM <- cbind(Ysel_envM, mem_select_M)  
YallselF <- cbind(Ysel_envF, mem_select_F)  

### Check correlation among variables
matrix_cor_M <- cor(YallselM)
cor(matrix_cor_M)
matrix_cor_F <- cor(YallselF)
cor(matrix_cor_F)

### Visualize the correlation matrix
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
plot_cor_M <- corrplot::corrplot(matrix_cor_M, method="color", col=col(200),  
                                 type="upper", 
                                 addCoef.col = "black", # Ajout du coefficient de corrélation
                                 tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
                                 # Cacher les coefficients de corrélation sur la diagonale
                                 diag=FALSE)

### Check each variable and keep the one with the highest % of explained variation
# for male

# MEM2 versus PC3 : MEM2
rdaMEM1<- rda(X_M, mem_select_M$MEM1)
RsquareAdj(rdaMEM1)
anova(rdaMEM1, perm=999)

rdaPC3<- rda(X_M, Ysel_envM[,1])
RsquareAdj(rdaPC3)
anova(rdaPC3, perm=999)


# for female

# MEM1 versus PC1 : MEM1
rdaMEM1<- rda(X_F, mem_select_F$MEM1)
RsquareAdj(rdaMEM1)
anova(rdaMEM1, perm=999)

rdaPC3<- rda(X_F, Ysel_envF[,1])
RsquareAdj(rdaPC3)
anova(rdaPC3, perm=999)


### Keep only non collinear variables
YfinalenvM <- dplyr::select(YallselM, PC2, PC3,MEM1)
YfinalenvF <- dplyr::select(YallselF, PC2,PC3, MEM1)

### Build the most significant model
rdaM<- rda(X_M ~., YfinalenvM)
rdaF<- rda(X_F ~.,YfinalenvF)

### Check the collinearity of the model
vif.cca(rdaF)
vif.cca(rdaM)

### Check the RDA summary
RsquareAdj(rdaF) # neutral : 0.0003974169 / out : 0.004625006
RsquareAdj(rdaM) # neutral : -4.414651e-05 / out : 0.003437396

### Test the significance
anova(rdaF, perm=999)
anova(rdaF, perm=999, by="margin")
anova(rdaF, perm=999, by="axis")
anova(rdaM, perm=999)
anova(rdaM, perm=999, by="margin")
anova(rdaM, perm=999, by="axis")

################# VISUALISE DB-RDA results ##################

###### MALE ####

### Getting the scores for plotting.
scrs <- scores(rdaM, display = c("species", "sites", "bp"), choices = 1:4, scaling = 2)
scrs

### Collect information on the pcao axes
species_centroids <- data.frame(scrs$species)
species_centroids
species_centroids$PC_names <- rownames(species_centroids)

### Add the PC components names.
species_centroids$PC_names <- rownames(species_centroids) 

### Add information of arrows
continuous_arrows <- data.frame(scrs$biplot)
continuous_arrows
rownames(continuous_arrows) <- c("PC2","PC3", "MEM1")

### Add names of variables
continuous_arrows$number <- rownames(continuous_arrows)
continuous_arrows
baseplot<-plot(rdaM, scaling = 2)
mult <- attributes(baseplot$biplot)$arrow.mul

### Add mpa names
sitenames<- mpa[which(mpa$IND %in% rownames(scrs$sites)), ] %>%
  arrange(match(IND, rownames(scrs$sites))) %>%
  pull(MPA)
sites_centroids <- data.frame(scrs$sites)
sites_centroids$SITE <- sitenames
head(sites_centroids)

### Make a ggplot
RDA_plotM <- ggplot(data = sites_centroids, aes(x = RDA1, y = RDA2))+
  geom_point(data = sites_centroids, pch = 21, size = 2, aes(fill = SITE))+
  scale_fill_manual(values = br_pal,name="Marine reserves")+
  geom_text(data = continuous_arrows,
            aes(x= (mult + mult/10) * RDA1, y = (mult + mult/10) * RDA2,
                label = number), size = 4, hjust = 0.5)+
  geom_segment(data = continuous_arrows,
               aes(x = 0, xend = mult * RDA1, y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, lty = "dotted") + geom_vline(xintercept = 0, lty = "dotted") +
  labs(x = paste("RDA1 (", round(summary(rdaM)$cont$importance[2,1]*100, 2), "%)", sep = ""), y = paste("RDA2 (", round(summary(rdaM)$cont$importance[2,2]*100, 2), "%)", sep = ""))+
  theme(axis.text = element_text(colour = "black", size = 12)) +
  theme(axis.title = element_text(size = 12, colour = "black", family = "Helvetica", face = "bold")) +
  theme_classic()
RDA_plotM

###### FEMALE ####

### Getting the scores for plotting.
scrs <- scores(rdaF, display = c("species", "sites", "bp"), choices = 1:4, scaling = 2)
scrs

### Collect information on the pcao axes
species_centroids <- data.frame(scrs$species)
species_centroids
species_centroids$PC_names <- rownames(species_centroids)

### Add the PC components names.
species_centroids$PC_names <- rownames(species_centroids) 

### Add information of arrows
continuous_arrows <- data.frame(scrs$biplot)
continuous_arrows
rownames(continuous_arrows) <- c("PC2","PC3", "MEM1")

### Add names of variables
continuous_arrows$number <- rownames(continuous_arrows)
continuous_arrows
baseplot<-plot(rdaF, scaling = 2)
mult <- attributes(baseplot$biplot)$arrow.mul

### Add mpa names
sitenames<- mpa[which(mpa$IND %in% rownames(scrs$sites)), ] %>%
  arrange(match(IND, rownames(scrs$sites))) %>%
  pull(MPA)
sites_centroids <- data.frame(scrs$sites)
sites_centroids$SITE <- sitenames
head(sites_centroids)

### Make a ggplot
RDA_plotF <- ggplot(data = sites_centroids, aes(x = RDA1, y = RDA2))+
  geom_point(data = sites_centroids, pch = 21, size = 2, aes(fill = SITE))+
  scale_fill_manual(values = br_pal, guide=FALSE)+
  geom_text(data = continuous_arrows,
            aes(x= (mult + mult/10) * RDA1, y = (mult + mult/10) * RDA2,
                label = number), size = 4, hjust = 0.5)+
  #  geom_text(data = sites_centroids, aes(x = RDA1, y = RDA2, label = SITE), 
  #            size = 2.5, col = "black", hjust = 1.2)+
  geom_segment(data = continuous_arrows,
               aes(x = 0, xend = mult * RDA1, y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, lty = "dotted") + geom_vline(xintercept = 0, lty = "dotted") +
  labs(x = paste("RDA1 (", round(summary(rdaF)$cont$importance[2,1]*100, 2), "%)", sep = ""), y = paste("RDA2 (", round(summary(rdaF)$cont$importance[2,2]*100, 2), "%)", sep = ""))+
  theme(axis.text = element_text(colour = "black", size = 12)) +
  theme(axis.title = element_text(size = 12, colour = "black", family = "Helvetica", face = "bold")) +
  theme_classic()
RDA_plotF

### Combine both RDA on male and female
pdf("04-dbrda/dbRDA_outliers.pdf",width = 10, height=5)
cowplot::plot_grid(RDA_plotF,RDA_plotM)
dev.off()
