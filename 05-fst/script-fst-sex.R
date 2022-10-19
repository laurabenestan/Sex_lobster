' --------------------------------------------------------------------------   @Header
#'
#' @title Genetic differentiation for Palinurus elephas according to gender
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
library(assigner)
library(radiator)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

### Pick up nice colorblind palette
col2 <- c("#D55E00", "#009E73")

###### STEP 1 : CALCULATE FST ON BOTH SEXES ####

#### Download data for diplodus
tidy_female <- tidy_genepop("../00-data-dryad/25230snps-79ind.txt", strata="../00-data/population_map_palinurus_241ind_mpa.txt")

#### Calculate fst on this vcf data
fst_female <- fst_WC84(
  data = tidy_female,
  holdout.samples = NULL, 
  pairwise = TRUE,
  ci = TRUE,
  iteration.ci = 10,
  quantiles.ci = c(0.025, 0.975),
  digits = 9,
  parallel.core = 8
)
write.table(fst_mullus$pairwise.fst,"Fst_female.txt",sep="\t", row.names = FALSE, quote=FALSE)

### Extract the FST matrix
fst.matrix <- fst_female$pairwise.fst.full.matrix

### Remove the _ character
fst_matrix <- gsub("_", " ", fst.matrix)

### Melt the FST matrix to produce the nice graph
melted_fst <- reshape::melt(fst_matrix,na.rm = FALSE)

### Verify that the Fst values is indeed a numerical value
melted_fst$value <- as.numeric(as.character(melted_fst$value))

### Remove NA
melted_fst <- na.omit(melted_fst)
melted_fst$value <- as.numeric(as.character(melted_fst$value))

### Make the graph
female <-ggplot2::ggplot(data = melted_fst, aes(X2, X1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "yellow", mid ="white", high = "#D55E00", space = "Lab", 
                       name="FST") +
  ylab("Sampling location A")+
  xlab("Sampling location B")+
  theme_bw()+ 
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 45, vjust = 1, 
                                                                size = 12, hjust = 1))+  
  coord_fixed()
female

### Save the graph
ggsave("fst_female.pdf", width=10, height=10)

#### Download data for diplodus
tidy_male <- tidy_genepop("../00-data-dryad/25230snps-101ind.txt", strata="../00-data/population_map_palinurus_241ind_mpa.txt")

#### Calculate fst on this vcf data
fst_male <- fst_WC84(
  data = tidy_male,
  holdout.samples = NULL, 
  pairwise = TRUE,
  ci = TRUE,
  iteration.ci = 10,
  quantiles.ci = c(0.025, 0.975),
  digits = 9,
  parallel.core = 8
)

### Extract the FST matrix
fst.male.matrix <- fst_male$pairwise.fst.full.matrix

### Remove the _ character
fst.male.matrix <- gsub("_", " ", fst.male.matrix)

### Write Fst table
write.table(fst_male$pairwise.fst,"Fst_male.txt",sep="\t", row.names = FALSE, quote=FALSE)

### Melt the FST matrix to produce the nice graph
melted_fst_male <- reshape::melt(fst.male.matrix,na.rm = FALSE)

### Verify that the Fst values is indeed a numerical value
melted_fst_male$value <- as.numeric(as.character(melted_fst_male$value))

### Remove NA
melted_fst_male <- na.omit(melted_fst_male)
melted_fst_male$value <- as.numeric(as.character(melted_fst_male$value))

### Make the graph
male <- ggplot(data = melted_fst_male, aes(X2, X1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "yellow", mid ="white", high = "#009E73", space = "Lab", 
                       name="FST") +
  ylab("Sampling location A")+
  xlab("Sampling location B")+
  theme_bw()+ 
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 45, vjust = 1, 
                                                                size = 12, hjust = 1))+  
  coord_fixed()
male

### Save the graph
ggsave("fst_male.pdf", width=10, height=10)

###### COMPARISON OF FST ####

### Combine both Fst values
male_female <- cbind(melted_fst_male,melted_fst)
colnames(male_female) <- c("ReserveA","ReserveB","Male","ReserveC","ReserveD","Female")

### Gather Fst info
male_female_fst <- select(male_female, ReserveA,ReserveB, Male,Female)
male_female_fst_long <- male_female_fst %>% gather(SEX, FST, Male:Female)

# Plot
boxplot <- male_female_fst_long %>%
  ggplot( aes(x=SEX, y=FST, fill=SEX)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_manual(values=col2,name=" ") +
  theme_bw() +
  ylab("Index of genetic differentiation (FST)")+
  xlab("Gender information")
boxplot

### Test the significance
wilcox.test(FST ~SEX, male_female_fst_long) 

### Save the plot
ggsave("Fst_differences_male.pdf", width=5,height=5)

### Combine plots
cowplot::plot_grid(female, male,nrow=2,labels="auto")

###data:  FST by SEX
##W = 872, p-value = 0.004925
##alternative hypothesis: true location shift is not equal to 0
save.image("FST-palinurus-all.Rdata")
