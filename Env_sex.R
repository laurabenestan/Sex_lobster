setwd("~/Documents/CEFE Langouste")

# Load Sex data
habitat <- read.table("00-Datasets/04-Environmental_data/sex_morpho_habitats_env_palinurus.txt", header=TRUE, sep="\t",dec=".")

# Keep only the columns sex and environmental variables
EnvMatrix_sex <- habitat[,c(1,2,5,7:30)]

# Calculate Moyenne et écart type entre males et femelles pour chaque variable
mean_F <- colMeans(EnvMatrix_sex[which(EnvMatrix_sex$Sex == "Female"), 3:27])
sd_F <- apply(EnvMatrix_sex[which(EnvMatrix_sex$Sex == "Female"), 3:27], 2, sd)
mean_M <- colMeans(EnvMatrix_sex[which(EnvMatrix_sex$Sex == "Male"), 3:27])
sd_M <- apply(EnvMatrix_sex[which(EnvMatrix_sex$Sex == "Male"), 3:27], 2, sd)

##### Perform t-tests to compare means
pval <- vector()
e=0
for (i in 3:27) {
  e=e+1
  var_F <- EnvMatrix_sex[which(EnvMatrix_sex$Sex == "Female"), i]
  var_M <- EnvMatrix_sex[which(EnvMatrix_sex$Sex == "Male"), i]
  res <- t.test(var_F, var_M, var.equal = TRUE)
  pval[e] <- res$p.value
}

## Plot the results (1 plot per variable)
library(ggplot2)

varnames <- colnames(EnvMatrix_sex)[3:27]

for (i in 1:24) {
  x = c("Female", "Male")
  y_i = c(mean_F[i], mean_M[i])
  sd_i = c(sd_F[i], sd_M[i])
  plot_i <- qplot(x,y_i)+geom_errorbar(aes(x=x, ymin=y_i-sd_i, ymax=y_i+sd_i), width=0.25)
  print(plot_i + labs(title= varnames[i],
                        y=varnames[i], x = ""))
  ggsave(paste(varnames[i],"plotSex.pdf", sep=""), height=7, width=7)
}
