' --------------------------------------------------------------------------   @Header
#'
#' @title Sampling for Palinurus elephas according to gender
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
library(ggplot2)
library(mapdata)
library(dplyr)
library(sf)
library(viridis)
library(colorBlindness)

### Pick up color-blind
safeColors
col2 <- c("#D55E00", "#009E73")

### Import geographic coordinates file
sites <-read.table("../00-data/env_geo_palinurus_241ind.txt",header=TRUE, dec=".",sep="\t",na.strings="NA",strip.white=T)
summary(sites)

### Download sex info
sex <- read.table("../00-data/256ind_palinurus_sex.txt")
colnames(sex) <- c("IND","GENDER")

### Merge sites and sex
sites_sex <- merge(sites, sex, by="IND")
length(which(is.na(sites_sex$GENDER)))

### Remove Na and Ind
target <- c("Female", "Male")
sites_sex2 <- sites_sex %>% filter(GENDER %in% target)

### Count the number of samples per lattitude and longitude points
sites_number <- sites_sex2 %>% group_by(LAT,LON,GENDER) %>%
  tally()

### Add reserve info
reserve <- read.table("../00-data/reserve_age.txt",header=TRUE,sep="\t",dec=".")

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
  geom_point(aes(x = LON, y = LAT,size=n, fill=GENDER), data=sites_number,shape = 21,alpha = 0.5)+
  geom_point(aes(x = LONGITUDE, y = LATITUDE), data=reserve,shape = 19)+
geom_text(aes(x=LONGITUDE, y=LATITUDE,label=NAME_EN),hjust=0, vjust=0, data=reserve)+
  xlab("Longitude") + ylab("Latitude") +
  theme_bw()+ 
  scale_fill_manual(values=col2)+
  facet_wrap(~GENDER)+
  theme(axis.text.x=element_text(colour="black",size=14))+
  theme(axis.text.y=element_text(colour="black", size=14))+
  theme(axis.title=element_text(colour="black",size=14))
map_lobster

### Save the map
ggsave("Figure1a_241.pdf",width=10, height=10)
