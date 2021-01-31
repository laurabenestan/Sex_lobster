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
### Download libraries
library(ggpplot2)
library(mapdata)
library(dplyr)
library(sf)
library(viridis)

### Import geographic coordinates file
sites <-read.table("../00-data/env_geo_palinurus_243ind.txt",header=TRUE, dec=".",sep="\t",na.strings="NA",strip.white=T)
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
  scale_color_viridis(discrete=TRUE)+
  facet_wrap(~GENDER)+
  theme(axis.text.x=element_text(colour="black",size=14))+
  theme(axis.text.y=element_text(colour="black", size=14))+
  theme(axis.title=element_text(colour="black",size=14))
map_lobster

### Save the map
ggsave("Map_lobster.pdf",width=10, height=10)
