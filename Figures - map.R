#########################################
####   Script to plot map of plots   ####
#########################################

### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("ggplot2","ggmap","plyr","rstudioapi")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### load up Callum's custom set of functions
f <- c("MultiplotFunction.R")
lapply(f, source)




### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.csv("Data/Coordinates.csv", header=TRUE)

dframe1$Treatment <- factor(ifelse(grepl("NF",dframe1$Plot),"NoFire","Fire"))
summary(dframe1)

means <- ddply(dframe1, .(Plot,Treatment), numcolwise(mean))


# check the polygons plot ok with a ggplot

gtest <- ggplot(dframe1, aes(x=LongDec,y=LatDec,group=Plot))+
  geom_polygon(aes(fill=Treatment))

gtest



# now we want to apply them to a map

map1a <- get_openstreetmap(bbox = c(left=-7.89, bottom = 37.158, right = -7.83, top = 37.202), scale=31000)

map1 <- get_map(location = c(lon=-7.86,lat=37.18), maptype="terrain", zoom=13,color="bw")


ggmap(map1)


# set things up for a scale bar
# get match attributes
bb <- attr(map1,"bb")


# figure out points to define scale bar
sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                   lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                   lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                   lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))

sbar$distance <- geosphere::distVincentyEllipsoid(c(sbar$lon.start,sbar$lat.start),
                                                  c(sbar$lon.end,sbar$lat.end))

ptspermm <- 2.83464567  # apparently, we "need this because geom_text uses mm, and themes use pts. Urgh."

# set the length of the scale bar - here 20km
scalebar.length <- 2
sbar$lon.end <- sbar$lon.start +
  ((sbar$lon.end-sbar$lon.start)/sbar$distance)*scalebar.length*1000


g1 <- ggmap(map1, extent = "device", crop = T) +
  geom_point(data = means, size=5,colour="black",fill="white",stroke=2,alpha=0.8,
             aes(x=LongDec,y=LatDec,shape=Treatment))+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  theme(legend.justification=c(0,0), legend.position=c(0.07,0.77))+ 
  theme(legend.text = element_text(size=20))+
  theme(legend.title = element_text(size=20))+
  theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"))+
  geom_segment(data = sbar, size=1.5,
               aes(x = lon.start+0.065,
                   xend = lon.end+0.065,
                   y = lat.start+0.06,
                   yend = lat.end+0.06),
               arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                           ends = "both", type = "open")) +
  geom_text(data = sbar,
            aes(x = ((lon.start + lon.end)/2)+0.065,
                y = lat.start + 0.025*(bb$ur.lat - bb$ll.lat) + 0.06,
                label = paste(format(scalebar.length),
                              'km')),
            hjust = 0.5,
            vjust = 0,
            size = 20/ptspermm)  +
  coord_map(projection = "mercator",
            xlim=c(bb$ll.lon, bb$ur.lon),
            ylim=c(bb$ll.lat, bb$ur.lat))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.key = element_blank()) 

               
g1
               
        


# we also need a big map to show where this map is


map2 <- get_map(location = c(lon=-7.86,lat=38.08), maptype="toner-lite",zoom=8,color="color")

ggmap(map2)



# set things up for a scale bar
# get match attributes
bb2 <- attr(map2,"bb")


# figure out points to define scale bar
sbar2 <- data.frame(lon.start = c(bb2$ll.lon + 0.1*(bb2$ur.lon - bb2$ll.lon)),
                   lon.end = c(bb2$ll.lon + 0.25*(bb2$ur.lon - bb2$ll.lon)),
                   lat.start = c(bb2$ll.lat + 0.1*(bb2$ur.lat - bb2$ll.lat)),
                   lat.end = c(bb2$ll.lat + 0.1*(bb2$ur.lat - bb2$ll.lat)))

sbar2$distance <- geosphere::distVincentyEllipsoid(c(sbar2$lon.start,sbar2$lat.start),
                                                  c(sbar2$lon.end,sbar2$lat.end))

ptspermm <- 2.83464567  # apparently, we "need this because geom_text uses mm, and themes use pts. Urgh."

# set the length of the scale bar - here 20km
scalebar.length2 <- 20
sbar2$lon.end <- sbar2$lon.start +
  ((sbar2$lon.end-sbar2$lon.start)/sbar2$distance)*scalebar.length2*1000


# construct the map


g2 <- ggmap(map2, extent = "device", crop = T) +
  geom_rect(aes(xmin=-7.92,xmax=-7.80,ymin=37.138,ymax=37.222),
            color="black",alpha=0,size=2)+
  geom_segment(data = sbar2, size=1.5,
               aes(x = lon.start,
                   xend = lon.end,
                   y = lat.start,
                   yend = lat.end),
               arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                           ends = "both", type = "open")) +
  geom_text(data = sbar2,
            aes(x = ((lon.start + lon.end)/2),
                y = lat.start + 0.025*(bb2$ur.lat - bb2$ll.lat),
                label = paste(format(scalebar.length2),
                              'km')),
            hjust = 0.5,
            vjust = 0,
            size = 20/ptspermm)  +
  coord_map(projection = "mercator",
            xlim=c(bb2$ll.lon, bb2$ur.lon),
            ylim=c(bb2$ll.lat, bb2$ur.lat))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.key = element_blank())


g2


m1 <- grid.arrange(g2,g1,ncol=2)

m1


ggsave("Fig1maps.svg", plot = m1, device = "svg", path = "Results/UpdatedFigs", width = 44, height = 24, units = "cm")


### try a slightly different approach to get just country outlines for map 1

portugal <- map_data("world", region = "Portugal")
spain <- map_data("world", region = "Spain")
france <- map_data("world", region = "France")
morocco <- map_data("world", region = "Morocco")




g3 <- ggplot() + 
  geom_polygon(data = portugal, aes(x=long, y = lat, group = group), fill = "white", color = "black") + 
  coord_fixed(1.3) +
  geom_polygon(data = spain, aes(x=long, y = lat, group = group), fill = "white", color = "black") + 
  coord_fixed(1.3)+
  geom_polygon(data = france, aes(x=long, y = lat, group = group), fill = "white", color = "black") + 
  coord_fixed(1.3)+
  geom_polygon(data = morocco, aes(x=long, y = lat, group = group), fill = "white", color = "black") + 
  coord_fixed(1.3)+
  geom_rect(aes(xmin=-7.92,xmax=-7.80,ymin=37.138,ymax=37.222),
            color="black",alpha=0,size=2)+
  geom_segment(data = sbar2, size=1.5,
               aes(x = lon.start,
                   xend = lon.end,
                   y = lat.start,
                   yend = lat.end),
               arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                           ends = "both", type = "open")) +
  geom_text(data = sbar2,
            aes(x = ((lon.start + lon.end)/2),
                y = lat.start + 0.025*(bb2$ur.lat - bb2$ll.lat),
                label = paste(format(scalebar.length2),
                              'km')),
            hjust = 0.5,
            vjust = 0,
            size = 20/ptspermm)+
  coord_map(projection = "mercator",
            xlim=c(bb2$ll.lon, bb2$ur.lon),
            ylim=c(bb2$ll.lat, bb2$ur.lat))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


g3




### try a slightly different approach to get a terrain map for map 2


map1b <- get_map(location = c(lon=-7.86,lat=37.18), maptype="terrain", zoom=13,color="bw")

ggmap(map1b)

g4 <- ggmap(map1b, extent = "device", crop = T) +
  geom_point(data = means, size=5,colour="black",fill="white",stroke=2,alpha=1,
             aes(x=LongDec,y=LatDec,shape=Treatment))+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  theme(legend.justification=c(0,0), legend.position=c(0.07,0.77))+ 
  theme(legend.text = element_text(size=20))+
  theme(legend.title = element_text(size=20))+
  theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"))+
  geom_segment(data = sbar, size=1.5,
               aes(x = lon.start+0.065,
                   xend = lon.end+0.065,
                   y = lat.start+0.06,
                   yend = lat.end+0.06),
               arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                           ends = "both", type = "open")) +
  geom_text(data = sbar,
            aes(x = ((lon.start + lon.end)/2)+0.065,
                y = lat.start + 0.025*(bb$ur.lat - bb$ll.lat) + 0.06,
                label = paste(format(scalebar.length),
                              'km')),
            hjust = 0.5,
            vjust = 0,
            size = 20/ptspermm)  +
  coord_map(projection = "mercator",
            xlim=c(bb$ll.lon, bb$ur.lon),
            ylim=c(bb$ll.lat, bb$ur.lat))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.key = element_blank()) 


g4




grid.arrange(g3,g4, ncol = 2)
