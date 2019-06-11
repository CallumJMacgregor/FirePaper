########################################
### Plot trap start and finish times ###
########################################

### Clear the workspace
rm(list=ls())

# set the working directory as the script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### install if necessary and then load the libraries you need

j <- c("ggplot2","plyr","reshape","gridExtra","scales")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/BetaDiversityTurnover.txt", header=TRUE)

summary(dframe1)
                

betadiv <- dframe1[,c(1:2,12:15)]
colnames(betadiv) <- c("ExactSeason","Order",
                       "Change in flowers and pollinators","Change in flowers",
                       "Change in pollinators","Interaction rewiring")

# need to make this data long-form
betadiv_melt <- melt(betadiv, id=1:2)

# now re-order some variables to make them appear in the best order
betadiv_melt$ExactSeason <- ordered(betadiv_melt$ExactSeason, levels=c("Spring_1","Summer_1",
                                                                       "Autumn_2","Winter_2","Spring_2","Summer_2",
                                                                       "Autumn_3","Winter_3","Spring_3"))

summary(betadiv_melt$variable)

betadiv_melt$variable <- factor(betadiv_melt$variable, levels = c("Interaction rewiring","Change in pollinators",
                                                                  "Change in flowers","Change in flowers and pollinators"))





g1 <- ggplot(betadiv_melt, aes(x=ExactSeason, y=value, fill=variable))+
              geom_col(colour="black")+
              labs(x="Sampling period", y="Jaccard beta-diversity", fill="Interaction change due to...")+
  scale_fill_grey(start = 1, end = 0, guide=guide_legend(nrow = 2))+
              theme(panel.background=element_rect(fill="white"),
                    strip.background = element_blank(),
                    panel.grid.major.y=element_line(colour="gray70"),
                    panel.grid.minor.y=element_line(colour="gray70"),
                    panel.grid.major.x=element_blank(),
                    panel.border=element_rect(color="black",fill=F,size=1),
                    text=element_text(size=15),
                    axis.text=element_text(color="black"),
                    legend.position = "bottom")
  

g1

# w/o legend

g1a <- ggplot(betadiv_melt, aes(x=ExactSeason, y=value, fill=variable))+
  geom_col(colour="black")+
  labs(x="Sampling period", y="Jaccard beta-diversity", fill="Interaction change due to...")+
  scale_fill_grey(start = 1, end = 0)+
  theme(panel.background=element_rect(fill="white"),
        strip.background = element_blank(),
        legend.position = "none",
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor.y=element_line(colour="gray70"),
        panel.grid.major.x=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=15),
        axis.text=element_text(color="black"))


g1a


# now do one showing the raw number of unique and shared interactions

rawinters <- dframe1[,c(1:2,5:7)]
colnames(rawinters) <- c("ExactSeason","Order","Both networks","Burned network only","Unburned network only")

# need to make this data long-form
rawinters_melt <- melt(rawinters, id=1:2)

# now re-order some variables to make them appear in the best order
rawinters_melt$ExactSeason <- ordered(rawinters_melt$ExactSeason, levels=c("Spring_1","Summer_1",
                                                                       "Autumn_2","Winter_2","Spring_2","Summer_2",
                                                                       "Autumn_3","Winter_3","Spring_3"))

summary(rawinters_melt$variable)

rawinters_melt$variable <- factor(rawinters_melt$variable, levels = c("Burned network only","Both networks","Unburned network only"))

summary(rawinters_melt)



g2 <- ggplot(rawinters_melt, aes(x=ExactSeason, y=value, fill=variable))+
  geom_col(colour="black")+
  labs(x="Sampling period", y="Number of interactions", fill="Interaction observed in...")+
  scale_fill_grey(start = 1, end = 0, guide=guide_legend(nrow=2))+
  theme(panel.background=element_rect(fill="white"),
        strip.background = element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor.y=element_line(colour="gray70"),
        panel.grid.major.x=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=15),
        axis.text=element_text(color="black"),
        legend.position = "bottom")


g2

# w/o legend
g2a <- ggplot(rawinters_melt, aes(x=ExactSeason, y=value, fill=variable))+
  geom_col(colour="black")+
  labs(x=" ", y="Number of interactions", fill="Interaction observed in...")+
  scale_fill_grey(start = 1, end = 0)+
  theme(panel.background=element_rect(fill="white"),
        legend.position = "none",
        strip.background = element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor.y=element_line(colour="gray70"),
        panel.grid.major.x=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=15),
        axis.text=element_text(color="black"))


g2a


# we want to save the versions with legends independently and the versions without as a grid.arrange (it's a bit of a fiddle!)

m1 <- grid.arrange(g2,g1, ncol=1)
m1a <- grid.arrange(g2a,g1a, ncol=1)



ggsave("Fig7a.svg", plot = g1, device = "svg", path = "Results/", width = 29.7, height = 21.6, units = "cm")
ggsave("Fig7b.svg", plot = g2, device = "svg", path = "Results/", width = 29.7, height = 21.6, units = "cm")

ggsave("Fig7.svg", plot = m1, device = "svg", path = "Results/UpdatedFigs", width = 22.8, height = 33.2, units = "cm")


# and redo this in a different orientation

m2 <- grid.arrange(g2,g1, ncol = 2)

ggsave("Fig7alt.svg", plot = m2, device = "svg", path = "Results/UpdatedFigs", width = 50, height = 20, units = "cm")
