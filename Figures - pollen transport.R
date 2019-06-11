##########################################################
####   Script for plotting pollen transport figures   ####
##########################################################

# we need to keep the necessary bits (wrangling and final models) from scripts 5 and 6 first

### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","arm","MASS","scales","AICcmodavg","svglite","effects","plyr","gridExtra","reshape2")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

# for some reason the package glmmADMB won't install via the usual methods, so:
#install.packages("R2admb")
#install.packages("glmmADMB", 
#                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                         getOption("repos")),
#                 type="source")
library(glmmADMB)



### load up Callum's custom set of functions
k <- c("CheckResidsFunction.R","CheckConvergenceFunction.R","MultiplotFunction.R")
lapply(k,source)

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

### tell R that SampleID and SlideNumber should be treated as factors
dframe1$SampleID <- factor(dframe1$SampleID)
dframe1$SlideNumber <- factor(dframe1$SlideNumber)
dframe1$Month<-ordered(dframe1$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
dframe1$Year <- factor(dframe1$Year)
dframe1$Treatment <- relevel(dframe1$Treatment, "NoFire")

summary(dframe1)
names(dframe1)


# we need a variable to allow the networks to be paired:

dframe1$ExactSeason <- do.call(paste, c(dframe1[c("Season","Year")], sep = "_"))
dframe1$ExactSeason <- as.factor(dframe1$ExactSeason)
summary(dframe1$ExactSeason)

# we also want this variable to be ordered so that we can look for any trend over time
# I've just crudely put this into a text file that we can import and merge in

SeasonOrder <- read.csv("Data/SeasonOrder.csv",header=T)

SeasonOrder$ExactSeason <- factor(SeasonOrder$ExactSeason, levels=SeasonOrder$ExactSeason[order(SeasonOrder$Order)], ordered=T)
summary(SeasonOrder)

dframe1 <- merge(SeasonOrder,dframe1)

summary(dframe1$Order)


### total up the pollen grains for each insect
dframe1$PollenLoad <- rowSums(dframe1[c(15:80)])

### total up the number of pollen species for each insect (i.e. how many columns do not contain 0)
dframe1$PollenTypes <- rowSums(dframe1[c(15:80)] != 0)

### create a binary (yes/no) variable for whether each insect is carrying any pollen
dframe1$PollenYN <- ifelse(dframe1$PollenTypes==0,0,1)
summary(dframe1$PollenYN)

dframe1B <- dframe1[dframe1$Treatment=="Fire",]
dframe1U <- dframe1[dframe1$Treatment=="NoFire",]

summary(dframe1B$PollenYN)
summary(dframe1U$PollenYN)


### create a subset dframe containing only the interactions
interactions <- subset(dframe1, select=-c(SampleID,Date,Site,SlideNumber,PollenCount,Treatment,SamplingDay,Sample,PollenTypes,PollenLoad,PollenYN,Season,Month,Year,ExactSeason,Order))
summary(interactions)

### create a subset dataframe containing only pollen-carrying moths
dframe1P <- dframe1[dframe1$PollenYN==1,]


dframe1Sum <- dframe1[dframe1$Season=="Summer",]
dframe1Spr <- dframe1[dframe1$Season=="Spring",]
dframe1Win <- dframe1[dframe1$Season=="Winter",]
dframe1Aut <- dframe1[dframe1$Season=="Autumn",]



### Let's first look at pollen load per-moth


model1PGa <- lmer(log(PollenLoad) ~ Treatment*Season + Order
                  + (1|Site),
                  data = dframe1P)

summary(model1PGa)
drop1(model1PGa, test = "Chi")



### Let's now look at pollen types per-moth

# Gaussian with log transformation

model2P <- glmer(PollenTypes ~ Treatment * Season + Treatment * Order # fixed effects
                 + (1|Site), # random effects
                 family = poisson (link = "log"),
                 data = dframe1P)


summary(model2P)
drop1(model2P, test = "Chi")


### finally, let's look at proportion of moths carrying pollen

### this data is definitely binomial, so we don't need to worry too much about model selection or residuals:

model3Ba <- glmer(PollenYN~Treatment*Season + Order
                  + (1|Site),
                  family = binomial (link = "logit"),
                  data = dframe1)


summary(model3Ba)
drop1(model3Ba, test="Chi")


### now prep the per-sample analyses

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe2<-read.table("Data/SamplesNoct.txt", header=TRUE)

summary(dframe2) # Check it's imported correctly

dframe2$Month<-ordered(dframe2$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
dframe2$Treatment <- relevel(dframe2$Treatment, "NoFire")

# we need a variable to allow the networks to be paired:

dframe2$ExactSeason <- do.call(paste, c(dframe2[c("Season","Year")], sep = "_"))
dframe2$ExactSeason <- as.factor(dframe2$ExactSeason)
summary(dframe2$ExactSeason)

# we also want this variable to be ordered so that we can look for any trend over time
# I've just crudely put this into a text file that we can import and merge in

dframe2 <- merge(SeasonOrder,dframe2)

summary(dframe2$Order)



### total up the pollen grains for each sample
dframe2$PollenCount <- rowSums(dframe2[c(14:79)])

### total up the number of pollen species for each sample (i.e. how many columns do not contain 0)
dframe2$PollenTypes <- rowSums(dframe2[c(14:79)] != 0)


### now you're ready to start looking for patterns!



### Let's first look at pollen load per-moth

model4Ga <- lmer(log(PollenCount+1) ~ Treatment + Season + Order
                 + (1|Site),
                 data = dframe2)

summary(model4Ga)
drop1(model4Ga, test = "Chi")



### Let's now look at pollen types per-moth

# Poisson

model5Pa <- glmer(PollenTypes ~ Treatment*Season  + Order # fixed effects
                  + (1|Site), # random effects
                  family = poisson (link = "log"),
                  data = dframe2)

summary(model5Pa)
drop1(model5Pa, test = "Chi")



### we have the necessary models so let's start plotting figures 

# we want a figure each from individual moth pollen load (carriers only) and types,
# moth probability of carrying
# sample level pollen load and types

# moth pollen load
summary(model1PGa)
drop1(model1PGa, test="Chi")

# create the framework for the plot
newdata1 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),PollenLoad=1)
errors1 <- data.frame(effect(c("Treatment","Season"),model1PGa))

newdata1 <- merge(newdata1,errors1)

newdata1$tfit <- exp(newdata1$fit)
newdata1$tlower <- exp(newdata1$lower)
newdata1$tupper <- exp(newdata1$upper)


newdata1$Treatment <- revalue(newdata1$Treatment, c("Fire"="Burned","NoFire"="Unburned"))



# construct the plot
g1 <- ggplot(newdata1,
             aes(x=Treatment, y=tfit, fill=Season))+
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  xlab(" ")+ ylab("Pollen load")+ 
  ylim(0,20)+
  facet_grid(.~Season, switch="x", space = "free_x", scales = "free_x")+
  guides(fill=F)+
  ggtitle("Individual moths\n")+
  theme(legend.justification=c(0,0), legend.position="none",
        legend.text = element_text(size=30),
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        panel.background=element_rect(fill="white"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.y=element_text(color="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(size = 5, fill="white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g1




## repeat this for the other three


# moth pollen types
summary(model2P)
drop1(model2P, test="Chi")

# create the framework for the plot
newdata2 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),PollenTypes=1)
errors2 <- data.frame(effect(c("Treatment","Season"),model2P))

newdata2 <- merge(newdata2,errors2)


newdata2$Treatment <- revalue(newdata2$Treatment, c("Fire"="Burned","NoFire"="Unburned"))



# construct the plot
g2 <- ggplot(newdata2,
             aes(x=Treatment, y=fit, fill=Season))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  ylim(0,4)+
  xlab("Season")+ ylab("Pollen species richness")+ 
  facet_grid(.~Season, switch="x", space = "free_x", scales = "free_x")+
  guides(fill=F)+
  theme(legend.justification=c(0,0), legend.position="none",
        legend.text = element_text(size=30),
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        panel.background=element_rect(fill="white"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.y=element_text(color="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(size = 5, fill="white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g2



# plant abundance
summary(model3Ba)
drop1(model3Ba, test = "Chi")

# create the framework for the plot
newdata3 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),PollenYN=1)
effects <- data.frame(effect(c("Treatment","Season"),model3Ba))

newdata3 <- merge(newdata3,effects)

newdata3$Treatment <- revalue(newdata3$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g3 <- ggplot(newdata3,
             aes(x=Treatment, y=fit, fill=Season))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  scale_y_continuous(limits = c(0,1), oob=squish)+
  xlab("Season")+ ylab("Proportion of moths transporting pollen")+
  facet_grid(.~Season, switch="x", space = "free_x", scales = "free_x")+
  guides(fill=F)+
  theme(legend.text = element_text(size=30),
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        panel.background=element_rect(fill="white"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.y=element_text(color="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill="white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g3



# sample pollen load
summary(model4Ga)
drop1(model4Ga, test="Chi")

# create the framework for the plot
newdata4 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),PollenCount=1)
errors4 <- data.frame(effect(c("Treatment","Season"), model4Ga))

newdata4 <- merge(newdata4,errors4)

newdata4$tfit <- exp(newdata4$fit)
newdata4$tlower <- exp(newdata4$lower)
newdata4$tupper <- exp(newdata4$upper)



newdata4$Treatment <- revalue(newdata4$Treatment, c("Fire"="Burned","NoFire"="Unburned"))



# construct the plot
g4 <- ggplot(newdata4,
             aes(x=Treatment, y=tfit, fill=Season))+
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  ylim(0,8000)+
  xlab(" ")+ ylab(" ")+ 
  facet_grid(.~Season, switch="x", space = "free_x", scales = "free_x")+
  guides(fill=F)+
  ggtitle("Accumulated samples\n")+
  theme(legend.justification=c(0,0), legend.position="none",
        legend.text = element_text(size=30),
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        panel.background=element_rect(fill="white"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.y=element_text(color="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(size = 5, fill="white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g4


# sample pollen load
summary(model5Pa)
drop1(model5Pa, test="Chi")

# create the framework for the plot
newdata5 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),PollenTypes=1)
errors5 <- data.frame(effect(c("Treatment","Season"), model5Pa))

newdata5 <- merge(newdata5,errors5)



newdata5$Treatment <- revalue(newdata5$Treatment, c("Fire"="Burned","NoFire"="Unburned"))



# construct the plot
g5 <- ggplot(newdata5,
             aes(x=Treatment, y=fit, fill=Season))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  ylim(0,25)+
  xlab("Season")+ ylab(" ")+ 
  facet_grid(.~Season, switch="x", space = "free_x", scales = "free_x")+
  guides(fill=F)+
  theme(legend.text = element_text(size=30),
        legend.background = element_rect(linetype="solid",colour ="black"),
        panel.background=element_rect(fill="white"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.y=element_text(color="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill="white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g5


# now stitch them into a multiplot

m1 <- grid.arrange(g1,g4,g2,g5, ncol=2, nrow=2)


# before we export this we want to move the legend to be outside the plots

# function to extract a legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(g5)

m1a <- grid.arrange(arrangeGrob(g1,g4,
                                g2,g5 + theme(legend.position="none"),
                                nrow=2, ncol=2),
                     mylegend, ncol =2, widths=c(10,2))




ggsave("Fig5.svg", plot = m1a, device = "svg", path = "Results/UpdatedFigs", width = 48, height = 40, units = "cm")



# and remind ourselves of the other one

g3

ggsave("Fig4.svg", plot = g3, device = "svg", path = "Results/UpdatedFigs", width = 22, height = 16, units = "cm")



## now try merging these together

m1b <- grid.arrange(g1,g4,
                    g2 + xlab(" "),g5 + theme(legend.position = "none"),
                    g3 + theme(legend.position = "none"),mylegend,
                    nrow = 3, ncol = 2)



ggsave("Fig4alt.svg", plot = m1b, device = "svg", path = "Results/UpdatedFigs", width = 48, height = 60, units = "cm")


### we also want to try and get figures showing effect over time rather than seasonal effect

# moth pollen load
summary(model1PGa)
drop1(model1PGa, test="Chi")

# create the framework for the plot
newdata6 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),PollenLoad=1)
errors6 <- data.frame(effect(c("Treatment","Order"),model1PGa,xlevels=list(Order=1:9)))

newdata6 <- merge(newdata6,errors6)

newdata6$tfit <- 10^newdata6$fit
newdata6$tlower <- 10^newdata6$lower
newdata6$tupper <- 10^newdata6$upper


newdata6$Treatment <- revalue(newdata6$Treatment, c("Fire"="Burned","NoFire"="Unburned"))



# construct the plot
g6 <- ggplot(newdata6,
             aes(x=Order, y=tfit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=tlower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=tupper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab(" ")+ ylab("Pollen load")+ 
  ylim(0,10)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
  ggtitle("Individual moths\n")+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text=element_text(size=17,color="black"),
        plot.title=element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.position="none")

# visualize the plot
g6




## repeat this for the other three


# moth pollen types
summary(model2P)
drop1(model2P, test="Chi")

# create the framework for the plot
newdata7 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),PollenLoad=1)
errors7 <- data.frame(effect(c("Treatment","Order"),model2P,xlevels=list(Order=1:9)))

newdata7 <- merge(newdata7,errors7)


newdata7$Treatment <- revalue(newdata7$Treatment, c("Fire"="Burned","NoFire"="Unburned"))



# construct the plot
g7 <- ggplot(newdata7,
             aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+  xlab(" ")+ ylab("Pollen types")+ 
  ylim(0,5)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
  xlab("Seasons since\nstudy began")+ ylab("Pollen species richness")+ 
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text=element_text(size=17,color="black"),
        plot.title=element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.position="none")



# visualize the plot
g7


# plant abundance
summary(model3Ba)
drop1(model3Ba, test = "Chi")

# create the framework for the plot
newdata8 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),PollenLoad=1)
errors8 <- data.frame(effect(c("Treatment","Order"),model3Ba,xlevels=list(Order=1:9)))


newdata8 <- merge(newdata8,errors8)

newdata8$Treatment <- revalue(newdata8$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g8 <- ggplot(newdata8,
             aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+  xlab(" ")+ ylab("Pollen types")+ 
  scale_y_continuous(limits = c(0,1), oob=squish)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
  xlab("Seasons since\nstudy began")+ ylab("Proportion of moths transporting pollen")+ 
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text=element_text(size=17,color="black"),
        plot.title=element_text(hjust=0.5),
        legend.key = element_blank())+
  scale_fill_discrete(name="Treatment")



# visualize the plot
g8





# sample pollen load
summary(model4Ga)
drop1(model4Ga, test="Chi")

# create the framework for the plot
newdata9 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),PollenLoad=1)
errors9 <- data.frame(effect(c("Treatment","Order"),model4Ga,xlevels=list(Order=1:9)))


newdata9 <- merge(newdata9,errors9)

newdata9$Treatment <- revalue(newdata9$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


newdata9$tfit <- 10^newdata9$fit
newdata9$tlower <- 10^newdata9$lower
newdata9$tupper <- 10^newdata9$upper



newdata9$Treatment <- revalue(newdata9$Treatment, c("Fire"="Burned","NoFire"="Unburned"))



# construct the plot
g9 <- ggplot(newdata9,
             aes(x=Order, y=tfit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=tlower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=tupper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+  xlab(" ")+ ylab("Pollen types")+ 
  scale_y_continuous(limits = c(0,2000), oob=squish)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
  xlab(" ")+ ylab(" ")+
  ggtitle("Accumulated samples\n")+ 
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text=element_text(size=17,color="black"),
        plot.title=element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.position="none")

# visualize the plot
g9



# sample pollen load
summary(model5Pa)
drop1(model5Pa, test="Chi")

# create the framework for the plot
newdata10 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),PollenLoad=1)
errors10 <- data.frame(effect(c("Treatment","Order"),model5Pa,xlevels=list(Order=1:9)))


newdata10 <- merge(newdata10,errors10)


newdata10$Treatment <- revalue(newdata10$Treatment, c("Fire"="Burned","NoFire"="Unburned"))



# construct the plot
g10 <- ggplot(newdata10,
             aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+  xlab(" ")+ ylab("Pollen types")+ 
  scale_y_continuous(limits = c(0,15), oob=squish)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=FALSE)+
  xlab("Seasons since\nstudy began")+ ylab(" ")+ 
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text=element_text(size=17,color="black"),
        plot.title=element_text(hjust=0.5),
        legend.position="none")


# visualize the plot
g10



# now stitch them into a multiplot
m2 <- grid.arrange(g6,g9,g7,g10, ncol=2, nrow=2)


mylegend <- g_legend(g8)

m2a <- grid.arrange(arrangeGrob(g6,g9,
                                g7,g10,
                                nrow=2, ncol=2),
                    mylegend, ncol =2, widths=c(10,2))


ggsave("FigS5.svg", plot=m2a, device = "svg", path = "Results/UpdatedFigs", width = 48, height = 40, units = "cm")


# and remind ourselves of the other one

g8

ggsave("FigS3.svg", plot = g8, device = "svg", path = "Results/UpdatedFigs", width = 33, height = 24, units = "cm")



#### finally, we want to try and plot the typical frequency density of number of pollen grains per plant species within each trap

# start with dframe1, keeping only SampleID, Sample and the pollen data
dframe1fd <- dframe1[,c(3,11,15:80)]

# melt into long form
dframe1fd.melt <- melt(dframe1fd, id = c("SampleID","Sample"))

# create a grouping variable for plant species in each sample
dframe1fd.melt$SampleSpecies <- factor(paste(dframe1fd.melt$Sample,dframe1fd.melt$variable,sep="-"))

max(dframe1fd.melt$value)

dframe1fd.melt <- dframe1fd.melt[!(dframe1fd.melt$value == 4058), ]
dframe1fd.melt <- dframe1fd.melt[!(dframe1fd.melt$value == 1437), ]

# now we want to create a geom_density overlaying all the different layers of SampleSpecies

g11 <- ggplot(dframe1fd.melt, aes(x = value))+
  geom_density()

g11

