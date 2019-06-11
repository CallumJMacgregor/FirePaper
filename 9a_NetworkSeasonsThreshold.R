###########################################################
####   Script for pollen transport analysis by season  ####
###########################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS","glmmADMB","scales","AICcmodavg","effects","svglite","plyr","gridExtra")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

# for some reason the package glmmADMB won't install via the usual methods, so if it's not already installed:
#install.packages("R2admb")
#install.packages("glmmADMB", 
#                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                         getOption("repos")),
#                 type="source")



### load up Callum's custom set of functions
k <- c("CheckResidsFunction.R","CheckConvergenceFunction.R","MultiplotFunction.R")
lapply(k,source)


### Flower visitation ###


### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/NetworkmetricsFVbySeasonThreshold.txt", header=TRUE)
dframe1$Treatment <- relevel(dframe1$Treatment, "NoFire")

summary(dframe1) # Check it's imported correctly

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
summary(dframe1$ExactSeason)


# finally, there is a row in the dataframe with NAs in it - this is Winter_2_Fire, where there was only a single interaction above the threshold
# the network metrics are therefore pretty meaningless so let's cut that row

dframe1 <- dframe1[complete.cases(dframe1), ]


### now you're ready to start looking for patterns!

### next - linkage density

hist(dframe1$linkage.density)   # this looks pretty good for a Poisson distribution, but has non-integers, so let's use Gaussian
plot(linkage.density ~ Treatment, data = dframe1)

# construct the model

model1LD <- lm(linkage.density ~ Treatment*Season + Treatment*Order,
               data=dframe1)

# inspect and test
summary(model1LD)
drop1(model1LD, test = "F")

chkres(model1LD,dframe1$Treatment,dframe1$Season) # these are fine
#
#
#
#
#
#

# try again without non-significant interaction term


model1LDa <- lm(linkage.density ~ Treatment + Season + Order,
                  data=dframe1)

# inspect and test
summary(model1LDa)
drop1(model1LDa, test = "F")

chkres(model1LDa,dframe1$Treatment,dframe1$Season)  # nothing concerning in these
#
#
#
#
#
#

### next - vulnerability

hist(dframe1$vulnerability.LL)  # not perfectly normal but not drastically non-normal
plot(vulnerability.LL ~ Treatment, data = dframe1)

# construct the model

model1V <- lm(vulnerability.LL ~ Treatment*Season + Treatment*Order,
                data = dframe1)

summary(model1V)
drop1(model1V, test = "F")

chkres(model1V, dframe1$Treatment, dframe1$Season)
#
#
#
#

# try again without non-significant interaction term


model1Va <- lm(vulnerability.LL ~ Treatment+Season+Order,
                  data=dframe1)

# inspect and test
summary(model1Va)
drop1(model1Va, test = "F")

chkres(model1Va,dframe1$Treatment,dframe1$Season)  # nothing concerning in these
#
#
#
#


### next - generality

hist(dframe1$generality.HL)
plot(generality.HL ~ Treatment, data = dframe1)

# construct the model

model1G <- lm(generality.HL ~ Treatment*Season + Treatment*Order,
                data=dframe1)

summary(model1G)
drop1(model1G, test = "F")

chkres(model1G,dframe1$Treatment,dframe1$Season)
#
#
#
#

# try again without non-significant interaction term

model1Ga <- lm(generality.HL ~ Treatment+Season+Order,
                 data=dframe1)

summary(model1Ga)
drop1(model1Ga, test = "F")

chkres(model1Ga,dframe1$Treatment,dframe1$Season)  # these are mostly ok
#
#
#
#

### next - robustness

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$robustness.HL.mean)                   # this looks roughly normal
plot(robustness.HL.mean ~ Treatment, data = dframe1)


# construct models using Date and Site as random effects

model1RB <- lm(robustness.HL.mean ~ Treatment*Season + Treatment*Order,
                 data = dframe1)

# inspect and test the model
summary(model1RB)
drop1(model1RB, test="F")  


# check the model's residuals
chkres(model1RB, dframe1$Treatment, dframe1$Season)  # these residuals aren't too bad - a hint of a positive trend
#
#
#
#

# drop non-sig interactions
model1RBa <- lm(robustness.HL.mean ~ Treatment + Season + Order,
                data = dframe1)

summary(model1RBa)
drop1(model1RBa, test="F")

# check the model's residuals
chkres(model1RBa, dframe1$Treatment, dframe1$Season)  # these residuals aren't too bad - a hint of a positive trend
#
#
#
#



# finally, because there is an effect on linkage density, but none on generality or vulnerability, we want to check niche overlap

# Plot it against treatment
plot(dframe1$niche.overlap.HL ~ dframe1$Treatment)
plot(dframe1$niche.overlap.HL ~ dframe1$Season)

hist(dframe1$niche.overlap.HL) # looks almost normal


# try a Gaussian

model1NO <- lm(niche.overlap.HL ~ Treatment * Season + Treatment * Order,
               data = dframe1)

summary(model1NO)
drop1(model1NO, test="F")

chkres(model1NO, dframe1$Treatment,dframe1$Season) # mostly ok, not perfect, but...
#
#
#
#

# try removing non-sig interactions
model1NOa <- lm(niche.overlap.HL ~ Treatment + Season + Order,
                data = dframe1)

summary(model1NOa)
drop1(model1NOa, test="F")

chkres(model1NOa, dframe1$Treatment,dframe1$Season) # better
#
#
#
#





### Pollen transport ###



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe2<-read.table("Data/NetworkMetricsPTbySeasonThreshold.txt", header=TRUE)
dframe2$Treatment <- relevel(dframe2$Treatment, "NoFire")

summary(dframe2) # Check it's imported correctly

# we need a variable to allow the networks to be paired:

dframe2$ExactSeason <- do.call(paste, c(dframe2[c("Season","Year")], sep = "_"))
dframe2$ExactSeason <- as.factor(dframe2$ExactSeason)
summary(dframe2$ExactSeason)

dframe2 <- merge(SeasonOrder,dframe2)
summary(dframe2$ExactSeason)

# and again, cut out Winter_2_Fire
dframe2 <- dframe2[complete.cases(dframe2), ]


### now you're ready to start looking for patterns!


### linkage density

hist(dframe2$linkage.density)   # this looks pretty good for a Gaussian distribution
plot(linkage.density ~ Treatment, data = dframe2)
plot(linkage.density ~ Season, data = dframe2)
plot(linkage.density ~ Order, data = dframe2)

# construct the model

model2LD <- lm(linkage.density ~ Treatment*Season + Treatment*Order,
                 data=dframe2)

# inspect and test
summary(model2LD)
drop1(model2LD, test = "F")

chkres(model2LD,dframe2$Treatment,dframe2$Season)
#
#
#
#


# try again without non-significant interaction term


model2LDa <- lm(linkage.density ~ Treatment + Season + Order,
                  data=dframe2)

# inspect and test
summary(model2LDa)
drop1(model2LDa, test = "F")

chkres(model2LDa,dframe2$Treatment,dframe2$Season)  # nothing too worrying in these
#
#
#
#


### next - vulnerability

hist(dframe2$vulnerability.LL)  # not perfectly normal but not drastically non-normal
plot(vulnerability.LL ~ Treatment, data = dframe2)

# construct the model

model2V <- lm(vulnerability.LL ~ Treatment*Season + Treatment*Order,
                data = dframe2)

summary(model2V)
drop1(model2V, test = "F")

chkres(model2V, dframe2$Treatment, dframe2$Season)
#
#
#
#


# try again without non-significant interaction term


model2Va <- lm(vulnerability.LL ~ Treatment+Season+Order,
                 data=dframe2)

# inspect and test
summary(model2Va)
drop1(model2Va, test = "F")

chkres(model2Va,dframe2$Treatment,dframe2$Season)  # nothing concerning in these
#
#
#
#

### next - generality

hist(dframe2$generality.HL)
plot(generality.HL ~ Treatment, data = dframe2)
plot(generality.HL ~ Season, data = dframe2)
plot(generality.HL ~ Order, data = dframe2)


# construct the model

model2G <- lm(generality.HL ~ Treatment*Season + Treatment*Order,
                data=dframe2)

summary(model2G)
drop1(model2G, test = "F")

chkres(model2G,dframe2$Treatment,dframe2$Season)
#
#
#
#


# try again without non-significant interaction term

model2Ga <- lm(generality.HL ~ Treatment+Season+Order,
                 data=dframe2)

summary(model2Ga)
drop1(model2Ga, test = "F")

chkres(model2Ga,dframe2$Treatment,dframe2$Season)  # these are mostly ok
#
#
#
#


### next - robustness

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$robustness.HL.mean)                   # this looks like a normal distribution
plot(robustness.HL.mean ~ Treatment, data = dframe2)
plot(robustness.HL.mean ~ Season, data = dframe2)
plot(robustness.HL.mean ~ Order, data = dframe2)


# construct models using Date and Site as random effects

model2RB <- lm(robustness.HL.mean ~ Treatment*Season + Treatment*Order,
                 data = dframe2)

# inspect and test the model
summary(model2RB)
drop1(model2RB, test="F")  


# check the model's residuals
chkres(model2RB, dframe2$Treatment, dframe2$Season)  # these residuals aren't too bad - a hint of a positive trend
#
#
#
#


# drop non-sig interaction
model2RBa <- lm(robustness.HL.mean ~ Treatment + Season + Order,
                 data = dframe2)

summary(model2RBa)
drop1(model2RBa, test = "F")

chkres(model2RBa, dframe2$Treatment, dframe2$Season)
#
#
#
#

# finally, because there is an effect on linkage density, but none on generality or vulnerability, we want to check niche overlap

# Plot it against treatment
plot(dframe2$niche.overlap.HL ~ dframe2$Treatment)
plot(dframe2$niche.overlap.HL ~ dframe2$Season)

hist(dframe2$niche.overlap.HL) # looks almost normal


# try a Gaussian

model2NO <- lm(niche.overlap.HL ~ Treatment * Season + Treatment * Order,
               data = dframe2)

summary(model2NO)
drop1(model2NO, test="F")

chkres(model2NO, dframe2$Treatment,dframe2$Season)
#
#
#
#

# try removing non-sig interactions
model2NOa <- lm(niche.overlap.HL ~ Treatment + Season + Order,
                data = dframe2)

summary(model2NOa)
drop1(model2NOa, test="F")

chkres(model2NOa, dframe2$Treatment,dframe2$Season)
#
#
#
#



### figures

# we want a 2x2 matrix of figures, showing (for pollen transport) the effects
# on linkage density, robustness, vulnerability and generality

# let's just do prediction +- 95CI

### pollen transport

## linkage density

summary(model2LDa)
drop1(model2LDa, test = "Chi")

# create the framework for the plot
newdata5 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),linkage.density=1)
effects5 <- data.frame(effect(c("Treatment","Season"),model2LDa))

newdata5 <- merge(newdata5,effects5)

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
  xlab(" ")+ ylab("Linkage density")+ 
  scale_y_continuous(limits=c(0,12),oob=squish)+
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
        legend.key = element_rect(size = 5, fill = "white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g5




## robustness
summary(model2RBa)
drop1(model2RBa, test = "Chi")

# create the framework for the plot
newdata6 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),robustness=1)
effects6 <- data.frame(effect(c("Treatment","Season"),model2RBa))

newdata6 <- merge(newdata6,effects6)
newdata6$Treatment <- revalue(newdata6$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g6 <- ggplot(newdata6,
             aes(x=Treatment, y=fit, fill=Season))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  xlab(" ")+ ylab("Robustness")+ 
  scale_y_continuous(limits = c(0,1), oob=squish)+
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
        legend.key = element_rect(size = 5, fill = "white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g6



## vulnerability

summary(model2Va)
drop1(model2Va, test="Chi")

# create the framework for the plot
newdata7 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),vulnerability.LL=1)
effects7 <- data.frame(effect(c("Treatment","Season"),model2Va))

newdata7 <- merge(newdata7,effects7)

newdata7$Treatment <- revalue(newdata7$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g7 <- ggplot(newdata7,
             aes(x=Treatment, y=fit, fill=Season))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  xlab("Season")+ ylab("Generality of plants")+ 
  scale_y_continuous(limits = c(0,20), oob=squish)+
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
        legend.key = element_rect(size = 5, fill = "white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g7



## generality

summary(model2Ga)
drop1(model2Ga, test="Chi")

# create the framework for the plot
newdata8 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),generality.HL=1)
effects8 <- data.frame(effect(c("Treatment","Season"),model2Ga))

newdata8 <- merge(newdata8,effects8)

newdata8$Treatment <- revalue(newdata8$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

# construct the plot
g8 <- ggplot(newdata8,
             aes(x=Treatment, y=fit, fill=Season))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  xlab("Season")+ ylab("Generality of pollinators")+ 
  scale_y_continuous(limits = c(0,10), oob=squish)+
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
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g8



### stitch these all together

m1 <- grid.arrange(g5,g6,g7,g8, ncol=2, nrow=2)

# before we export this we want to move the legend to be outside the plots

# function to extract a legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(g8)

m1a <- grid.arrange(arrangeGrob(g5,g6,
                                g7,g8 + theme(legend.position="none"),
                                nrow=2, ncol=2),
                    mylegend, ncol =2, widths=c(10,2))




ggsave("FigS6.svg", plot = m1a, device = "svg", path = "Results/UpdatedFigs", width = 44, height = 40, units = "cm")


## linkage density

summary(model2LDa)
drop1(model2LDa, test = "Chi")

# create the framework for the plot
newdata9 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),linkage.density=1)
effects9 <- data.frame(effect(c("Treatment","Order"),model2LDa,xlevels=list(Order=1:9)))

newdata9 <- merge(newdata9,effects9)

newdata9$Treatment <- revalue(newdata9$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g9 <- ggplot(newdata9,
             aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab(" ")+ ylab("Linkage density")+ 
  ylim(0,10)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
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



## robustness

# create the framework for the plot
newdata10 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),robustness=1)
effects10 <- data.frame(effect(c("Treatment","Order"),model2RBa,xlevels=list(Order=1:9)))

newdata10 <- merge(newdata10,effects10)
newdata10$Treatment <- revalue(newdata10$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g10 <- ggplot(newdata10,
             aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab(" ")+ ylab("Robustness")+ 
  ylim(0,1)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
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
g10


## vulnerability

# create the framework for the plot
newdata11 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),vulnerability.LL=1)
effects11 <- data.frame(effect(c("Treatment","Order"),model2Va,xlevels=list(Order=1:9)))

newdata11 <- merge(newdata11,effects11)

newdata11$Treatment <- revalue(newdata11$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g11 <- ggplot(newdata11,
              aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab("Seasons since study began")+ ylab("Generality of plants")+ 
  ylim(0,20)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
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
g11




## generality

# create the framework for the plot
newdata12 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),generality.HL=1)
effects12 <- data.frame(effect(c("Treatment","Order"),model2Ga,xlevels=list(Order=1:9)))

newdata12 <- merge(newdata12,effects12)

newdata12$Treatment <- revalue(newdata12$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

# construct the plot
g12 <- ggplot(newdata12,
              aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab("Seasons since study began")+ ylab("Generality of pollinators")+ 
  ylim(0,10)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
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
g12




### stitch these all together

m2 <- grid.arrange(g9,g10,g11,g12, ncol=2, nrow=2)

ggsave("FigS8.svg", plot = m2, device = "svg", path = "Results/UpdatedFigs", width = 44, height = 40, units = "cm")





### interactions

## linkage density

summary(model1LDa)
drop1(model1LDa, test = "Chi")

# create the framework for the plot
newdata1 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),linkage.density=1)
effects1 <- data.frame(effect(c("Treatment","Season"),model1LDa))

newdata1 <- merge(newdata1,effects1)

newdata1$Treatment <- revalue(newdata1$Treatment, c("Fire"="Burned","NoFire"="Unburned"))



# construct the plot
g1 <- ggplot(newdata1,
             aes(x=Treatment, y=fit, fill=Season))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  xlab(" ")+ ylab("Linkage density")+ 
  scale_y_continuous(limits=c(0,20),oob=squish)+
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
        legend.key = element_rect(size = 5, fill = "white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g1






## robustness
summary(model1RBa)
drop1(model1RBa, test = "Chi")

# create the framework for the plot
newdata2 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),robustness=1)
effects2 <- data.frame(effect(c("Treatment","Season"),model1RBa))

newdata2 <- merge(newdata2,effects2)
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
  xlab(" ")+ ylab("Robustness")+ 
  scale_y_continuous(limits = c(0,1), oob=squish)+
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
        legend.key = element_rect(size = 5, fill = "white"),
        legend.key.size = unit(2, 'lines'))

# visualize the plot
g2


## vulnerability

summary(model1Va)
drop1(model1Va, test="Chi")

# create the framework for the plot
newdata3 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),vulnerability.LL=1)
effects3 <- data.frame(effect(c("Treatment","Season"),model1Va))

newdata3 <- merge(newdata3,effects3)

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
  xlab("Season")+ ylab("Generality of plants")+ 
  scale_y_continuous(limits = c(0,30), oob=squish)+
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
        legend.key = element_rect(size = 5, fill = "white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g3



## generality

summary(model1Ga)
drop1(model1Ga, test="Chi")

# create the framework for the plot
newdata4 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),generality.HL=1)
effects4 <- data.frame(effect(c("Treatment","Season"),model1Ga))

newdata4 <- merge(newdata4,effects4)

newdata4$Treatment <- revalue(newdata4$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

# construct the plot
g4 <- ggplot(newdata4,
             aes(x=Treatment, y=fit, fill=Season))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  xlab("Season")+ ylab("Generality of pollinators")+ 
  scale_y_continuous(limits = c(0,10), oob=squish)+
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
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(2, 'lines'))


# visualize the plot
g4




### stitch these all together

m3 <- grid.arrange(g1,g2,g3,g4, ncol=2, nrow=2)

# before we export this we want to move the legend to be outside the plots

# function to extract a legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(g4)

m3a <- grid.arrange(arrangeGrob(g1,g2,
                                g3,g4 + theme(legend.position="none"),
                                nrow=2, ncol=2),
                    mylegend, ncol =2, widths=c(10,2))



ggsave("Fig6.svg", plot = m3a, device = "svg", path = "Results/UpdatedFigs", width = 44, height = 40, units = "cm")


## linkage density

summary(model1LDa)
drop1(model1LDa, test = "Chi")

# create the framework for the plot
newdata13 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),linkage.density=1)
effects13 <- data.frame(effect(c("Treatment","Order"),model1LDa,xlevels=list(Order=1:9)))

newdata13 <- merge(newdata13,effects13)

newdata13$Treatment <- revalue(newdata13$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g13 <- ggplot(newdata13,
             aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab(" ")+ ylab("Linkage density")+ 
  ylim(0,20)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
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
g13



## robustness

# create the framework for the plot
newdata14 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),robustness=1)
effects14 <- data.frame(effect(c("Treatment","Order"),model1RBa,xlevels=list(Order=1:9)))

newdata14 <- merge(newdata14,effects14)
newdata14$Treatment <- revalue(newdata14$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g14 <- ggplot(newdata14,
              aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab(" ")+ ylab("Robustness")+ 
  ylim(0,1)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
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
g14


## vulnerability

# create the framework for the plot
newdata15 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),vulnerability.LL=1)
effects15 <- data.frame(effect(c("Treatment","Order"),model1Va,xlevels=list(Order=1:9)))

newdata15 <- merge(newdata15,effects15)

newdata15$Treatment <- revalue(newdata15$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g15 <- ggplot(newdata15,
              aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab("Seasons since study began")+ ylab("Generality of plants")+ 
  ylim(0,25)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
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
g15




## generality

# create the framework for the plot
newdata16 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),generality.HL=1)
effects16 <- data.frame(effect(c("Treatment","Order"),model1Ga,xlevels=list(Order=1:9)))

newdata16 <- merge(newdata16,effects16)

newdata16$Treatment <- revalue(newdata16$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

# construct the plot
g16 <- ggplot(newdata16,
              aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab("Seasons since study began")+ ylab("Generality of pollinators")+ 
  ylim(0,10)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
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
g16




### stitch these all together

m4 <- grid.arrange(g13,g14,g15,g16, ncol=2, nrow=2)

ggsave("FigS7.svg", plot = m4, device = "svg", path = "Results/UpdatedFigs", width = 44, height = 40, units = "cm")










