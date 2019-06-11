###################################################################################
####   Script for basic insect abundance and assemblage composition analysis   ####
###################################################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS","plyr","reshape2","vegan","scales","effects","svglite","gridExtra")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

# for some reason the package glmmADMB won't install via the usual methods, so:
#install.packages("R2admb")
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),type="source")
library(glmmADMB)



### load up Callum's custom set of functions
f <- c("CheckResidsFunction.R","CheckConvergenceFunction.R","SpecAccumFunctions.R","MultiplotFunction.R")
lapply(f, source)


### Moths

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

### tell R that SampleID and SamplingDay should be treated as factors
dframe1$SamplingDay <- factor(dframe1$SamplingDay)
dframe1$Month <- ordered(dframe1$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
dframe1$Year <- factor(dframe1$Year)
dframe1$Treatment <- relevel(dframe1$Treatment, "NoFire")
dframe1$Date <- as.Date(dframe1$Date, format="%d/%m/%Y")

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


# We don't need the pollen data or the individual SampleID for this analysis, so get rid of it:

dframe1r <- dframe1[,c(2:4,7:12,79)]
summary(dframe1r)

# each row currently contains 1 insect, so add a column for "Count"; 
# this is 1 in every instance except the one sample with zero insects
dframe1r$Count <- ifelse(dframe1r$Family_Species=="none",0,1)

# summarise the dataframe to how many individuals of each species in each sample
dframe2 <- ddply(dframe1r, .(Family_Species,Site,Treatment,Date,Sample,Season,Month,Year,ExactSeason), numcolwise(sum))

# and additionally summarise it into how many insects (of any species) in each sample
dframe2r <- dframe2[,-1]
dframe2r$SpeciesRichness <- ifelse(dframe2r$Count==0,0,1)
dframe3 <- ddply(dframe2r, .(Site,Treatment,Date,Sample,Season,Month,Year,ExactSeason), numcolwise(sum))

# now introduce the season order data

# we also want this variable to be ordered so that we can look for any trend over time
# I've just crudely put this into a text file that we can import and merge in

dframe3 <- merge(SeasonOrder,dframe3)
summary(dframe3$ExactSeason)
summary(dframe3$Order)


# now we're ready to start analysing 
summary(dframe3)

# plot abundance against Treatment and Season
plot(Count ~ Treatment, data = dframe3)
plot(Count ~ Season, data = dframe3)
plot(Count ~ Month, data = dframe3)
plot(Count ~ ExactSeason, data = dframe3)
plot(Count ~ Order, data = dframe3)
plot(Count ~ Date, data = dframe3)

# test effects of Treatment and Season on abundance
# first plot distribution of abundance
hist(dframe3$Count)
summary(dframe3$Count)

# data appear Poisson and possibly overdispersed
# a quick check for overdispersion; in a regular Poisson distribution, mean=variance:
mean(dframe3$Count)
var(dframe3$Count)  # var > mean therefore data are overdispersed


### Data appear overdispersed, so we will try two types of model in addition to Poisson: quasi-Poisson, and negative binomial
### However not many zeroes in data, so zero-inflated models would be inappropriate

# construct models using Date and Site as random effects
# random effects are factors that might affect the output variable, but not in a way that is interesting to us
# you might see variation between different sampling days, or between different fire or non-fire sites...
# ...but we are only really interested in variation due to Treatment


# Poisson model using lme4

model1P <- glmer(Count ~ Treatment * Season + Treatment * Order # fixed effects
                 + (1|Site), # random effects
                 family = poisson (link = "log"),
                 data = dframe3)

# inspect and test the model
summary(model1P)
 
# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model1P, test="Chisq")

chkconv(model1P)

# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you need to pass a model and you can optionally pass up to two explanatory variables as well
# you want residuals that are roughly normally distributed around zero with no obvious trends
# you want residuals that are roughly the same for each level of your treatments
chkres(model1P,dframe3$Treatment,dframe3$Season)  # these residuals are fine

# but that second interaction is non-significant so let's try without

model1Pa <- glmer(Count ~ Treatment * Season + Order
                  + (1|Site),
                  family = poisson (link = "log"),
                  data = dframe3)

chkconv(model1Pa)

summary(model1Pa)
drop1(model1Pa, test = "Chi")

chkres(model1Pa,dframe3$Treatment,dframe3$Season) # these are fine so that's our final model.



### Species Richness ###

# plot richness against Treatment and Season
plot(SpeciesRichness ~ interaction(Treatment,Season), data = dframe3)
plot(SpeciesRichness ~ Order, data = dframe3)


# test effects of Treatment and Season on abundance
# first plot distribution of abundance
hist(dframe3$SpeciesRichness)
hist(log(dframe3$SpeciesRichness+1)) # There are zeroes in this data which can't be log transformed
summary(dframe3$SpeciesRichness)


# data are Poisson and possibly overdispersed
mean(dframe3$SpeciesRichness)
var(dframe3$SpeciesRichness)  # var > mean therefore data are overdispersed


# Poisson model using lme4

model2P <- glmer(SpeciesRichness ~ Treatment * Season + Treatment * Order # fixed effects
                 + (1|Site), # random effects
                 family = poisson (link = "log"),
                 data = dframe3)

# warning about model convergence so check it again, using:

chkconv(model2P) # model is fine

# inspect and test the model
summary(model2P)
drop1(model2P, test="Chisq")  

# check the model's residuals
chkres(model2P,dframe3$Treatment,dframe3$Season)  # these residuals are fine

# therefore let's retry without some of the interactions

model2Pa <- glmer(SpeciesRichness ~ Treatment*Season + Order # fixed effects
                  + (1|Site), # random effects
                  family = poisson (link = "log"),
                  data = dframe3)

chkconv(model2Pa)
  
# inspect and test the model
summary(model2Pa)
drop1(model2Pa, test="Chisq")  

# check the model's residuals
chkres(model2Pa,dframe3$Treatment,dframe3$Season)  # these residuals are fine






### Plants

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe4<-read.table("Data/PlantTransects.txt", header=TRUE)

dframe4$Year <- factor(dframe4$Year)
dframe4$Treatment <- relevel(dframe4$Treatment, "NoFire")
summary(dframe4)

# we need a variable to allow the networks to be paired:

dframe4$ExactSeason <- do.call(paste, c(dframe4[c("Season","Year")], sep = "_"))
dframe4$ExactSeason <- as.factor(dframe4$ExactSeason)
summary(dframe4$ExactSeason)

# We don't need the individual TransectID for this analysis, so get rid of it:

dframe4r <- dframe4[,c(2:15)]
summary(dframe4r)

# each row currently contains 1 plant species, so add a column for "Count"; 
# this is 1 in every instance except the transect with zero plants
dframe4r$SpeciesRichness <- ifelse(dframe4r$PlantSpecies=="none",0,1)

# summarise the dataframe to how many individuals of each species in each sample
dframe5 <- ddply(dframe4r, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year,ExactSeason), numcolwise(mean))

dframe5SR <- ddply(dframe4r, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year,ExactSeason), numcolwise(sum))


dframe5$Month <- ordered(dframe5$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

# merge Order in
dframe5 <- merge(SeasonOrder,dframe5)
summary(dframe5$ExactSeason)
summary(dframe5$Order)

summary(dframe5)


dframe5SR <- merge(SeasonOrder,dframe5SR)
summary(dframe5SR$ExactSeason)
summary(dframe5SR$Order)

summary(dframe5SR)

### finally, I want to add a quick check of the relative roles of different life-histories 
# (annuals, perennials, shrubs (/trees) and bulbs) in any patterns in abundance & species richness

# I have gathered this data in a file, so read that in now
strategies <- read.csv("Data/AllPlantSpecies_strategy.csv", header=T)

# we want to generate equivalents to dframe5 but with the additional level of the different strategies
# first merge the strategy data into dframe4

dframe4.strat <- merge(dframe4, strategies, all=T)

# get rid of individual Transect ID - and also the piggybacking family and no. transects columns


dframe4r.strat <- dframe4.strat[,c(1,3:15,18)]
summary(dframe4r.strat)

# each row currently contains 1 plant species, so add a column for "Count"; 
# this is 1 in every instance except the transect with zero plants
dframe4r.strat$SpeciesRichness <- ifelse(dframe4r$PlantSpecies=="none",0,1)

# summarise the dataframe to how many individuals of each species in each sample
dframe5.strat <- ddply(dframe4r.strat, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year,ExactSeason,Strategy), numcolwise(mean))

dframe5SR.strat <- ddply(dframe4r.strat, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year,ExactSeason,Strategy), numcolwise(sum))


dframe5.strat$Month <- ordered(dframe5.strat$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

# merge Order in
dframe5.strat <- merge(SeasonOrder,dframe5.strat)
summary(dframe5.strat$ExactSeason)
summary(dframe5.strat$Order)

summary(dframe5.strat)


dframe5SR.strat <- merge(SeasonOrder,dframe5SR.strat)
summary(dframe5SR.strat$ExactSeason)
summary(dframe5SR.strat$Order)

summary(dframe5SR.strat)



### now to analyse...

# Floral abundance
plot(PlantCoverage ~ Treatment, dframe5)
plot(PlantCoverage ~ Season, dframe5)
plot(PlantCoverage ~ Month, dframe5)
plot(PlantCoverage ~ Year, dframe5)
plot(PlantCoverage ~ interaction(Treatment,Season), dframe5)
plot(PlantCoverage ~ Order, dframe5)


hist(dframe5$PlantCoverage)
hist(log(dframe5$PlantCoverage+1))
summary(dframe5$PlantCoverage)

# possibly overdispersion, possibly zero-inflation (1st quartile is 0)
# because it's average data it's non-integer so Poisson won't work
# but the log looks close to Gaussian so let's try that first
mean(dframe5$PlantCoverage)
var(dframe5$PlantCoverage)  # var > mean therefore data are overdispersed

# Gaussian

dframe5$tPlantCoverage <- log(dframe5$PlantCoverage + 1)
summary(dframe5$tPlantCoverage)
hist(dframe5$tPlantCoverage)


model5G <- lmer(tPlantCoverage ~ Treatment*Season + Treatment*Order
                 + (1|Site),
                 data = dframe5)

summary(model5G)
drop1(model5G, test = "Chi")

chkres(model5G,dframe5$Treatment,dframe5$Season) # these look ok-ish 
# but have a fairly clear negative trend and problems with the fit of season

# check it without non-sig interaction
model5Ga <- lmer(tPlantCoverage ~ Treatment*Season + Order
                 + (1|Site),
                 data=dframe5)

summary(model5Ga)
drop1(model5Ga, test = "Chi")

chkres(model5Ga,dframe5$Treatment,dframe5$Season) # same issue - mostly pretty good except a negative trend



# Poisson
# can't deal with non-integers so we'll have to use total rather than mean coverage

model5P <- glmer(PlantCoverage ~ Treatment*Season + Treatment*Order
                 + (1|Site),
                 family = poisson (link = "log"),
                 data = dframe5SR)

chkconv(model5P)

summary(model5P)
drop1(model5P, test = "Chi")

chkres(model5P,dframe5SR$Treatment,dframe5SR$Season) # these are worse than the Gaussian

# nbinom

#model5NB <- glmer.nb(PlantCoverage ~ Treatment*Season + Treatment*Order
#                  + (1|Site),
#                  data = dframe5SR)

#chkconv(model5NB)

#summary(model5NB)
#drop1(model5NB, test = "Chi")
#
#chkres(model5NB,dframe5SR$Treatment,dframe5SR$Season) # these are the worst yet

# we need to bite the bullet and use zero-inflated models
# zero-inflated poisson

#model5ZIP <- glmmadmb(PlantCoverage ~  Treatment*Season + Treatment*Order
#                      + (1|Site), #Random effects
#                      zeroInflation=TRUE,
#                      family = "poisson",
#                      link = "log",
#                      data = dframe5)

#summary(model5ZIP)

#Anova(model5ZIP, type="III")
#drop1(model5ZIP, test="Chisq")
#
#chkres.zi(model5ZIP,dframe5SR$Treatment,dframe5SR$Season) # these are pretty horrible


# zero-inflated nbinom

#model5ZINB <- glmmadmb(PlantCoverage ~  Treatment*Season + Treatment* Order
#                       + (1|Site), #Random effects
#                       zeroInflation=TRUE,
#                       family = "nbinom",
#                       link = "log",
#                       data = dframe5SR)

#summary(model5ZINB)

#Anova(model5ZINB, type = "III")
#drop1(model5ZINB, test="Chisq")

#chkres.zi(model5ZINB,dframe5SR$Treatment,dframe5SR$Season)



# actually everything is substantially worse than that log-Gaussian model at the start
# let's just use that one.

summary(model5Ga)
drop1(model5Ga, test = "Chi")




# Floral Species Richness
plot(SpeciesRichness ~ Treatment, dframe5SR)
plot(SpeciesRichness ~ Season, dframe5SR)
plot(SpeciesRichness ~ Month, dframe5SR)
plot(SpeciesRichness ~ interaction(Treatment,Season), dframe5SR)
plot(SpeciesRichness ~ Order, dframe5SR)

hist(dframe5SR$SpeciesRichness)

# possibly overdispersion (but not sure...)
mean(dframe5SR$SpeciesRichness)
var(dframe5SR$SpeciesRichness)  # var > mean therefore data are overdispersed but not too badly here


# Poisson model
model6P <- glmer(SpeciesRichness ~ Treatment*Season+Treatment*Order
                 + (1|Site),
                 family = poisson (link="log"),
                 data = dframe5SR)

chkconv(model6P)

summary(model6P)
drop1(model6P, test = "Chi")

chkres(model6P,dframe5SR$Treatment,dframe5SR$Season) # these aren't fabulous

# try removing non-significant interaction

model6Pa <- glmer(SpeciesRichness ~ Treatment*Season + Order
                  + (1|Site),
                  family = poisson (link="log"),
                  data = dframe5SR)

summary(model6Pa)
drop1(model6Pa)

chkres(model6Pa,dframe5SR$Treatment,dframe5SR$Season) # these aren't much better!


# negative binomial
#model6NB <- glmer.nb(SpeciesRichness ~ Treatment*Season + Treatment*Order
#                     + (1|Site),
#                     data = dframe5SR)

#chkconv(model6NB)

#summary(model6NB)
#drop1(model6NB, test="Chi")

#chkres(model6NB,dframe5SR$Treatment,dframe5SR$Season)

# non-sig interaction
#model6NBa <- glmer.nb(SpeciesRichness ~ Treatment*Season + Order
#                      + (1|Site),
#                      data = dframe5SR)

#chkconv(model6NBa)

#summary(model6NBa)
#drop1(model6NBa, test="Chi")

#chkres(model6NBa,dframe5SR$Treatment,dframe5SR$Season) # still no better!



# try a log-transformed gaussian as that worked previously
dframe5SR$tSpeciesRichness <- log(dframe5SR$SpeciesRichness+1)
hist(dframe5SR$tSpeciesRichness)

model6G <- lmer(tSpeciesRichness ~ Treatment*Season + Treatment*Order
                + (1|Site),
                data = dframe5SR)

summary(model6G)
drop1(model6G, test = "Chi")

chkres(model6G,dframe5SR$Treatment,dframe5SR$Season) # these are not terrible, though the zeroes in summer continue to have a visible effect

# Gaussian is the best yet anyway



### Figures
# we want a figure each from moth abundance, moth estimated SR, plant abundance and plant estimated SR

# moth abundance
summary(model1Pa)
drop1(model1Pa, test="Chi")

# create the framework for the plot
newdata1 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),Order=9,Count=1)
# predict the model outputs
effects <- data.frame(effect(c("Treatment","Season"),model1Pa))
newdata1 <- merge(newdata1,effects)
newdata1$Treatment <- revalue(newdata1$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g1 <- ggplot(newdata1,
             aes(x=Treatment, y=fit, fill=Season))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  ggtitle("Moths")+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  xlab(" ")+ ylab("Abundance")+ 
  scale_y_continuous(limits = c(0,200), oob=squish)+
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
        text=element_text(size=30,color="black"),
        axis.text.y=element_text(color="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5))


# visualize the plot
g1



## repeat this for the other three


# extrapolated moth SR
# I've already prepped the data here in a different script where I carried out this analysis
newdata2 <- expand.grid(Treatment=c("Burned","Unburned"),Season=c("Spring","Summer","Autumn","Winter"),SR=1)
effects <- read.table("Results/mothSRmodel.txt", header = T)
newdata2 <- merge(newdata2, effects)


# construct the plot
g2 <- ggplot(newdata2,
             aes(x=Treatment, y=tfit, fill=Season))+
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  scale_y_continuous(limits = c(0,100), oob=squish)+
  xlab("Season")+ ylab("Species richness")+ 
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
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5))


# visualize the plot
g2



# plant abundance
summary(model5Ga)
drop1(model5Ga, test = "Chi")

# create the framework for the plot
newdata3 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),tPlantCoverage=1)
effects <- data.frame(effect(c("Treatment","Season"),model5Ga))
newdata3 <- merge(newdata3,effects)

newdata3$tfit <- (exp(newdata3$fit))-1
newdata3$tlower <- (exp(newdata3$lower))-1
newdata3$tupper <- (exp(newdata3$upper))-1


newdata3$Treatment <- revalue(newdata3$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g3 <- ggplot(newdata3,
             aes(x=Treatment, y=tfit, fill=Season))+
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  scale_y_continuous(limits=c(0,15),oob=squish)+
  xlab(" ")+ ylab(" ")+
  facet_grid(.~Season, switch="x", space = "free_x", scales = "free_x")+
  guides(fill=F)+
  ggtitle("Flowers")+
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
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5))


# visualize the plot
g3





# plant SR - again I've prepared the data externally
newdata4 <- expand.grid(Treatment=c("Burned","Unburned"),Season=c("Spring","Summer","Autumn","Winter"),tPlantCoverage=1)
effects <- read.table("Results/plantSRmodel.txt", header = T)
newdata4 <- merge(newdata4, effects)

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
  scale_y_continuous(limits = c(0,8), oob=squish)+
  xlab("Season")+ ylab(" ")+ 
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



# now stitch them into a multiplot

m1 <- grid.arrange(g1,g3,g2,g4, ncol=2, nrow=2)

# before we export this we want to move the legend to be outside the plots

# function to extract a legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(g4)

m1a <- grid.arrange(arrangeGrob(g1,g3,
                                g2,g4 + theme(legend.position="none"),
                                nrow=2, ncol=2),
                    mylegend, ncol =2, widths=c(10,2))






#ggsave("Fig3.svg", plot = m1a, device = "svg", path = "Results/UpdatedFigs", width = 44, height = 40, units = "cm")






### now we want to repeat those figures but for order

# we want a figure each from moth abundance, moth SR, plant abundance and plant SR

# moth abundance
summary(model1Pa)
drop1(model1Pa, test="Chi")

# create the framework for the plot
newdata5 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),Count=1)
# predict the model outputs
effects5 <- data.frame(effect(c("Treatment","Order"),model1Pa,xlevels=list(Order=1:9)))
newdata5 <- merge(newdata5,effects5)
newdata5$Treatment <- revalue(newdata5$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g5 <- ggplot(newdata5,
             aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab(" ")+ ylab("Abundance")+ 
  scale_y_continuous(limits = c(0,150), oob=squish)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
  ggtitle("Moths")+
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
g5




## repeat this for the other three


# moth SR
# create the framework for the plot
newdata6 <- expand.grid(Treatment=c("Burned","Unburned"),Order=c(1:9),Count=1)
# predict the model outputs
effects6 <- read.table("Results/mothSRmodelorder.txt", header = T)
newdata6 <- merge(newdata6,effects6)
newdata6


# construct the plot
g6 <- ggplot(newdata6,
             aes(x=Order, y=tfit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=tlower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=tupper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab("Seasons since\nstudy began")+ ylab("Species richness")+ 
  scale_y_continuous(limits = c(0,60), oob=squish)+
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
g6





# plant abundance
summary(model5Ga)
drop1(model5Ga, test = "Chi")

# create the framework for the plot
newdata7 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),Count=1)
# predict the model outputs
effects7 <- data.frame(effect(c("Treatment","Order"),model5G,xlevels=list(Order=1:9)))
newdata7 <- merge(newdata7,effects7)
newdata7$Treatment <- revalue(newdata7$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

newdata7$tfit <- (exp(newdata7$fit))-1
newdata7$tlower <- (exp(newdata7$lower))-1
newdata7$tupper <- (exp(newdata7$upper))-1



# construct the plot
g7 <- ggplot(newdata7,
             aes(x=Order, y=tfit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=tlower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=tupper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab(" ")+ ylab(" ")+ 
  scale_y_continuous(limits = c(0,10), oob=squish)+
  scale_x_continuous(limits=c(1,9),breaks=c(1:9))+
  guides(fill=F)+
  ggtitle("Flowers")+
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





# plant SR
# create the framework for the plot
newdata8 <- expand.grid(Treatment=c("Burned","Unburned"),Order=c(1:9),Count=1)
# predict the model outputs
effects8 <- read.table("Results/plantSRmodelorder.txt", header = T)
newdata8 <- merge(newdata8,effects8)

# construct the plot
g8 <- ggplot(newdata8,
             aes(x=Order, y=fit, group=Treatment))+
  geom_line(aes(linetype=Treatment),
            colour="black",stat="identity",size=2)+
  geom_line(aes(x=Order,y=lower,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  geom_line(aes(x=Order,y=upper,linetype=Treatment),
            colour="black",stat="identity",size=0.75)+
  xlab("Seasons since\nstudy began")+ ylab(" ")+ 
  scale_y_continuous(limits = c(0,6), oob=squish)+
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
g8



# now stitch them into a multiplot

m2 <- grid.arrange(g5,g7,g6,g8, ncol=2, nrow=2)

ggsave("FigS2.svg", plot = m2, device = "svg", path = "Results/", width = 44, height = 40, units = "cm")




### Arbutus unedo/Ulex argenteus

# A. unedo reported by Paula to be the dominant pollen source in autumn/winter. 
# From my examination of data (below), this appears to be true for autumn, but overtaken by U. argenteus in winter
# We want to see whether their floral abundance differs substantially between burned and unburned sites,
# as this could help explain the observed patterns

# recall dframe4
summary(dframe4)

arbutus <- dframe4[which(dframe4$PlantSpecies=="Arbutus unedo"),]
summary(arbutus)

# evidently there is almost no data on this species from floral transects
# try Ulex

ulex <- dframe4[which(dframe4$PlantSpecies=="Ulex argenteus"),]
summary(ulex)

# that's more like it! Let's check if there's any effect of burning on winter-time abundance of U. australis flowers...

ulex <- ulex[which(ulex$Season=="Winter"),]

# merge Order in
ulex <- merge(SeasonOrder,ulex)
summary(ulex$ExactSeason)
summary(ulex$Order)

summary(ulex)


plot(PlantCoverage ~ Treatment, ulex)
plot(PlantCoverage ~ Order, ulex)

# Season doesn't exist

modulex <- glmer(PlantCoverage ~ Treatment * Order
                 + (1|Site),
                 family=poisson (link="log"),
                 data = ulex)

chkconv(modulex)

summary(modulex)
drop1(modulex, test = "Chi")

chkres(modulex, ulex$Treatment, ulex$Season) # these are fine


# and repeat the analysis with everything but Ulex
noulex <- dframe4[which(dframe4$PlantSpecies!="Ulex argenteus"),]
summary(noulex)

# that's more like it! Let's check if there's any effect of burning on winter-time abundance of U. australis flowers...

noulex <- noulex[which(noulex$Season=="Winter"),]

# merge Order in
noulex <- merge(SeasonOrder,noulex)
summary(noulex$ExactSeason)
summary(noulex$Order)

summary(noulex)

# for this we need to collapse everything by sample
# each row currently contains 1 plant species, so add a column for "Count"; 
# this is 1 in every instance except the transect with zero plants
noulex$SpeciesRichness <- ifelse(noulex$PlantSpecies=="none",0,1)

# summarise the dataframe to how many individuals of each species in each sample
noulex1 <- ddply(noulex, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year,ExactSeason), numcolwise(mean))

plot(PlantCoverage ~ Treatment, noulex1)
plot(PlantCoverage ~ Order, noulex1)

# Season doesn't exist

modnoulex <- glmer(PlantCoverage ~ Treatment * Order
                 + (1|Site),
                 family=poisson (link="log"),
                 data = noulex)

chkconv(modnoulex)

summary(modnoulex)
drop1(modnoulex, test = "Chi")

chkres(modnoulex, noulex$Treatment, noulex$Season) # these aren't great - negative trend

# before we used a Gaussian + transformation for this (with Ulex included) so let's replicate that now it's not

# Gaussian

noulex$tPlantCoverage <- log(noulex$PlantCoverage + 1)
summary(noulex$tPlantCoverage)
hist(noulex$tPlantCoverage)


modnoulexG <- lmer(tPlantCoverage ~ Treatment*Order
                + (1|Site),
                data = noulex)

summary(modnoulexG)
drop1(modnoulexG, test = "Chi")

chkres(modnoulexG,noulex$Treatment,noulex$Season) # these look a bit better 

# check it without non-sig interaction
modnoulexGa <- lmer(tPlantCoverage ~ Treatment + Order
                   + (1|Site),
                   data = noulex)

summary(modnoulexGa)
drop1(modnoulexGa, test = "Chi")

chkres(modnoulexGa,noulex$Treatment,noulex$Season) # these look a bit better 



# let's also look into the patterns of its pollen.
# first let's check what the dominant plant species were in all seasons, just in case we want to check them too
winter <- dframe4[which(dframe4$Season=="Winter"),]
autumn <- dframe4[which(dframe4$Season=="Autumn"),]
spring <- dframe4[which(dframe4$Season=="Spring"),]
summer <- dframe4[which(dframe4$Season=="Summer"),]

summary(winter)
summary(autumn)
summary(spring)
summary(summer)

# Ulex in winter is interesting and we may come back to it

# look instead at the pollen data
summary(dframe1)

# check dominant pollen species
poll_autumn <- dframe1[which(dframe1$Season=="Autumn"),]
poll_winter <- dframe1[which(dframe1$Season=="Winter"),]
poll_spring <- dframe1[which(dframe1$Season=="Spring"),]
poll_summer <- dframe1[which(dframe1$Season=="Summer"),]


summary(poll_autumn) # Arbutus, Ulex, Solanum
summary(poll_winter) # Ulex, possibly Gladiolus, not Arbutus
summary(poll_spring) # diverse, but nothing really stands out (Anagallis? has a 3rd quartile!)
summary(poll_summer) # not much at all


# let's return focus to Ulex australis in winter
pollen_ulex <- poll_winter[,c(1:12,71,79)]
pollen_no_ulex <- poll_winter[,-71]


### total up the pollen grains and types for each insect
pollen_ulex$PollenLoad <- rowSums(pollen_ulex[c(13:(length(pollen_ulex)-1))])
pollen_no_ulex$PollenLoad <- rowSums(pollen_no_ulex[c(13:(length(pollen_no_ulex)-1))])

pollen_ulex$PollenTypes <- rowSums(pollen_ulex[c(13:(length(pollen_ulex)-2))] != 0)  # shouldn't produce anything larger than 1!
pollen_no_ulex$PollenTypes <- rowSums(pollen_no_ulex[c(13:(length(pollen_no_ulex)-2))] != 0)

### create a binary (yes/no) variable for whether each insect is carrying any pollen
pollen_ulex$PollenYN <- ifelse(pollen_ulex$PollenTypes==0,0,1)
pollen_no_ulex$PollenYN <- ifelse(pollen_no_ulex$PollenTypes==0,0,1)

summary(pollen_ulex)
summary(pollen_no_ulex)


# now introduce the season order data

# we also want this variable to be ordered so that we can look for any trend over time
# I've just crudely put this into a text file that we can import and merge in

pollen_ulex <- merge(SeasonOrder,pollen_ulex)
summary(pollen_ulex$ExactSeason)

pollen_no_ulex <- merge(SeasonOrder,pollen_no_ulex)
summary(pollen_no_ulex$ExactSeason)

# now we're ready to start analysing 


# likelihood of pollen transport
# first Ulex only
plot(pollen_ulex$PollenYN ~ pollen_ulex$Treatment)
plot(pollen_ulex$PollenYN ~ pollen_ulex$Season)

modelulexYN <- glmer(PollenYN ~ Treatment * Order
                   + (1|Site),
                   family = binomial (link = "logit"),
                   data = pollen_ulex)

chkconv(modelulexYN)

summary(modelulexYN)
drop1(modelulexYN, test = "Chi")

# non-convergent but let's try dropping the non-sig interaction

modelulexYNa <- glmer(PollenYN ~ Treatment + Order
                     + (1|Site),
                     family = binomial (link = "logit"),
                     data = pollen_ulex)

chkconv(modelulexYNa)

summary(modelulexYNa)
drop1(modelulexYNa, test = "Chi")


# and non-Ulex only
plot(pollen_no_ulex$PollenYN ~ pollen_no_ulex$Treatment)
plot(pollen_no_ulex$PollenYN ~ pollen_no_ulex$Season)

modelnoulexYN <- glmer(PollenYN ~ Treatment * Order
                   + (1|Site),
                   family = binomial (link = "logit"),
                   data = pollen_no_ulex)

chkconv(modelnoulexYN)

summary(modelnoulexYN)
drop1(modelnoulexYN, test = "Chi")

# drop non-sig interaction

modelnoulexYNa <- glmer(PollenYN ~ Treatment + Order
                       + (1|Site),
                       family = binomial (link = "logit"),
                       data = pollen_no_ulex)

chkconv(modelnoulexYNa)

summary(modelnoulexYNa)
drop1(modelnoulexYNa, test = "Chi")





# pollen load per indiv

# first Ulex only
pollen_ulexP <- pollen_ulex[which(pollen_ulex$PollenYN!=0),]
pollen_no_ulexP <- pollen_no_ulex[which(pollen_no_ulex$PollenYN!=0),]


plot(pollen_ulexP$PollenLoad ~ pollen_ulexP$Treatment)
plot(pollen_ulexP$PollenLoad ~ pollen_ulexP$Order)

hist(pollen_ulexP$PollenLoad)

# there's an obvious outlier here which I'm going to remove for safety's sake

pollen_ulexPa <- pollen_ulexP[which(pollen_ulexP$PollenLoad<1000),]

plot(pollen_ulexPa$PollenLoad ~ pollen_ulexPa$Treatment)
plot(pollen_ulexPa$PollenLoad ~ pollen_ulexPa$Order)

hist(pollen_ulexPa$PollenLoad)

modelulexLoad <- lmer(log(PollenLoad) ~ Treatment * Order
                      + (1|Site),
                      data = pollen_ulexPa)

chkconv(modelulexLoad)

summary(modelulexLoad)
drop1(modelulexLoad, test = "Chi")

# overfitted but let's try dropping the non-sig interaction

modelulexLoada <- lmer(log(PollenLoad) ~ Treatment + Order
                      + (1|Site),
                      data = pollen_ulexPa)

chkconv(modelulexLoada)

summary(modelulexLoada)
drop1(modelulexLoada, test = "Chi")


# and non-Ulex only
plot(pollen_no_ulexP$PollenLoad ~ pollen_no_ulexP$Treatment)
plot(pollen_no_ulexP$PollenLoad ~ pollen_no_ulexP$Order)

hist(pollen_no_ulex$PollenLoad)

modelno_ulexLoad <- lmer(log(PollenLoad) ~ Treatment * Order
                      + (1|Site),
                      data = pollen_no_ulexP)

chkconv(modelno_ulexLoad)

summary(modelno_ulexLoad)
drop1(modelno_ulexLoad, test = "Chi")


# finally we want to look at total pollen count per sample, for which we'll need to wrangle the data a bit further
pollen_ulexS <- ddply(pollen_ulex, .(Treatment,Sample,Order,Site), numcolwise(sum))
pollen_no_ulexS <- ddply(pollen_no_ulex, .(Treatment,Sample,Order,Site), numcolwise(sum))


# first Ulex only
plot(pollen_ulexS$PollenLoad ~ pollen_ulexS$Treatment)
plot(pollen_ulexS$PollenLoad ~ pollen_ulexS$Order)

hist(pollen_ulexS$PollenLoad)

modelulexSLoad <- lmer(log(PollenLoad+1) ~ Treatment * Order
                      + (1|Site),
                      data = pollen_ulexS)

chkconv(modelulexSLoad)

summary(modelulexSLoad)
drop1(modelulexSLoad, test = "Chi")

# try dropping the non-sig interaction

modelulexSLoada <- lmer(log(PollenLoad+1) ~ Treatment + Order
                       + (1|Site),
                       data = pollen_ulexS)

chkconv(modelulexSLoada)

summary(modelulexSLoada)
drop1(modelulexSLoada, test = "Chi")


# and non-Ulex only
plot(pollen_no_ulexS$PollenLoad ~ pollen_no_ulexS$Treatment)
plot(pollen_no_ulexS$PollenLoad ~ pollen_no_ulexS$Order)

hist(pollen_no_ulexS$PollenLoad)

modelno_ulexSLoad <- lmer(log(PollenLoad+1) ~ Treatment * Order
                         + (1|Site),
                         data = pollen_no_ulexS)

chkconv(modelno_ulexSLoad)

summary(modelno_ulexSLoad)
drop1(modelno_ulexSLoad, test = "Chi")


# drop non-sig interaction
modelno_ulexSLoada <- lmer(log(PollenLoad+1) ~ Treatment + Order
                          + (1|Site),
                          data = pollen_no_ulexS)

chkconv(modelno_ulexSLoada)

summary(modelno_ulexSLoada)
drop1(modelno_ulexSLoada, test = "Chi")



### figures
# plot some figures based on these analyses

# floral abundance
# Ulex only
summary(modulex)

newdata9 <- expand.grid(Treatment=c("Fire","NoFire"),tPlantCoverage=1)
effects9 <- data.frame(effect(c("Treatment"),modulex))
newdata9 <- merge(newdata9,effects9)

newdata9$Treatment <- revalue(newdata9$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g9 <- ggplot(newdata9,
             aes(x=Treatment, y=fit))+
  geom_point(cex=3,position=position_dodge(width=0.5))+
  scale_y_continuous(limits = c(0,30), oob=squish)+
  xlab(" ")+ ylab("Floral abundance")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  guides(fill=F)+
  ggtitle(expression(paste(italic("Ulex australis "))))+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.x=element_text(size=25, color="black"),
        axis.text.y=element_text(size=25, color="black"))


# visualize the plot
g9



# non-Ulex only
summary(modnoulex)

newdata10 <- expand.grid(Treatment=c("Fire","NoFire"),tPlantCoverage=1)
effects10 <- data.frame(effect(c("Treatment"),modnoulex))
newdata10 <- merge(newdata10,effects10)

newdata10$Treatment <- revalue(newdata10$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g10 <- ggplot(newdata10,
             aes(x=Treatment, y=fit))+
  geom_point(cex=3,position=position_dodge(width=0.5))+
  scale_y_continuous(limits = c(0,30), oob=squish)+
  xlab(" ")+ ylab(" ")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  guides(fill=F)+
  ggtitle(expression(paste(italic("Ulex australis "), "  excluded")))+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.x=element_text(size=25, color="black"),
        axis.text.y=element_text(size=25, color="black"))


# visualize the plot
g10



# proportion of pollen carriers
# Ulex only
summary(modelulexYNa)

newdata11 <- expand.grid(Treatment=c("Fire","NoFire"),tPlantCoverage=1)
effects11 <- data.frame(effect(c("Treatment"),modelulexYNa))
newdata11 <- merge(newdata11,effects11)

newdata11$Treatment <- revalue(newdata11$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g11 <- ggplot(newdata11,
             aes(x=Treatment, y=fit))+
  geom_point(cex=3,position=position_dodge(width=0.5))+
  scale_y_continuous(limits = c(0,1), oob=squish)+
  xlab(" ")+ ylab("Probability of detecting pollen")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  guides(fill=F)+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.x=element_text(size=25, color="black"),
        axis.text.y=element_text(size=25, color="black"))


# visualize the plot
g11



# non-Ulex only
summary(modelnoulexYNa)

newdata12 <- expand.grid(Treatment=c("Fire","NoFire"),tPlantCoverage=1)
effects12 <- data.frame(effect(c("Treatment"),modelnoulexYNa))
newdata12 <- merge(newdata12,effects12)

newdata12$Treatment <- revalue(newdata12$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g12 <- ggplot(newdata12,
              aes(x=Treatment, y=fit))+
  geom_point(cex=3,position=position_dodge(width=0.5))+
  scale_y_continuous(limits = c(0,1), oob=squish)+
  xlab(" ")+ ylab(" ")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  guides(fill=F)+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.x=element_text(size=25, color="black"),
        axis.text.y=element_text(size=25, color="black"))


# visualize the plot
g12



# individual pollen load
# Ulex only
summary(modelulexLoada)

newdata13 <- expand.grid(Treatment=c("Fire","NoFire"),tPlantCoverage=1)
effects13 <- data.frame(effect(c("Treatment"),modelulexLoada))
newdata13 <- merge(newdata13,effects13)

newdata13$tfit <- exp(newdata13$fit)
newdata13$tlower <- exp(newdata13$lower)
newdata13$tupper <- exp(newdata13$upper)

newdata13$Treatment <- revalue(newdata13$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g13 <- ggplot(newdata13,
              aes(x=Treatment, y=tfit))+
  geom_point(cex=3,position=position_dodge(width=0.5))+
  scale_y_continuous(limits = c(0,150), oob=squish)+
  xlab(" ")+ ylab("Pollen load per individual")+ 
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  guides(fill=F)+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.x=element_text(size=25, color="black"),
        axis.text.y=element_text(size=25, color="black"))


# visualize the plot
g13



# non-Ulex only
summary(modelno_ulexLoad)

newdata14 <- expand.grid(Treatment=c("Fire","NoFire"),tPlantCoverage=1)
effects14 <- data.frame(effect(c("Treatment"),modelno_ulexLoad))
newdata14 <- merge(newdata14,effects14)

newdata14$tfit <- exp(newdata14$fit)
newdata14$tlower <- exp(newdata14$lower)
newdata14$tupper <- exp(newdata14$upper)

newdata14$Treatment <- revalue(newdata14$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g14 <- ggplot(newdata14,
              aes(x=Treatment, y=tfit))+
  geom_point(cex=3,position=position_dodge(width=0.5))+
  scale_y_continuous(limits = c(0,150), oob=squish)+
  xlab(" ")+ ylab(" ")+ 
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  guides(fill=F)+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.x=element_text(size=25, color="black"),
        axis.text.y=element_text(size=25, color="black"))


# visualize the plot
g14






# sample pollen load
# Ulex only
summary(modelulexSLoada)

newdata15 <- expand.grid(Treatment=c("Fire","NoFire"),tPlantCoverage=1)
effects15 <- data.frame(effect(c("Treatment"),modelulexSLoada))
newdata15 <- merge(newdata15,effects15)

newdata15$tfit <- exp(newdata15$fit)
newdata15$tlower <- exp(newdata15$lower)
newdata15$tupper <- exp(newdata15$upper)

newdata15$Treatment <- revalue(newdata15$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g15 <- ggplot(newdata15,
              aes(x=Treatment, y=tfit))+
  geom_point(cex=3,position=position_dodge(width=0.5))+
  scale_y_continuous(limits = c(0,100), oob=squish)+
  xlab("Treatment")+ ylab("Pollen load per sample")+ 
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  guides(fill=F)+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.x=element_text(size=25, color="black"),
        axis.text.y=element_text(size=25, color="black"))


# visualize the plot
g15



# non-Ulex only
summary(modelno_ulexSLoada)

newdata16 <- expand.grid(Treatment=c("Fire","NoFire"),tPlantCoverage=1)
effects16 <- data.frame(effect(c("Treatment"),modelno_ulexSLoada))
newdata16 <- merge(newdata16,effects16)

newdata16$tfit <- exp(newdata16$fit)
newdata16$tlower <- exp(newdata16$lower)
newdata16$tupper <- exp(newdata16$upper)

newdata16$Treatment <- revalue(newdata16$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g16 <- ggplot(newdata16,
              aes(x=Treatment, y=tfit))+
  geom_point(cex=3,position=position_dodge(width=0.5))+
  scale_y_continuous(limits = c(0,100), oob=squish)+
  xlab("Treatment")+ ylab(" ")+ 
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  guides(fill=F)+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text.x=element_text(size=25, color="black"),
        axis.text.y=element_text(size=25, color="black"))


# visualize the plot
g16




# now stitch them into a multiplot

m2 <- grid.arrange(g9,g10,g11,g12,g13,g14,g15,g16, ncol=2, nrow=4)

#ggsave("FigS6.svg", plot = m2, device = "svg", path = "Results/", width = 40, height = 70, units = "cm")






#### different strategies

# we know know that there is an overall increase in plant coverage and richness at burned sites - i.e.

plot(PlantCoverage ~ Treatment, dframe5)
plot(SpeciesRichness ~ Treatment, dframe5SR)

# but does this vary for each type of plant? first a quick(ish) visual check

plot(PlantCoverage ~ Treatment, dframe5.strat[which(dframe5.strat$Strategy=="Annual"),])
plot(PlantCoverage ~ Treatment, dframe5.strat[which(dframe5.strat$Strategy=="Perennial"),])
plot(PlantCoverage ~ Treatment, dframe5.strat[which(dframe5.strat$Strategy=="Bulb"),])
plot(PlantCoverage ~ Treatment, dframe5.strat[which(dframe5.strat$Strategy=="Tree/shrub"),])


plot(SpeciesRichness ~ Treatment, dframe5SR.strat[which(dframe5SR.strat$Strategy=="Annual"),])
plot(SpeciesRichness ~ Treatment, dframe5SR.strat[which(dframe5SR.strat$Strategy=="Perennial"),])
plot(SpeciesRichness ~ Treatment, dframe5SR.strat[which(dframe5SR.strat$Strategy=="Bulb"),])
plot(SpeciesRichness ~ Treatment, dframe5SR.strat[which(dframe5SR.strat$Strategy=="Tree/shrub"),])

# it's quite clear that *annuals* drive the patterns in plant species richness


# let's try to look at this analytically, with the same re-analysis as we've just done for Ulex

# first generate some sub-datasets

annuals <- dframe5.strat[which(dframe5.strat$Strategy=="Annual"),]
nonannuals <- dframe5.strat[which(!(dframe5.strat$Strategy=="Annual")),]

summary(annuals)
summary(nonannuals)

annualsSR <- dframe5SR.strat[which(dframe5SR.strat$Strategy=="Annual"),]
nonannualsSR <- dframe5SR.strat[which(!(dframe5SR.strat$Strategy=="Annual")),]

summary(annualsSR)
summary(nonannualsSR)


# try some analyses - first the annuals


plot(PlantCoverage ~ Treatment, annuals)
plot(PlantCoverage ~ Season, annuals)
plot(PlantCoverage ~ Order, annuals)



hist(annuals$PlantCoverage)
hist(log(annuals$PlantCoverage+1))
summary(annuals$PlantCoverage)

# possibly overdispersion
# because it's average data it's non-integer so Poisson won't work
# but the log looks close to Gaussian so let's try that first
mean(annuals$PlantCoverage)
var(annuals$PlantCoverage)  # var > mean therefore data are overdispersed

# Gaussian

annuals$tPlantCoverage <- log(annuals$PlantCoverage + 1)
summary(annuals$tPlantCoverage)
hist(annuals$tPlantCoverage)


modannuals <- lmer(tPlantCoverage ~ Treatment*Season + Treatment*Order
                + (1|Site),
                data = annuals)

summary(modannuals)
drop1(modannuals, test = "Chi")

chkres(modannuals,annuals$Treatment,annuals$Season) # these look ok-ish 

# check it without non-sig interaction
modannualsa <- lmer(tPlantCoverage ~ Treatment*Season + Order
                 + (1|Site),
                 data=annuals)

summary(modannualsa)
drop1(modannualsa, test = "Chi")

chkres(modannualsa,annuals$Treatment,annuals$Season) # mostly pretty good



# try Poisson
# can't deal with non-integers so we'll have to use total rather than mean coverage

modannualsP <- glmer(PlantCoverage ~ Treatment*Season + Treatment*Order
                 + (1|Site),
                 family = poisson (link = "log"),
                 data = annualsSR)

chkconv(modannualsP)

summary(modannualsP)
drop1(modannualsP, test = "Chi")

chkres(modannualsP,annualsSR$Treatment,annualsSR$Season) # these are worse than the Gaussian

summary(modannualsa)
drop1(modannualsa, test = "Chi")


# and repeat the analysis with everything but annuals


plot(PlantCoverage ~ Treatment, nonannuals)
plot(PlantCoverage ~ Season, nonannuals)
plot(PlantCoverage ~ Order, nonannuals)



hist(nonannuals$PlantCoverage)
hist(log(nonannuals$PlantCoverage+1))
summary(nonannuals$PlantCoverage)

# possibly overdispersion
# because it's average data it's non-integer so Poisson won't work
# but the log looks close to Gaussian so let's try that first
mean(nonannuals$PlantCoverage)
var(nonannuals$PlantCoverage)  # var > mean therefore data are overdispersed

# Gaussian

nonannuals$tPlantCoverage <- log(nonannuals$PlantCoverage + 1)
summary(nonannuals$tPlantCoverage)
hist(nonannuals$tPlantCoverage)


modnonannuals <- lmer(tPlantCoverage ~ Treatment*Season + Treatment*Order
                   + (1|Site),
                   data = nonannuals)

summary(modnonannuals)
drop1(modnonannuals, test = "Chi")

chkres(modnonannuals,nonannuals$Treatment,nonannuals$Season) # these look fine 

# check it without non-sig interaction
modnonannualsa <- lmer(tPlantCoverage ~ Treatment*Season + Order
                    + (1|Site),
                    data=nonannuals)

summary(modnonannualsa)
drop1(modnonannualsa, test = "Chi")

chkres(modnonannualsa,nonannuals$Treatment,nonannuals$Season) # mostly pretty good


### now let's check species richness as well

plot(SpeciesRichness ~ Treatment, annualsSR)
plot(SpeciesRichness ~ Season, annualsSR)
plot(SpeciesRichness ~ Order, annualsSR)



# try a log-transformed gaussian as that worked previously
annualsSR$tSpeciesRichness <- log(annualsSR$SpeciesRichness+1)
hist(annualsSR$tSpeciesRichness)

modannualsSR <- lmer(tSpeciesRichness ~ Treatment*Season + Treatment*Order
                + (1|Site),
                data = annualsSR)

summary(modannualsSR)
drop1(modannualsSR, test = "Chi")

chkres(modannualsSR,annualsSR$Treatment,annualsSR$Season) # these are not terrible


# try without the non-sig interaction

modannualsSRa <- lmer(tSpeciesRichness ~ Season + Treatment*Order
                      + (1|Site),
                      data = annualsSR)

summary(modannualsSRa)
drop1(modannualsSRa, test = "Chi")

chkres(modannualsSRa,annualsSR$Treatment,annualsSR$Season) # these are not terrible


# nonannuals

plot(SpeciesRichness ~ Treatment, nonannualsSR)
plot(SpeciesRichness ~ Season, nonannualsSR)
plot(SpeciesRichness ~ Order, nonannualsSR)



# try a log-transformed gaussian as that worked previously
nonannualsSR$tSpeciesRichness <- log(nonannualsSR$SpeciesRichness+1)
hist(nonannualsSR$tSpeciesRichness)

modnonannualsSR <- lmer(tSpeciesRichness ~ Treatment*Season + Treatment*Order
                     + (1|Site),
                     data = nonannualsSR)

summary(modnonannualsSR)
drop1(modnonannualsSR, test = "Chi")

chkres(modnonannualsSR,nonannualsSR$Treatment,nonannualsSR$Season) # these are not terrible


# nothing significant, so 
# try without the non-sig interactions

modnonannualsSRa <- lmer(tSpeciesRichness ~ Season + Treatment + Order
                      + (1|Site),
                      data = nonannualsSR)

summary(modnonannualsSRa)
drop1(modnonannualsSRa, test = "Chi")

chkres(modnonannualsSRa,nonannualsSR$Treatment,nonannualsSR$Season) # these are not terrible


# SO - conclusions
# effect of burning on *annuals* is positive - both abundance and species richness
# effect of burning on *everything else* (i.e. perennials, bulbs and shrubs) is negative (slightly) on abundance and neutral on richness
# therefore the flush of floral abundance is entirely driven by annuals

# let's make some...
### figures





### Figures
# we want a figure each from annual abundance, annual SR, non-annual abundance and non-annual SR

# plant abundance
summary(modannualsa)
drop1(modannualsa, test = "Chi")

# create the framework for the plot
newdata17<- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),tPlantCoverage=1)
effects <- data.frame(effect(c("Treatment","Season"),modannualsa))
newdata17 <- merge(newdata17,effects)

newdata17$tfit <- (exp(newdata17$fit))-1
newdata17$tlower <- (exp(newdata17$lower))-1
newdata17$tupper <- (exp(newdata17$upper))-1


newdata17$Treatment <- revalue(newdata17$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g17 <- ggplot(newdata17,
             aes(x=Treatment, y=tfit, fill=Season))+
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  scale_y_continuous(limits=c(0,20),oob=squish)+
  xlab(" ")+ ylab("Abundance")+
  facet_grid(.~Season, switch="x", space = "free_x", scales = "free_x")+
  guides(fill=F)+
  ggtitle("Annuals")+
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
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5))


# visualize the plot
g17



## repeat this for the other three


# annual SR
summary(modannualsSRa)
drop1(modannualsSRa, test="Chi")

# create the framework for the plot
newdata18 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),tSpeciesRichness=1)
effects <- data.frame(effect(c("Treatment","Season"),modannualsSRa))
newdata18 <- merge(newdata18,effects)

newdata18$tfit <- (exp(newdata18$fit))-1
newdata18$tlower <- (exp(newdata18$lower))-1
newdata18$tupper <- (exp(newdata18$upper))-1



newdata18$Treatment <- revalue(newdata18$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g18 <- ggplot(newdata18,
             aes(x=Treatment, y=tfit, fill=Season))+
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  scale_y_continuous(limits = c(0,2), oob=squish)+
  xlab("Season")+ ylab("Species richness")+ 
  facet_grid(.~Season, switch="x", space = "free_x", scales = "free_x")+
  guides(fill=F)+
  theme(legend.text = element_text(size=30),
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        legend.position = "none",
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
g18




# non-annuals abundance
summary(modnonannualsa)
drop1(modnonannualsa, test = "Chi")

# create the framework for the plot
newdata19 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),tPlantCoverage=1)
effects <- data.frame(effect(c("Treatment","Season"),modnonannualsa))
newdata19 <- merge(newdata19,effects)

newdata19$tfit <- (exp(newdata19$fit))-1
newdata19$tlower <- (exp(newdata19$lower))-1
newdata19$tupper <- (exp(newdata19$upper))-1


newdata19$Treatment <- revalue(newdata19$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g19 <- ggplot(newdata19,
             aes(x=Treatment, y=tfit, fill=Season))+
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  scale_y_continuous(limits=c(0,20),oob=squish)+
  xlab(" ")+ ylab(" ")+
  facet_grid(.~Season, switch="x", space = "free_x", scales = "free_x")+
  guides(fill=F)+
  ggtitle("Perennials")+
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
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5))


# visualize the plot
g19





# plant SR
summary(modnonannualsSRa)
drop1(modnonannualsSRa, test="Chi")

# create the framework for the plot
newdata20 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),tSpeciesRichness=1)
effects <- data.frame(effect(c("Treatment","Season"),modnonannualsSRa))
newdata20 <- merge(newdata20,effects)

newdata20$tfit <- (exp(newdata20$fit))-1
newdata20$tlower <- (exp(newdata20$lower))-1
newdata20$tupper <- (exp(newdata20$upper))-1



newdata20$Treatment <- revalue(newdata20$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


# construct the plot
g20 <- ggplot(newdata20,
             aes(x=Treatment, y=tfit, fill=Season))+
  geom_errorbar(aes(ymin = tlower, ymax = tupper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.5))+
  geom_point(cex=3,position=position_dodge(width=0.5),aes(shape = Treatment),
             colour="black", fill="white",stroke=1.2)+
  scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
  scale_y_continuous(limits = c(0,2), oob=squish)+
  xlab("Season")+ ylab(" ")+ 
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
g20



# now stitch them into a multiplot

m.strat <- grid.arrange(g17,g18,g19,g20, ncol=2, nrow=2)

# before we export this we want to move the legend to be outside the plots

# function to extract a legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(g20)

m.strata <- grid.arrange(arrangeGrob(g17,g19,
                                g18,g20 + theme(legend.position="none"),
                                nrow=2, ncol=2),
                    mylegend, ncol =2, widths=c(10,2))






#ggsave("Fig_strategies.svg", plot = m.strata, device = "svg", path = "Results/UpdatedFigs", width = 44, height = 40, units = "cm")
