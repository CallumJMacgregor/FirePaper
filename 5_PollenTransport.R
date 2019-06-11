########################################################
####   Script for basic pollen transport analysis   ####
########################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS")

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
k <- c("CheckResidsFunction.R","CheckConvergenceFunction.R")
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
dframe1$PollenLoad <- rowSums(dframe1[c(15:length(dframe1))])

### total up the number of pollen species for each insect (i.e. how many columns do not contain 0)
dframe1$PollenTypes <- rowSums(dframe1[c(15:80)] != 0)


### create a binary (yes/no) variable for whether each insect is carrying any pollen
dframe1$PollenYN <- ifelse(dframe1$PollenTypes==0,0,1)
summary(dframe1$PollenYN)

dframe1B <- dframe1[dframe1$Treatment=="Fire",]
dframe1U <- dframe1[dframe1$Treatment=="NoFire",]

summary(dframe1B$PollenYN)
summary(dframe1U$PollenYN)



dframe1spr <- dframe1[dframe1$Season=="Spring",]
dframe1sum <- dframe1[dframe1$Season=="Summer",]
dframe1aut <- dframe1[dframe1$Season=="Autumn",]
dframe1win <- dframe1[dframe1$Season=="Winter",]

summary(dframe1spr$PollenYN)
summary(dframe1sum$PollenYN)
summary(dframe1aut$PollenYN)
summary(dframe1win$PollenYN)


### create a subset dframe containing only the interactions
interactions <- subset(dframe1, select=-c(SampleID,Date,Site,SlideNumber,PollenCount,Treatment,SamplingDay,Sample,PollenTypes,PollenLoad,PollenYN,Season,Month,Year,ExactSeason,Order))
summary(interactions)

### create a subset dataframe containing only pollen-carrying moths
dframe1P <- dframe1[dframe1$PollenYN==1,]


### now you're ready to start looking for patterns!
  
### let's do some exploratory plots first

dframe1Sum <- dframe1[dframe1$Season=="Summer",]
dframe1Spr <- dframe1[dframe1$Season=="Spring",]
dframe1Win <- dframe1[dframe1$Season=="Winter",]
dframe1Aut <- dframe1[dframe1$Season=="Autumn",]

plot(PollenLoad ~ Treatment, dframe1Sum, outline=F)
plot(PollenLoad ~ Treatment, dframe1Spr, outline=F)
plot(PollenLoad ~ Treatment, dframe1Win, outline=F)
plot(PollenLoad ~ Treatment, dframe1Aut, outline=F)

plot(PollenTypes ~ Treatment, dframe1Sum)
plot(PollenTypes ~ Treatment, dframe1Spr)
plot(PollenTypes ~ Treatment, dframe1Win)
plot(PollenTypes ~ Treatment, dframe1Aut)

plot(PollenYN ~ Treatment, dframe1Sum)
plot(PollenYN ~ Treatment, dframe1Spr)
plot(PollenYN ~ Treatment, dframe1Win)
plot(PollenYN ~ Treatment, dframe1Aut)

plot(PollenLoad ~ Order, dframe1)
plot(PollenTypes ~ Order, dframe1)
plot(PollenYN ~ Order, dframe1)


### Let's first look at pollen load per-moth
### Plot it against treatment so you have an idea of what to expect
plot(PollenLoad ~ Treatment, data = dframe1, outline=F)
plot(PollenLoad ~ Season, data = dframe1, outline=FALSE)
plot(PollenLoad ~ Month, data = dframe1, outline=FALSE)
hist(dframe1$PollenLoad)
hist(log(dframe1$PollenLoad+1,10))

### Data clearly have lots of skew, so it's worth checking for overdispersion
### Simplest test is to compare the mean and variance. 
### In a regular Poisson distribution, mean ~= variance; if variance is much larger, data are overdispersed
mean(dframe1$PollenLoad)
var(dframe1$PollenLoad)

### Data are overdispersed, so we will try two types of model in addition to Poisson: quasi-Poisson, and negative binomial
### Zero-inflated models might be appropriate to try as well

# construct models using Date and Site as random effects
# random effects are factors that might affect the output variable, but not in a way that is interesting to us
# you might see variation between different sampling days, or between different fire or non-fire sites...
# ...but we are only really interested in variation due to Treatment


# Poisson model using lme4

model1P <- glmer(PollenLoad ~ Treatment*Season + Treatment*Order # fixed effects
                 + (1|Site), # random effects
                 family = poisson (link = "log"),
                 data = dframe1)

chkconv(model1P)

# inspect and test the model
summary(model1P)

# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model1P, test="Chisq")  



# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you want residuals that are roughly normally distributed around zero with no obvious trends
chkres(model1P, dframe1$Treatment, dframe1$Season)  # these residuals do appear to have a negative trend so they are not ideal; this might be driven by zeroes


# QuasiPoisson model using MASS

model1Q <- glmmPQL(PollenLoad ~ Treatment*Season + Treatment*Order,
                   random = list(~1|Site),
                   family = quasipoisson (link = "log"),
                   data = dframe1)


summary(model1Q)
Anova(model1Q, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
# this function produces a subset of the previous plots that are available for this model class
chkres.PQL(model1Q, dframe1$Treatment, dframe1$Season) # these are bad


# Zero-inflated Poisson

#model1ZIP <- glmmadmb(PollenLoad ~  Treatment*Season + Treatment*Order
#                      + (1|Site), #Random effects
#                      zeroInflation=TRUE,
#                      family = "poisson",
#                      data = dframe1)

#summary(model1ZIP)

#Anova(model1ZIP, type="III")
#drop1(model1ZIP, test="Chisq")

# chkres.zi(model1ZIP, dframe1$Treatment, dframe1$Season) # these aren't great either

# negative binomial

#model1NB <- glmer.nb(PollenLoad ~ Treatment*Season + Treatment*Order # fixed effects
#                     + (1|Site), # random effects
#                     data = dframe1)

#chkconv(model1NB)

#summary(model1NB)
#drop1(model1NB, test= "Chisq") # glmer.nb produces the same model class as glmer so we can treat it the same

#chkres(model1NB, dframe1$Treatment, dframe1$Season) # these are still worse than just the Poisson

# Gaussian with log transformation

model1G <- lmer(log(PollenLoad+1) ~ Treatment*Season + Treatment*Order
                + (1|Site),
                data = dframe1)

summary(model1G)
drop1(model1G, test = "Chi")

chkres(model1G, dframe1$Treatment, dframe1$Season)   # these are still affected by the zeroes but are probably the most balanced yet

### choose from these candidate error families
### we can see from the residuals that model1G is the best option so we use this result

summary(model1G)
drop1(model1G, test = "Chi")




### retry this using pollen carriers only

plot(PollenLoad ~ interaction(Treatment,Season), data = dframe1P)
plot(PollenTypes ~ interaction(Treatment,Season), data = dframe1P)

hist(dframe1P$PollenLoad)
hist(log(dframe1P$PollenLoad,10))

# Poisson

model1PP <- glmer(PollenLoad ~ Treatment*Season + Treatment*Order
                  + (1|Site),
                  family=poisson (link="log"),
                  data = dframe1P)


summary(model1PP)
drop1(model1PP, test = "Chi")

chkres(model1PP, dframe1P$Treatment, dframe1P$Season)  # these are bad

# Gaussian + log transformation

model1PG <- lmer(log(PollenLoad) ~ Treatment*Season + Treatment*Order
                 + (1|Site),
                 data = dframe1P)

summary(model1PG)
drop1(model1PG, test = "Chi")


chkres(model1PG, dframe1P$Treatment, dframe1P$Season)   # these still have a trend from having many 1s but are way better than anything else we've seen

# remove non-significant interaction

model1PGa <- lmer(log(PollenLoad) ~ Treatment*Season + Order
                  + (1|Site),
                  data = dframe1P)

summary(model1PGa)
drop1(model1PGa, test = "Chi")

chkres(model1PGa, dframe1P$Treatment, dframe1P$Season)


### Let's now look at pollen types per-moth
### Plot it against treatment so you have an idea of what to expect
# start with the pollen-carriers only since that worked well last time
plot(PollenTypes ~ Treatment, data = dframe1P)
plot(PollenTypes ~ Season, data = dframe1P)
hist(dframe1P$PollenTypes)
hist(log(dframe1P$PollenTypes))

### Again data clearly have lots of skew, though not so bad, so it's worth checking for overdispersion
### Simplest test is to compare the mean and variance. 
### In a regular Poisson distribution, mean ~= variance; if variance is much larger, data are overdispersed
mean(dframe1P$PollenTypes)
var(dframe1P$PollenTypes)

### Data appear possibly overdispersed, so we will be careful

# construct models using Date and Site as random effects


# Poisson model using lme4

model2P <- glmer(PollenTypes ~ Treatment * Season + Treatment * Order # fixed effects
                 + (1|Site), # random effects
                 family = poisson (link = "log"),
                 data = dframe1P)

chkconv(model2P)

# inspect and test the model
summary(model2P)

# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model2P, test="Chisq")  



# check the model's residuals
chkres(model2P, dframe1P$Treatment, dframe1P$Season)  # these residuals do appear to have a slight positive trend but probably nothing to worry about too much



# Gaussian with log transformation

model2G <- lmer(log(PollenTypes) ~ Treatment * Season + Treatment * Order
                + (1|Site),
                data = dframe1P)

summary(model2G)
drop1(model2G, test = "Chi")

chkres(model2G, dframe1P$Treatment, dframe1P$Season)  # these don't improve sufficiently on the Poisson to justify using this inferior model type



### choose from these candidate error families
summary(model2P)
drop1(model2P, test = "Chi")


# remove a non-sig interaction at a time - neither is significant so remove both...

model2Pi <- glmer(PollenTypes ~ Treatment + Season + Order # fixed effects
                 + (1|Site), # random effects
                 family = poisson (link = "log"),
                 data = dframe1P)

chkres(model2Pi, dframe1P$Treatment, dframe1P$Season)

# inspect and test the model
summary(model2Pi)

# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model2Pi, test="Chisq")  




### finally, let's look at proportion of moths carrying pollen

### Plot it against treatment so you have an idea of what to expect
plot(PollenYN ~ Treatment, data = dframe1)
plot(PollenYN ~ Season, data = dframe1)
plot(PollenYN ~ interaction(Treatment,Season), data = dframe1)
plot(PollenYN ~ Order, data = dframe1)
hist(dframe1$PollenYN)

### this data is definitely binomial, so we don't need to worry too much about model selection or residuals:

model3B <- glmer(PollenYN~Treatment*Season + Treatment*Order
                 + (1|Site),
                 family = binomial (link = "logit"),
                 data = dframe1)

summary(model3B)
drop1(model3B, test="Chi")

# interaction is non-sig so let's try subsetting out

model3Ba <- glmer(PollenYN~Treatment*Season + Order
                  + (1|Site),
                  family = binomial (link = "logit"),
                  data = dframe1)

chkconv(model3Ba)

summary(model3Ba)
drop1(model3Ba, test="Chi")






