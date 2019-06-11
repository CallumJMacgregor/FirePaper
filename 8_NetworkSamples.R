###########################################################
####   Script for pollen transport analysis by sample  ####
###########################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS","glmmADMB")

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
source("CheckResidsFunction.R")



### Flower visitation ###


### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/NetworkMetricsFVbySample.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly


### now you're ready to start looking for patterns!

### let's start simple - links per species

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$links.per.species)                   # this looks like a Poisson distribution, but...
hist(log(dframe1$links.per.species,10))           # ... the logarithm looks like a normal distribution
plot(links.per.species ~ Treatment, data = dframe1)


# construct models using Date and Site as random effects
# random effects are factors that might affect the output variable, but not in a way that is interesting to us
# you might see variation between different sampling days, or between different fire or non-fire sites...
# ...but we are only really interested in variation due to Treatment

# although the data look Poisson, technically Poisson is only meant to handle integer values, which these are not...
# instead, let's use Gaussian distribution on the logarithm of the values

model1LPS <- lmer(log(links.per.species,10) ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 data = dframe1)

# inspect and test the model
summary(model1LPS)


# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model1LPS, test="Chi")  


# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you want residuals that are roughly normally distributed around zero (in the histogram)
# and with no obvious trends (in the scatterplots)
chkres(model1LPS)  # these residuals look fine




### next - linkage density

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$linkage.density)                   # this looks like a Poisson distribution, but...
hist(log(dframe1$linkage.density,10))           # ... the logarithm looks like a normal distribution
plot(linkage.density ~ Treatment, data = dframe1)


# construct models using Date and Site as random effects

model1LD <- lmer(log(linkage.density,10) ~ Treatment # fixed effects
                  + (1|Date) + (1|Site), # random effects
                  data = dframe1)

# inspect and test the model
summary(model1LD)
drop1(model1LD, test="Chi")  


# check the model's residuals
chkres(model1LD)  # these residuals look fine




### next - selectivity

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$H2)                   # this looks like a Poisson distribution, but...
hist(log(dframe1$H2,10))           # ... the logarithm looks like a normal distribution
plot(H2 ~ Treatment, data = dframe1)


# construct models using Date and Site as random effects
# here the data include zeroes. Log(0) doesn't exist, so we can't use our previous approach.
# nor can we use either Poisson or negative binomial, as the data is non-integer
# therefore, use a QuasiPoisson model using MASS

model1H2 <- glmmPQL(H2 ~ Treatment,
                   random = list(~1|Date, ~1|Site),
                   family = quasipoisson (link = "log"),
                   data = dframe1)


summary(model1H2)
Anova(model1H2, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
# this function produces a subset of the previous plots that are available for this model class
chkres.PQL(model1H2) # these look ok



### next - interaction evenness

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$interaction.evenness)                   # this looks like a normal distribution
plot(interaction.evenness ~ Treatment, data = dframe1)


# construct models using Date and Site as random effects

model1IE <- lmer(interaction.evenness ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 data = dframe1)

# inspect and test the model
summary(model1IE)
drop1(model1IE, test="Chi")  


# check the model's residuals
chkres(model1IE)  # these residuals look ok




### next - weighted nestedness (NODF)

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$weighted.nestedness)                   # this looks roughly like a normal distribution
plot(weighted.nestedness ~ Treatment, data = dframe1)


# construct models using Date and Site as random effects

model1NODF <- lmer(weighted.nestedness ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 data = dframe1)

# inspect and test the model
summary(model1NODF)
drop1(model1NODF, test="Chi")  


# check the model's residuals
chkres(model1NODF)  # these residuals look mostly ok



### next - robustness

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$robustness)                   # this looks roughly like a Poisson distribution but non-integer
hist(log(dframe1$robustness,10))           # this looks roughly like a normal distribution
plot(robustness ~ Treatment, data = dframe1)


# construct models using Date and Site as random effects

model1RB <- lmer(log(robustness,10) ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 data = dframe1)

# inspect and test the model
summary(model1RB)
drop1(model1RB, test="Chi")  


# check the model's residuals
chkres(model1RB)  # these residuals aren't too bad - a hint of a positive trend


# try different error families

model1RBa <- glmmPQL(robustness ~ Treatment,
                    random = list(~1|Date, ~1|Site),
                    family = quasipoisson (link = "log"),
                    data = dframe1)

summary(model1RBa)
Anova(model1RBa, type="III")

chkres.PQL(model1RBa)


# Try a Gamma distribution?


model1RBb <- glmer(robustness ~ Treatment # fixed effects
                  + (1|Date) + (1|Site), # random effects
                  family = Gamma(link="log"),
                  data = dframe1)

# inspect and test the model
summary(model1RBb)
drop1(model1RBb, test="Chi")  


# check the model's residuals
chkres(model1RBb)  # these residuals are definitely worse



# inverse Gaussian (we don't often use this but conceptually it could work probably right for robustness)


model1RBc <- glmer(robustness ~ Treatment # fixed effects
                  + (1|Date) + (1|Site), # random effects
                  family = inverse.gaussian(link="log"),
                  data = dframe1)

# inspect and test the model
summary(model1RBc)
drop1(model1RBc, test="Chi")  


# check the model's residuals
chkres(model1RBc)  # these residuals are surprisingly pretty good!






### Pollen transport ###



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe2<-read.table("Data/NetworkMetricsPTbySample.txt", header=TRUE)

summary(dframe2) # Check it's imported correctly


### now you're ready to start looking for patterns!

### let's start simple - links per species

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$links.per.species)                   # this looks like a Poisson distribution, but...
hist(log(dframe2$links.per.species,10))           # ... the logarithm looks like a normal distribution
plot(log(links.per.species,10) ~ Treatment, data = dframe2)


# construct models using Date and Site as random effects
# random effects are factors that might affect the output variable, but not in a way that is interesting to us
# you might see variation between different sampling days, or between different fire or non-fire sites...
# ...but we are only really interested in variation due to Treatment

# although the data look Poisson, technically Poisson is only meant to handle integer values, which these are not...
# instead, let's use Gaussian distribution on the logarithm of the values

model2LPS <- lmer(log(links.per.species,10) ~ Treatment # fixed effects
                  + (1|Date) + (1|Site), # random effects
                  data = dframe2)

# inspect and test the model
summary(model2LPS)


# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model2LPS, test="Chi")  


# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you want residuals that are roughly normally distributed around zero (in the histogram)
# and with no obvious trends (in the scatterplots)
chkres(model2LPS)  # these residuals look fine




### next - linkage density

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$linkage.density)                   # this looks like a Poisson distribution, but...
hist(log(dframe2$linkage.density,10))           # ... the logarithm looks like a normal distribution
plot(log(linkage.density,10) ~ Treatment, data = dframe2)


# construct models using Date and Site as random effects

model2LD <- lmer(log(linkage.density,10) ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 data = dframe2)

# inspect and test the model
summary(model2LD)
drop1(model2LD, test="Chi")  


# check the model's residuals
chkres(model2LD)  # these residuals look not great but not too bad




### next - selectivity

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$H2)                   # this looks like a normal distribution
plot(H2 ~ Treatment, data = dframe2)


# construct models using Date and Site as random effects
# here the data include zeroes. Log(0) doesn't exist, so we can't use our previous approach.
# nor can we use either Poisson or negative binomial, as the data is non-integer
# therefore, use a QuasiPoisson model using MASS

model2H2 <- lmer(H2 ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 data = dframe2)


summary(model2H2)
drop1(model2H2, test = "Chi")

# check residuals
chkres(model2H2) # these have a really clear trend so we need to try other error distributions

# metric bounded at 0 and 1 so let's try inverse-Gaussian first
# we need rid of the zeroes for this

dframe2$H2 <- 0.00001+dframe2$H2
summary(dframe2$H2)

model2H2a <- glmer(H2 ~ Treatment # fixed effects
                   + (1|Date) + (1|Site), # random effects
                   family = inverse.gaussian(link = "log"),
                   data = dframe2)

summary(model2H2a)
drop1(model2H2a, test = "Chi")

chkres(model2H2a)  # much better but still not perfect


# Gamma?

model2H2b <- glmer(H2 ~ Treatment # fixed effects
                   + (1|Date) + (1|Site), # random effects
                   family = Gamma(link = "identity"),
                   data = dframe2)

summary(model2H2b)
drop1(model2H2b, test = "Chi")

chkres(model2H2b)  # definite trend - worse than inverse-Gaussian




### next - interaction evenness

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$interaction.evenness)                   # this looks like a normal distribution
plot(interaction.evenness ~ Treatment, data = dframe2)


# construct models using Date and Site as random effects

model2IE <- lmer(interaction.evenness ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 data = dframe2)

# inspect and test the model
summary(model2IE)
drop1(model2IE, test="Chi")  


# check the model's residuals
chkres(model2IE)  # these residuals look ok




### next - weighted nestedness (NODF)

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$weighted.nestedness)                   # this looks roughly like a normal distribution
plot(weighted.nestedness ~ Treatment, data = dframe2)


# construct models using Date and Site as random effects

model2NODF <- lmer(weighted.nestedness ~ Treatment # fixed effects
                   + (1|Date) + (1|Site), # random effects
                   data = dframe2)

# inspect and test the model
summary(model2NODF)
drop1(model2NODF, test="Chi")  


# check the model's residuals
chkres(model2NODF)  # these residuals look mostly ok



### next - robustness

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$robustness)                   # this looks roughly like a Poisson distribution but non-integer
hist(log(dframe2$robustness,10))           # this looks roughly like a normal distribution
plot(robustness ~ Treatment, data = dframe2)


# construct models using Date and Site as random effects

model2RB <- lmer(log(robustness,10) ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 data = dframe2)

# inspect and test the model
summary(model2RB)
drop1(model2RB, test="Chi")  


# check the model's residuals
chkres(model2RB)  # these residuals aren't too bad - a hint of a positive trend


# try different error families

model2RB <- glmmPQL(robustness ~ Treatment,
                    random = list(~1|Date, ~1|Site),
                    family = quasipoisson (link = "log"),
                    data = dframe1)

# Quasi-poisson model fails to converge! Eek.


# Try a Gamma distribution?


model3RB <- glmer(robustness ~ Treatment # fixed effects
                  + (1|Date) + (1|Site), # random effects
                  family = Gamma(link="log"),
                  data = dframe1)

# inspect and test the model
summary(model3RB)
drop1(model3RB, test="Chi")  


# check the model's residuals
chkres(model3RB)  # these residuals are definitely worse



# inverse Gaussian (we don't often use this but conceptually it could work)


model4RB <- glmer(robustness ~ Treatment # fixed effects
                  + (1|Date) + (1|Site), # random effects
                  family = inverse.gaussian(link="log"),
                  data = dframe1)

# inspect and test the model
summary(model4RB)
drop1(model4RB, test="Chi")  


# check the model's residuals
chkres(model4RB)  # these residuals are surprisingly pretty good!





