####   Script for analysis of functionally extrapolated species richness   ####
###############################################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS","plyr","reshape2","vegan","tidyr","RColorBrewer")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

# for some reason the package glmmADMB won't install via the usual methods, so:
#install.packages("R2admb")
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),type="source")
library(glmmADMB)



### load up Callum's custom set of functions
f <- c("CheckResidsFunction.R","CheckConvergenceFunction.R","SpecAccumFunctions.R")
lapply(f, source)


### Moths

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

### tell R that SampleID and SamplingDay should be treated as factors
dframe1$SamplingDay <- factor(dframe1$SamplingDay)
dframe1$Month <- ordered(dframe1$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
dframe1$Year <- factor(dframe1$Year)
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


# We don't need the pollen data or the individual SampleID for this analysis, so get rid of it:

dframe1r <- dframe1[,c(2:4,7:12,79)]
summary(dframe1r)

# each row currently contains 1 insect, so add a column for "Count"; 
# this is 1 in every instance except the one sample with zero insects
dframe1r$Count <- ifelse(dframe1r$Family_Species=="none",0,1)

# summarise the dataframe to how many individuals of each species in each sample
dframe2 <- ddply(dframe1r, .(Family_Species,Site,Treatment,Date,Sample,Season,Month,Year,ExactSeason), numcolwise(sum))

summary(dframe2)

### Functional extrapolation of species richness


# we have two options - either estimate the site-based richness based on occurrance in multiple samples,
# or estimate the sample-based richness based on abundance in the sample

# for the moth data, the latter is probably better due to high seasonal turnover
# however, we could do the former for site_season-based richness as well - ie accumulation for each site across years but within seasons
# doing both has the benefit of letting us see if we get roughly equivalent results

dframe2$SiteSeason <- do.call(paste, c(dframe2[c("Site","Season")], sep = "_"))
dframe2$SiteSeason <- as.factor(dframe2$SiteSeason)

dframe2$ExactNetwork <- do.call(paste, c(dframe2[c("Season","Year","Treatment")], sep = "_"))
dframe2$ExactNetwork <- as.factor(dframe2$ExactNetwork)


# dframe2 contains information of how many of each species are in each sample
# we'll need this information in a matrix format though


matrix1 <- dcast(dframe2, Date + Site + Treatment
                 + Sample + Season + Month + Year + SiteSeason + ExactSeason + ExactNetwork
                 ~ Family_Species,
                 value.var = "Count",
                 fun.aggregate = sum)

# we need to remove the column for "none", which is the 242nd
matrix1 <- matrix1[,c(1:241,243:length(matrix1))]


# first, let's inspect the species accumulation curves for the data
rownames(matrix1) <- matrix1[,4]
matrix1a <- matrix1[,-c(1:10)]

# ...accumulation for each site within each season
matrices1 <- split(matrix1, list(matrix1$SiteSeason))  # this creates a list of smaller dframes, one for each level of sample

# first let's use another of Callum's custom functions to plot all the species accumulation curves
lapply(matrices1, sacplot, cols=10)

# let's also try for each site irrespective of season
matrices2 <- split(matrix1, list(matrix1$Site))  # this creates a list of smaller dframes, one for each level of sample
lapply(matrices2, sacplot, cols=10)



# we can see that very few of these show any signs of nearing an asymptote, even with larger numbers of samples

# therefore we need to extrapolate species richness
# first try sample-based - doing it for every sample based on abundance within the sample
matrices3 <- split(matrix1, list(matrix1$Sample))
SampleSR <- lapply(matrices3, samplebased, cols=10)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
SampleSR.merge <- do.call("cbind", SampleSR)    # merge the data with one sample per column
colnames(SampleSR.merge) <- names(SampleSR)     # assign the sample names to each column
SampleSR.merge <- data.frame(t(SampleSR.merge))
SampleSR.merge


### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the sample-level pollen dataframe:
dframeP<-read.table("Data/SamplesNoct.txt", header=TRUE)
summary(dframeP)


# we need a variable to allow the networks to be paired and ordered:

dframeP$ExactSeason <- do.call(paste, c(dframeP[c("Season","Year")], sep = "_"))
dframeP$ExactSeason <- as.factor(dframeP$ExactSeason)
summary(dframeP$ExactSeason)

dframeP <- merge(SeasonOrder,dframeP)


# we don't need any of the pollen data, so:
dframeP <- dframeP[,1:10]


# for merge to work, we need to set the row names
rownames(dframeP) <- dframeP$Sample

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
SampleSR.full <- merge(dframeP, SampleSR.merge, by=0)





# now try site-based - doing it for every site based on repeated sampling at the site
SiteSR <- lapply(matrices1, sitebased, cols=10)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
SiteSR.merge <- do.call("cbind", SiteSR)    # merge the data with one sample per column
colnames(SiteSR.merge) <- names(SiteSR)     # assign the sample names to each column
SiteSR.merge <- data.frame(t(SiteSR.merge))
SiteSR.merge


### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the sample-level pollen dataframe:
# but it needs collapsing to the same number of sites
dframePS <- matrix1[,c(1:11)]
dframePS <- ddply(dframePS, .(Treatment,Season,SiteSeason), numcolwise(sum))
dframePS <- dframePS[,c(1:3)]

# for merge to work, we need to set the row names
rownames(dframePS) <- dframePS$SiteSeason

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
SiteSR.full <- merge(dframePS, SiteSR.merge, by=0)




# finally try by ExactSeason (18 paired samples)
matrices4 <- split(matrix1, list(matrix1$ExactNetwork))
TreatmentSR <- lapply(matrices4, sitebased, cols=10)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
TreatmentSR.merge <- do.call("cbind", TreatmentSR)    # merge the data with one sample per column
colnames(TreatmentSR.merge) <- names(TreatmentSR)     # assign the sample names to each column
TreatmentSR.merge <- data.frame(t(TreatmentSR.merge))

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the sample-level pollen dataframe:
# but it needs collapsing to the same number of sites
dframePST <- matrix1[,c(1:11)]
dframePST <- ddply(dframePST, .(Treatment,Season,ExactNetwork), numcolwise(sum))
dframePST <- dframePST[,c(1:3)]

# for merge to work, we need to set the row names
rownames(dframePST) <- dframePST$ExactNetwork

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
TreatmentSR.full <- merge(dframePST, TreatmentSR.merge, by=0)


### as an aside, we want to calculate sampling completeness here
# this is just observed species/extrapolated species *100

TreatmentSR.full$completeness <- (TreatmentSR.full$Species*100)/TreatmentSR.full$chao

# we want an average sampling completeness across all networks as well
summary(TreatmentSR.full$completeness)
mean(TreatmentSR.full$completeness, na.rm=T)

# let's also do a weighted interaction completeness - more attention is paid to more species-rich networks?

weighted.mean(TreatmentSR.full$completeness, TreatmentSR.full$chao, na.rm=T)



### now, back to the analysis... let's analyse each of these in turn

# first by sample, using the Chao1 estimated SR
summary(SampleSR.full)

SampleSR.full$Treatment <- relevel(SampleSR.full$Treatment, "NoFire")

hist(SampleSR.full$S.chao1)
hist(log(SampleSR.full$S.chao1+1))

plot(S.chao1 ~ Treatment, SampleSR.full)
plot(S.chao1 ~ Season, SampleSR.full)
plot(S.chao1 ~ interaction(Treatment,Season), SampleSR.full)
plot(S.chao1 ~ Order, SampleSR.full)


model1G <- lmer(log(S.chao1+1) ~ Treatment * Season + Treatment * Order
                + (1|Site),
                data = SampleSR.full)

summary(model1G)
drop1(model1G, test = "Chi")

chkres(model1G,SampleSR.full$Treatment,SampleSR.full$Season)  # these are actually not at all bad

# let's try dropping non-sig interactions
model1Ga <- lmer(log(S.chao1+1) ~ Treatment + Season + Order 
                 + (1|Site),
                 data = SampleSR.full)

summary(model1Ga)
drop1(model1Ga, test = "Chi")

chkres(model1Ga,SampleSR.full$Treatment,SampleSR.full$Season)  # these are actually not at all bad




### Plants

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe3<-read.table("Data/PlantTransects.txt", header=TRUE)

dframe3$Year <- factor(dframe3$Year)
dframe3$Treatment <- relevel(dframe3$Treatment, "NoFire")
summary(dframe3)


# We don't need the individual TransectID for this analysis, so get rid of it:

dframe3r <- dframe3[,c(2:14)]
summary(dframe3r)

# each row currently contains 1 plant species, so add a column for "Count"; 
# this is 1 in every instance except the transect with zero insects
dframe3r$Transects <- ifelse(dframe3r$PlantSpecies=="none",0,1)

# summarise the dataframe to how much total coverage of each species in each sample (non-integers could be a problem so we won't use mean)
dframe4 <- ddply(dframe3r, .(Date,Site,Treatment,Sample,Season,Month,Year,PlantSpecies), numcolwise(sum))

# set up the SiteSeason and TreatmentSeason variables for condensing the data later

dframe4$SiteSeason <- do.call(paste, c(dframe4[c("Site","Season")], sep = "_"))
dframe4$SiteSeason <- as.factor(dframe4$SiteSeason)


dframe4$ExactNetwork <- do.call(paste, c(dframe4[c("Season","Year","Treatment")], sep = "_"))
dframe4$ExactNetwork <- as.factor(dframe4$ExactNetwork)

# put the data out into a matrix format

matrix2 <- dcast(dframe4, Date + Site + Treatment
                 + Sample + Season + Month + Year + SiteSeason + ExactNetwork
                 ~ PlantSpecies,
                 value.var = "PlantCoverage",
                 fun.aggregate = sum)

# we need to remove the column for "none", which is the 52nd and the mystery blank, which is the 10th
matrix2 <- matrix2[,c(1:9,11:51,53:length(matrix2))]

# first, let's inspect the species accumulation curves for the data
rownames(matrix2) <- matrix2[,4]
matrix2a <- matrix2[,-c(1:9)]

# ...accumulation for each site within each season
matricesP1 <- split(matrix2, list(matrix2$SiteSeason))  # this creates a list of smaller dframes, one for each level of sample


# first let's use another of Callum's custom functions to plot all the species accumulation curves
lapply(matricesP1, sacplot,9)

# let's also try for each site irrespective of season
matricesP2 <- split(matrix2, list(matrix2$Site))  # this creates a list of smaller dframes, one for each level of sample
lapply(matricesP2, sacplot,9)


# some of these might be nearing an asymptote but none are there, even with larger numbers of samples

# therefore we need to extrapolate species richness

# first try sample-based - doing it for every sample based on abundance within the sample
matricesP3 <- split(matrix2, list(matrix2$Sample))
SamplePSR <- lapply(matricesP3, samplebased, cols=9)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
SamplePSR.merge <- do.call("cbind", SamplePSR)    # merge the data with one sample per column
colnames(SamplePSR.merge) <- names(SamplePSR)     # assign the sample names to each column
SamplePSR.merge <- data.frame(t(SamplePSR.merge))
SamplePSR.merge


### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the sample-level pollen dataframe:

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
SamplePSR.full <- merge(dframeP, SamplePSR.merge, by=0)


# now try site-based - doing it for every site based on repeated sampling at the site
SitePSR <- lapply(matricesP1, sitebased,9)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
SitePSR.merge <- do.call("cbind", SitePSR)    # merge the data with one sample per column
colnames(SitePSR.merge) <- names(SitePSR)     # assign the sample names to each column
SitePSR.merge <- data.frame(t(SitePSR.merge))

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is previously prepared in dframePS:


# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
SitePSR.full <- merge(dframePS, SitePSR.merge, by=0)



# finally try by Treatment + Season (18 paired samples)
matricesP3 <- split(matrix2, list(matrix2$ExactNetwork))
TreatmentPSR <- lapply(matricesP3, sitebased,9)

# combine all dataframes together
TreatmentPSR.merge <- do.call("cbind", TreatmentPSR)    # merge the data with one sample per column
colnames(TreatmentPSR.merge) <- names(TreatmentPSR)     # assign the sample names to each column
TreatmentPSR.merge <- data.frame(t(TreatmentPSR.merge))


### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is previously prepared in dframePST:

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
TreatmentPSR.full <- merge(dframePST, TreatmentPSR.merge, by=0)



### as an aside, we want to calculate sampling completeness here
# this is just observed species/extrapolated species *100

TreatmentPSR.full$completeness <- (TreatmentPSR.full$Species*100)/TreatmentPSR.full$chao

# we want an average sampling completeness across all networks as well
summary(TreatmentPSR.full$completeness)
mean(TreatmentPSR.full$completeness, na.rm=T)

# let's also do a weighted interaction completeness - more attention is paid to more species-rich networks?

weighted.mean(TreatmentPSR.full$completeness, TreatmentPSR.full$chao, na.rm=T)








### let's analyse each of these in turn, using the estimated SR

# first, samples
summary(SamplePSR.full)

SamplePSR.full$Treatment <- relevel(SamplePSR.full$Treatment, "NoFire")

plot(S.chao1 ~ Treatment, SamplePSR.full)
plot(S.chao1 ~ Season, SamplePSR.full)
plot(S.chao1 ~ interaction(Treatment,Season), SamplePSR.full)
plot(S.chao1 ~ Order, SamplePSR.full)

hist(SamplePSR.full$S.chao1)
hist(log(SamplePSR.full$S.chao1))

model6P <- glmer(S.chao1 ~ Treatment * Season + Treatment * Order
                 + (1|Site),
                 family = poisson (link = "log"),
                 data = SamplePSR.full)

chkconv(model6P)

summary(model6P)
drop1(model6P, test = "Chi")

chkres(model6P, SamplePSR.full$Treatment, SamplePSR.full$Season) # these are fine

# so let's try dropping variables til I find a model that fits

model6Pa <- glmer(S.chao1 ~ Treatment + Season + Order
                 + (1|Site),
                 family = poisson (link = "log"),
                 data = SamplePSR.full)

chkconv(model6Pa)

summary(model6Pa)
drop1(model6Pa, test = "Chi")

chkres(model6Pa, SamplePSR.full$Treatment, SamplePSR.full$Season) # these are fine


####

# in another script I need to plot a multiplot figure using the parameters of model6Pa and model1Ga
# I'm going to predict from them here, then export the predictions in .txt files to read in over there


# moth species richness
summary(model1Ga)
drop1(model1Ga, test="Chi")

# create the framework for the plot
newdata1 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),Order=9,Count=1)
# predict the model outputs
effects <- data.frame(effect(c("Treatment","Season"),model1Ga))
newdata1 <- merge(newdata1,effects)
newdata1$Treatment <- revalue(newdata1$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


newdata1$tfit <- (exp(newdata1$fit))-1
newdata1$tlower <- (exp(newdata1$lower))-1
newdata1$tupper <- (exp(newdata1$upper))-1

newdata1


# floral species richness
summary(model6Pa)
drop1(model6Pa, test="Chi")

# create the framework for the plot
newdata2 <- expand.grid(Treatment=c("Fire","NoFire"),Season=c("Spring","Summer","Autumn","Winter"),Order=9,Count=1)
# predict the model outputs
effects <- data.frame(effect(c("Treatment","Season"),model6Pa))
newdata2 <- merge(newdata2,effects)
newdata2$Treatment <- revalue(newdata2$Treatment, c("Fire"="Burned","NoFire"="Unburned"))
newdata2


# write them both out
write.table(newdata1, file = "Results/mothSRmodel.txt", row.names = F)
write.table(newdata2, file = "Results/plantSRmodel.txt", row.names = F)


# and repeat for order

# moth species richness
summary(model1Ga)
drop1(model1Ga, test="Chi")


# create the framework for the plot
newdata3 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),Count=1)
# predict the model outputs
effects3 <- data.frame(effect(c("Treatment","Order"),model1Ga,xlevels=list(Order=1:9)))
newdata3 <- merge(newdata3,effects3)
newdata3$Treatment <- revalue(newdata3$Treatment, c("Fire"="Burned","NoFire"="Unburned"))


newdata3$tfit <- (exp(newdata3$fit))-1
newdata3$tlower <- (exp(newdata3$lower))-1
newdata3$tupper <- (exp(newdata3$upper))-1

newdata3


# floral species richness
summary(model6Pa)
drop1(model6Pa, test="Chi")

# create the framework for the plot
newdata4 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),Count=1)
# predict the model outputs
effects4 <- data.frame(effect(c("Treatment","Order"),model6Pa,xlevels=list(Order=1:9)))
newdata4 <- merge(newdata4,effects4)
newdata4$Treatment <- revalue(newdata4$Treatment, c("Fire"="Burned","NoFire"="Unburned"))
newdata4


# write them both out
write.table(newdata3, file = "Results/mothSRmodelorder.txt", row.names = F)
write.table(newdata4, file = "Results/plantSRmodelorder.txt", row.names = F)




# create the framework for the plot
newdata5 <- expand.grid(Treatment=c("Fire","NoFire"),Order=c(1:9),Count=1)
# predict the model outputs
effects5 <- data.frame(effect(c("Treatment","Order"),model1Pa,xlevels=list(Order=1:9)))
newdata5 <- merge(newdata5,effects5)
newdata5$Treatment <- revalue(newdata5$Treatment, c("Fire"="Burned","NoFire"="Unburned"))





####


# next, sites

summary(SitePSR.full)

plot(chao ~ Treatment, SitePSR.full)
plot(chao ~ Season, SitePSR.full)
plot(chao ~ interaction(Treatment,Season), SitePSR.full)


hist(SitePSR.full$chao)
hist(log(SitePSR.full$chao))


model4QP <- glm(chao ~ Treatment + Season,
                family = quasipoisson (link = "log"),
                data=SitePSR.full)

summary(model4QP)
drop1(model4QP, test = "Chi")

chkres(model4QP,SitePSR.full$Treatment,SitePSR.full$Season)  # really decent fit here





# second, paired treatment-seasons

summary(TreatmentPSR.full)

# in order for this to be treated as a paired test we need to create a factor with a unique level for each exact season

# first break the ExactSeason variable up at the _ breaks
TreatmentPSR.full <- separate(data = TreatmentPSR.full, col = ExactNetwork, into = c("ExactSeason","Year","Treatment1"), sep = "_")

# then stitch the Season and Year components back together
TreatmentPSR.full$ExactSeason <- do.call(paste, c(TreatmentPSR.full[c("ExactSeason","Year")], sep = "_"))
TreatmentPSR.full$ExactSeason <- as.factor(TreatmentPSR.full$ExactSeason)

# get rid of the columns we don't need
TreatmentPSR.full <- TreatmentPSR.full[,c(2:6,8:length(TreatmentPSR.full))]



plot(chao ~ Treatment, TreatmentPSR.full)
plot(chao ~ Season, TreatmentPSR.full)
plot(chao ~ interaction(Treatment,Season), TreatmentPSR.full)


hist(TreatmentPSR.full$chao)
hist(log(TreatmentPSR.full$chao))

model5QP <- glmmPQL(chao ~ Treatment + Season,
                    random = ~1|ExactNetwork,
                    family = quasipoisson (link = "log"),
                    data=TreatmentPSR.full)

summary(model5QP)
Anova(model5QP, type = "III")

chkres.PQL(model5QP,TreatmentPSR.full$Treatment,TreatmentPSR.full$Season)  # not tooo bad!



model5G <- lmer(log(chao) ~ Treatment * Season
                + (1|ExactNetwork),
                data=TreatmentPSR.full)

summary(model5G)
drop1(model5G, test = "Chi")

chkres(model5G,TreatmentSR.full$Treatment,TreatmentSR.full$Season) # not bad


# try it as a Poisson, so let's try rounding SR estimates to the nearest integer

TreatmentPSR.full$chao.int <- round(TreatmentPSR.full$chao, digits = 0)

plot(chao.int ~ Treatment, TreatmentPSR.full)
plot(chao.int ~ Season, TreatmentPSR.full)
plot(chao.int ~ interaction(Treatment,Season), TreatmentPSR.full)


hist(TreatmentPSR.full$chao.int) # this definitely looks like Poisson - and values are now integer
mean(TreatmentPSR.full$chao.int)
var(TreatmentPSR.full$chao.int)

model5P <- glmer(chao.int ~ Treatment * Season
                 +(1|ExactSeason),
                 family = poisson (link = "log"),
                 data = TreatmentPSR.full)

summary(model5P)
drop1(model5P, test = "Chi")

chkres(model5P,TreatmentPSR.full$Treatment,TreatmentPSR.full$Season) # these aren't terrible but not perfect either









### Community dissimilarity using Adonis (vegan)



### Moths

# we need the data in the original matrix format
summary(matrix1)

rownames(matrix1) <- matrix1$Sample

# but let's trim off some of the excess columns - we only need Treatment, Season and Site
matrix1d <- matrix1[,c(2:3,5,11:length(matrix1))]

# we need to get rid of any rows with zeroes in this time
matrix1d$Total <- rowSums(matrix1d[,c(4:length(matrix1d))])

summary(matrix1d$Total)

matrix1ds <- subset(matrix1d,matrix1d$Total>0)

# we also need a version with all factor columns removed
community1 <- matrix1ds[,c(4:329)]

# first we want to visualise things so as to know what to expect
mod1 <- metaMDS(community1)

MothsTreat <- plot(mod1)
ordihull(mod1, group=matrix1ds$Treatment, show="Fire")
ordihull(mod1, group=matrix1ds$Treatment, show="NoFire")

MothsSeason <- plot(mod1)
ordihull(mod1, group=matrix1ds$Season, show="Spring")
ordihull(mod1, group=matrix1ds$Season, show="Summer")
ordihull(mod1, group=matrix1ds$Season, show="Autumn")
ordihull(mod1, group=matrix1ds$Season, show="Winter")


# based on these plots there may be an effect of Treatment and there almost certainly will be an effect of Season

# now we use adonis to test the dissimilarities between communities in different samples
# we will ask it to test for effects of Treatment and Season (and an interaction), constraining permutations to within Sites
# we'll set the method as Bray-Curtis

AdMoths <- adonis(community1 ~ Treatment*Season,
                  data = matrix1ds,
                  strata = matrix1ds$Site,
                  method = "bray",
                  perm=1e5)

AdMoths

# interaction non-significant but repeating without doesn't really change things

AdMothsa <- adonis(community1 ~ Treatment+Season,
                  data = matrix1ds,
                  strata = matrix1ds$Site,
                  method = "bray",
                  perm=1e5)

AdMothsa

# so may as well use the one with the interaction in



### Plants


# we need the data in the original matrix format
summary(matrix2)


# but let's trim off some of the excess columns - we only need Treatment, Season and Site
matrix2d <- matrix2[,c(2:3,5,10:length(matrix2))]





# we need to get rid of any rows with zeroes in this time
matrix2d$Total <- rowSums(matrix2d[,c(4:length(matrix2d))])

summary(matrix2d$Total)

matrix2ds <- subset(matrix2d,matrix2d$Total>0)





# we also need a version with all factor columns removed
community2 <- matrix2ds[,c(4:74)]
community2[is.na(community2)] <- 0


# first we want to visualise things so as to know what to expect
mod2 <- metaMDS(community2)

FlowersTreat <- plot(mod2)
ordihull(mod2, group=matrix2ds$Treatment, show="Fire")
ordihull(mod2, group=matrix2ds$Treatment, show="NoFire")

FlowersSeason <- plot(mod2)
ordihull(mod2, group=matrix2ds$Season, show="Spring")
ordihull(mod2, group=matrix2ds$Season, show="Summer")
ordihull(mod2, group=matrix2ds$Season, show="Autumn")
ordihull(mod2, group=matrix2ds$Season, show="Winter")

# these plots aren't that informative! Some communities are wildly different from the majority (probably winter!)

# now we use adonis to test the dissimilarities between communities in different samples
# we will ask it to test for effects of Treatment and Season (and an interaction), constraining permutations to within Sites
# we'll set the method as Bray-Curtis

AdFlowers <- adonis(community2 ~ Treatment*Season,
                  data = matrix2ds,
                  strata = matrix2ds$Site,
                  method = "bray",
                  perm=1e5)

AdFlowers



### re-do these at family level
# this is a basic check for any sampling effect in the data that might influence the likelihood of pollen transport between sites

# dframe2 contains information of how many of each species are in each sample
# we'll first need to add in details of which family each species comes from
# this is in a file called AllSpeciesEdited.csv

fams.raw <- read.csv("Data/AllSpeciesEdited.csv", header = T)
summary(fams.raw)

# this file contains the original and the corrected names - the original names share a column name with dframe2, but we need to merge by the corrected names
# therefore let's first discard the original names and rename the columns

fams.raw <- fams.raw[,c(2,4)]

colnames(fams.raw) <- c("Family_Species","Family")

# now we can merge the family data in

fams.merge <- merge(fams.raw,dframe2, all = T)

# discard the things that weren't identified to family
fams.merge <- na.omit(fams.merge)


# summarise the dataframe to how many individuals of each family in each sample
families <- ddply(fams.merge, .(Family,Site,Treatment,Date,Sample,Season,Month,Year,ExactSeason), numcolwise(sum))

summary(families)



# we'll need this information in a matrix format - and we only need Treatment, Season and Site as explanatory variables but Sample too, to maintain the separate rows

matrixF <- dcast(families, Sample + Site + Treatment + Season
                 ~ Family,
                 value.var = "Count",
                 fun.aggregate = sum)


# let's trim off the excess sample column now - we only need Treatment, Season and Site
matrixFd <- matrixF[,c(2:length(matrixF))]

# we need to get rid of any rows with zeroes in this time
matrixFd$Total <- rowSums(matrixFd[,c(4:length(matrixFd))])
summary(matrixFd$Total) # there are none


# we also need a version with all factor columns removed
communityF <- matrixFd[,c(4:34)]

# first we want to visualise things so as to know what to expect
modF <- metaMDS(communityF)

MothsFTreat <- plot(modF)
ordihull(modF, group=matrixFd$Treatment, show="Fire")
ordihull(modF, group=matrixFd$Treatment, show="NoFire")

MothsFSeason <- plot(modF)
ordihull(modF, group=matrixFd$Season, show="Spring")
ordihull(modF, group=matrixFd$Season, show="Summer")
ordihull(modF, group=matrixFd$Season, show="Autumn")
ordihull(modF, group=matrixFd$Season, show="Winter")


# based on these plots there may be an effect of Treatment and there almost certainly will be an effect of Season

# now we use adonis to test the dissimilarities between communities in different samples
# we will ask it to test for effects of Treatment and Season (and an interaction), constraining permutations to within Sites
# we'll set the method as Bray-Curtis

AdMothsF <- adonis(communityF ~ Treatment*Season,
                  data = matrixFd,
                  strata = matrixFd$Site,
                  method = "bray",
                  perm=1e5)

AdMothsF

# interaction non-significant so try repeating without

AdMothsFa <- adonis(communityF ~ Treatment+Season,
                   data = matrixFd,
                   strata = matrixFd$Site,
                   method = "bray",
                   perm=1e5)

AdMothsFa


### do this for plants as well

# dframe4 contains information of how many of each species are in each sample
# we'll first need to add in details of which family each species comes from
# this is in a file called AllPlantSpeciesStrategy.csv

famsP.raw <- read.csv("Data/AllPlantSpecies_strategy.csv", header = T)
summary(famsP.raw)

# this file contains the names and the 'strategy', but let's discard the latter

famsP.raw <- famsP.raw[,c(1,3)]

colnames(famsP.raw) <- c("PlantSpecies","Family")

# now we can merge the family data in

famsP.merge <- merge(famsP.raw,dframe4, all = T)

# discard the things that weren't identified to family, and any transects with 'none'
famsP.merge <- na.omit(famsP.merge)


# summarise the dataframe to how many individuals of each family in each sample
familiesP <- ddply(famsP.merge, .(Family,Site,Treatment,Date,Sample,Season,Month,Year,SiteSeason,ExactNetwork), numcolwise(sum))

summary(familiesP)

# get rid of some excess columns here
familiesP <- familiesP[,c(1:10,14)]



# we'll need this information in a matrix format - and we only need Treatment, Season and Site as explanatory variables but Sample too, to maintain the separate rows

matrixFP <- dcast(familiesP, Sample + Site + Treatment + Season
                 ~ Family,
                 value.var = "PlantCoverage",
                 fun.aggregate = sum)


# let's trim off the excess sample column now - we only need Treatment, Season and Site
matrixFPd <- matrixFP[,c(2:length(matrixFP))]

# we need to get rid of any rows with zeroes in this time
matrixFPd$Total <- rowSums(matrixFPd[,c(4:length(matrixFPd))])
summary(matrixFPd$Total) # there are none


# we also need a version with all factor columns removed
communityFP <- matrixFPd[,c(4:31)]

# first we want to visualise things so as to know what to expect
modFP <- metaMDS(communityFP)

FlowersFTreat <- plot(modFP)
ordihull(modFP, group=matrixFPd$Treatment, show="Fire")
ordihull(modFP, group=matrixFPd$Treatment, show="NoFire")

FlowersFSeason <- plot(modFP)
ordihull(modFP, group=matrixFPd$Season, show="Spring")
ordihull(modFP, group=matrixFPd$Season, show="Summer")
ordihull(modFP, group=matrixFPd$Season, show="Autumn")
ordihull(modFP, group=matrixFPd$Season, show="Winter")


# based on these plots there may be an effect of Treatment and there almost certainly will be an effect of Season

# now we use adonis to test the dissimilarities between communities in different samples
# we will ask it to test for effects of Treatment and Season (and an interaction), constraining permutations to within Sites
# we'll set the method as Bray-Curtis

AdFlowersF <- adonis(communityFP ~ Treatment*Season,
                   data = matrixFPd,
                   strata = matrixFPd$Site,
                   method = "bray",
                   perm=1e5)

AdFlowersF





### check family composition by Treatment + Season (visually)
# just out of interest we want to have a look at the family-level composition in each season
# this is a further collapse of the 'families' dataframe

fams_seasons <- ddply(families, .(Family,Season,Treatment), numcolwise(sum))


matrixFS <- dcast(fams_seasons, Treatment + Season
                 ~ Family,
                 value.var = "Count",
                 fun.aggregate = sum)

# that's not that easy to see what's going on, so let's make everything into proportions of that season's total

matrixFSa <- matrixFS

matrixFS$Total <- rowSums(matrixFS[,3:33])

matrixFSa[,3:33] <- round((matrixFS[,3:33])*100/(matrixFS[,34]), digits=2)


# actually, I want a figure of this. I want to keep every family which contributes >10% of moths in any given row separate, and bundle the rest as "other"
# manually, I have found which columns to keep separate

matrixFSabridged <- matrixFSa[,c(1:2)]

matrixFSabridged$Others <- rowSums(matrixFSa[,-c(1:2,9,13,15,17,21,29)])


matrixFSabridged <- cbind(matrixFSabridged,matrixFSa[,c(29,21,17,15,13,9)])


# now create a column adding all the others

# check this has worked by adding along the rows - this should add to 100 in all cases (but there may be some tiny rounding errors)
test <- rowSums(matrixFSabridged[,-c(1:2)])
summary(test)

# finally we want to melt this into long form
FSabridged <- melt(matrixFSabridged,
                   id=c("Treatment","Season"))

colnames(FSabridged) <- c("Treatment","Season","Family","Percentage")


# wrangle a few variables to control names, order etc.
FSabridged$Treatment <- relevel(FSabridged$Treatment,"Fire")
FSabridged$Treatment <- revalue(FSabridged$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

FSabridged$Season <- factor(FSabridged$Season, levels=c("Spring","Summer","Autumn","Winter"))
FSabridged$Family <- factor(FSabridged$Family, levels=c("Crambidae","Erebidae","Geometridae","Lasiocampidae","Noctuidae","Pyralidae","Others"))

# now produce a plot
g1 <- ggplot(FSabridged, aes(x=Treatment,y=Percentage,fill=Family,group=Season))+
  geom_bar(position = "stack",stat="identity")+
  facet_wrap(~Season)+
  ylab("Percentage composition")+
  scale_fill_brewer(palette = "PRGn")+
  theme(panel.background=element_rect(fill="white"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=15),
        axis.text=element_text(color="black"))


g1




#### repeat for plants

# just out of interest we want to have a look at the family-level composition in each season
# this is a further collapse of the 'families' dataframe

fams_seasonsP <- ddply(familiesP, .(Family,Season,Treatment), numcolwise(sum))


matrixFSP <- dcast(fams_seasonsP, Treatment + Season
                  ~ Family,
                  value.var = "PlantCoverage",
                  fun.aggregate = sum)

# that's not that easy to see what's going on, so let's make everything into proportions of that season's total

matrixFSPa <- matrixFSP

matrixFSP$Total <- rowSums(matrixFSP[,3:30])

matrixFSPa[,3:30] <- round((matrixFSP[,3:30])*100/(matrixFSP[,31]), digits=2)


# actually, I want a figure of this. I want to keep every family which contributes >10% of moths in any given row separate, and bundle the rest as "other"
# manually, I have found which columns to keep separate

matrixFSPabridged <- matrixFSPa[,c(1:2)]

matrixFSPabridged$Others <- rowSums(matrixFSPa[,-c(1:2,7,13,14,16,20,29,30)])


matrixFSPabridged <- cbind(matrixFSPabridged,matrixFSPa[,c(30,29,20,16,14,13,7)])


# now create a column adding all the others

# check this has worked by adding along the rows - this should add to 100 in all cases (but there may be some tiny rounding errors)
test <- rowSums(matrixFSPabridged[,-c(1:2)])
summary(test)

# finally we want to melt this into long form
FSPabridged <- melt(matrixFSPabridged,
                   id=c("Treatment","Season"))

colnames(FSPabridged) <- c("Treatment","Season","Family","Percentage")


# wrangle a few variables to control names, order etc.
FSPabridged$Treatment <- relevel(FSPabridged$Treatment,"Fire")
FSPabridged$Treatment <- revalue(FSPabridged$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

FSPabridged$Season <- factor(FSPabridged$Season, levels=c("Spring","Summer","Autumn","Winter"))
FSPabridged$Family <- factor(FSPabridged$Family, levels=c("Asteraceae","Cistaceae","Ericaceae","Fabaceae","Lamiaceae","Solanaceae","Thymelaeaceae","Others"))

# now produce a plot
g2 <- ggplot(FSPabridged, aes(x=Treatment,y=Percentage,fill=Family,group=Season))+
  geom_bar(position = "stack",stat="identity")+
  facet_wrap(~Season)+
  ylab("Percentage composition")+
  scale_fill_brewer(palette = "PRGn")+
  theme(panel.background=element_rect(fill="white"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=15),
        axis.text=element_text(color="black"))


g2


############################ development ##################################



