#############################################################################
####   Script for plotting figures from the sampling completeness data   ####
#############################################################################


### We need the end summary tables from scripts 3 & 4, but not much else from either script, so we can trim a lot out

### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS","plyr","reshape2","vegan","tidyr","data.table","svglite")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

# for some reason the package glmmADMB won't install via the usual methods, so:
#install.packages("R2admb")
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),type="source")
library(glmmADMB)



### load up Callum's custom set of functions
f <- c("CheckResidsFunction.R","CheckConvergenceFunction.R","SpecAccumFunctions.R","NetworkFunction.R")
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


# We don't need the pollen data or the individual SampleID for this analysis, so get rid of it:

dframe1r <- dframe1[,c(2:4,7:12)]
summary(dframe1r)

# each row currently contains 1 insect, so add a column for "Count"; 
# this is 1 in every instance except the one sample with zero insects
dframe1r$Count <- ifelse(dframe1r$Family_Species=="none",0,1)

# summarise the dataframe to how many individuals of each species in each sample
dframe2 <- ddply(dframe1r, .(Family_Species,Site,Treatment,Date,Sample,Season,Month,Year), numcolwise(sum))



### Functional extrapolation of species richness


# we have two options - either estimate the site-based richness based on occurrance in multiple samples,
# or estimate the sample-based richness based on abundance in the sample

# for the moth data, the latter is probably better due to high seasonal turnover
# however, we could do the former for site_season-based richness as well - ie accumulation for each site across years but within seasons
# doing both has the benefit of letting us see if we get roughly equivalent results

dframe2$SiteSeason <- do.call(paste, c(dframe2[c("Site","Season")], sep = "_"))
dframe2$SiteSeason <- as.factor(dframe2$SiteSeason)

dframe2$ExactSeason <- do.call(paste, c(dframe2[c("Season","Year","Treatment")], sep = "_"))
dframe2$ExactSeason <- as.factor(dframe2$ExactSeason)


# dframe2 contains information of how many of each species are in each sample
# we'll need this information in a matrix format though


matrix1 <- dcast(dframe2, Date + Site + Treatment
                 + Sample + Season + Month + Year + SiteSeason + ExactSeason
                 ~ Family_Species,
                 value.var = "Count",
                 fun.aggregate = sum)

# we need to remove the column for "none", which is the 242nd
matrix1 <- matrix1[,c(1:241,243:length(matrix1))]


# first, let's inspect the species accumulation curves for the data
rownames(matrix1) <- matrix1[,4]
matrix1a <- matrix1[,-c(1:9)]

# ...accumulation for each site within each season
matrices1 <- split(matrix1, list(matrix1$SiteSeason))  # this creates a list of smaller dframes, one for each level of sample

# first let's use another of Callum's custom functions to plot all the species accumulation curves
lapply(matrices1, sacplot, cols=9)

# let's also try for each site irrespective of season
matrices2 <- split(matrix1, list(matrix1$Site))  # this creates a list of smaller dframes, one for each level of sample
lapply(matrices2, sacplot, cols=9)



# we can see that very few of these show any signs of nearing an asymptote, even with larger numbers of samples

# therefore we need to extrapolate species richness
# first try sample-based - doing it for every sample based on abundance within the sample
matrices3 <- split(matrix1, list(matrix1$Sample))
SampleSR <- lapply(matrices3, samplebased, cols=9)

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

# we don't need any of the pollen data, so:
dframeP <- dframeP[,1:8]

# but we do need a Treatment variable, so:
dframeP$Treatment <- ifelse(grepl("NF",dframeP$Site),"NoFire","Fire")
dframeP$Treatment <- factor(dframeP$Treatment)

# for merge to work, we need to set the row names
rownames(dframeP) <- dframeP$Sample

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
SampleSR.full <- merge(dframeP, SampleSR.merge, by=0)





# now try site-based - doing it for every site based on repeated sampling at the site
SiteSR <- lapply(matrices1, sitebased, cols=9)

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
dframePS <- matrix1[,c(1:10)]
dframePS <- ddply(dframePS, .(Treatment,Season,SiteSeason), numcolwise(sum))
dframePS <- dframePS[,c(1:3)]

# for merge to work, we need to set the row names
rownames(dframePS) <- dframePS$SiteSeason

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
SiteSR.full <- merge(dframePS, SiteSR.merge, by=0)




# finally try by ExactSeason (18 paired samples)
matrices4 <- split(matrix1, list(matrix1$ExactSeason))
TreatmentSR <- lapply(matrices4, sitebased, cols=9)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
TreatmentSR.merge <- do.call("cbind", TreatmentSR)    # merge the data with one sample per column
colnames(TreatmentSR.merge) <- names(TreatmentSR)     # assign the sample names to each column
TreatmentSR.merge <- data.frame(t(TreatmentSR.merge))

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the sample-level pollen dataframe:
# but it needs collapsing to the same number of sites
dframePST <- matrix1[,c(1:10)]
dframePST <- ddply(dframePST, .(Treatment,Season,ExactSeason), numcolwise(sum))
dframePST <- dframePST[,c(1:3)]

# for merge to work, we need to set the row names
rownames(dframePST) <- dframePST$ExactSeason

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


dframe4$ExactSeason <- do.call(paste, c(dframe4[c("Season","Year","Treatment")], sep = "_"))
dframe4$ExactSeason <- as.factor(dframe4$ExactSeason)

# put the data out into a matrix format

matrix2 <- dcast(dframe4, Date + Site + Treatment
                 + Sample + Season + Month + Year + SiteSeason + ExactSeason
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

# sample-based won't work here as it's dependent on the number of singletons
# as this is percent cover, 1% cover would be a singleton, which isn't really appropriate


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
matricesP3 <- split(matrix2, list(matrix2$ExactSeason))
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






################################################################
####   Script for network analysis by season and treatment  ####
################################################################



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe5<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe5) # Check it's imported correctly

dframe5$ExactSeason <- do.call(paste, c(dframe5[c("Season","Year","Treatment")], sep = "_"))
dframe5$ExactSeason <- as.factor(dframe5$ExactSeason)


### prepare the data for network analysis

# trim off extra columns and reorder
dframe5r <- dframe5[,c(2,7,10,length(dframe5),12:(length(dframe5)-1))]


# break up dataframe per species and per Season-Year-Treatment
# now we have a massive number of dframes, each containing every individual sampled of a given species in a given season/treatment combo

dframes <- split(dframe5r, list(dframe5r$ExactSeason,dframe5r$Family_Species))  # this creates a list of smaller dframes, one for each level of sample
summary(dframes)


# now apply another custom function to extract interaction richness and extrapolated Chao richness for every species-season-treatment combo
interaction.completeness <- lapply(dframes, interaction.complete, cols=5, threshold=0)

# this function only gets applied where 10 or more individuals of a species are sampled, so there will be a whole bunch of warnings where this isn't the case
warnings()

# now we have a list of every species-season-treatment combo, with interaction richness estimates where possible and error messages everywhere else
# let's turn it into a dataframe
ICmerge <- do.call("cbind", interaction.completeness)    # merge the data with one sample per column
colnames(ICmerge) <- names(interaction.completeness)     # assign the sample names to each column
ICmerge <- data.frame(t(ICmerge))

# and get rid of all the ones with too few individuals to reliably estimate
ICgood <- ICmerge[ which( ! ICmerge$Species %in% "Fewer than 10 individuals sampled") , ]

# now we want to calculate the % interaction completeness
# first we need to make sure the outputs from the loop are stored as numbers rather than factors
ICgood$nSpecies <- as.numeric(as.character(ICgood$Species))
ICgood$nchao <- as.numeric(as.character(ICgood$chao))
ICgood$nn <- as.numeric(as.character(ICgood$n))


# check that's worked correctly - if so this should return 100% TRUE (some NAs are acceptable for Chao only, but no FALSE)
ICgood$testSpec <- ifelse(ICgood$Species==ICgood$nSpecies,T,F)
ICgood$testChao <- ifelse(ICgood$chao==ICgood$nchao,T,F)
ICgood$testn <- ifelse(ICgood$n==ICgood$nn,T,F)

summary(ICgood$testSpec)  
summary(ICgood$testChao)
summary(ICgood$testn)

# next we want to get rid of those NAs (these are instances where no pollen has been sampled despite 5+ individuals being examined - IR=0 so Chao fails)
ICgood$nchao <- ifelse(ICgood$nSpecies==0,0,ICgood$nchao)

# now we can calculate the % completeness:
ICgood$completeness <- (ICgood$nSpecies*100)/ICgood$nchao

# wherever no pollen was sampled, it's produced NaNs in the completeness column
# we need to remove these in order to calculate means
ICgood <- ICgood[ which( ! ICgood$completeness %in% "NaN") , ]


summary(ICgood$completeness)
mean(ICgood$completeness, na.rm=T)

# let's also do a weighted interaction completeness - more attention is paid to more generalist insects?

weighted.mean(ICgood$completeness, ICgood$nchao, na.rm=T)

## LATER EDIT - let's re-do the full network completeness following Macgregor et al. 2017 bioRxiv.
# the necessary code is in Appendix S1.2 of that paper - copied here

source("AppendixS1.2.R")
SCw2(dframe5r, cols = 5, species.col = "Family_Species")

# remember, though, that full network completeness here is not consistent with the estimates for moths and plants -
# we want to calculate completeness for each of the 18 networks and then take a mean

##

# now we want to do this for each of the 18 networks individually - so for that we need to extract the network titles
# this is actually easier to do bit-by-bit than all at once

# season
ICgood$Season <- ifelse(grepl("Summer",rownames(ICgood)),"Summer",
                        ifelse(grepl("Spring",rownames(ICgood)),"Spring",
                               ifelse(grepl("Winter",rownames(ICgood)),"Winter",
                                      ifelse(grepl("Autumn",rownames(ICgood)),"Autumn","Fail"))))

ICgood$Season <- factor(ICgood$Season)
summary(ICgood$Season)

# year
ICgood$Year <- ifelse(grepl("1",rownames(ICgood)),1,
                      ifelse(grepl("2",rownames(ICgood)),2,
                             ifelse(grepl("3",rownames(ICgood)),3,"Fail")))

ICgood$Year <- factor(ICgood$Year)
summary(ICgood$Year)

# treatment
ICgood$Treatment <- ifelse(grepl("NoFire",rownames(ICgood)),"NoFire","Fire")

ICgood$Treatment <- factor(ICgood$Treatment)
summary(ICgood$Treatment)


# stitch them back together
ICgood$ExactSeason <- do.call(paste, c(ICgood[c("Season","Year","Treatment")], sep = "_"))
ICgood$ExactSeason <- as.factor(ICgood$ExactSeason)
summary(ICgood$ExactSeason)



# now, calculate means for each level of ExactSeason
means <- aggregate(ICgood[,'completeness'],list(ICgood$ExactSeason),mean, na.rm=T)

colnames(means) <- c("ExactSeason","Mean")

# do the same for weighted mean
DT <- data.table(ICgood)
weighted.means <- DT[,list(wret = weighted.mean(completeness,nchao, na.rm=T)),by=ExactSeason]

colnames(weighted.means) <- c("ExactSeason","WeightedMean")


# we also want information about how many species are used to calculate these values for each network:
ICgood$n.species <- 1

n.species <- aggregate(ICgood[,'n.species'],list(ICgood$ExactSeason),sum)
colnames(n.species) <- c("ExactSeason","n.species")  


# merge these three together to give final values
means.both <- merge(means,weighted.means)
means.full <- merge(n.species,means.both)

means.full$Treatment <- factor(ifelse(grepl("NoFire",means.full$ExactSeason),"NoFire","Fire"))

# now we can take the full network mean in a way that's consistent with that for moths and plants

mean(means.full$WeightedMean, na.rm=T)



# we now have all the necessary information for plotting the figure, held in two tables with a bunch of extra stuff we don't need

summary(means.full)
summary(TreatmentSR.full)
summary(TreatmentPSR.full)


insects <- TreatmentSR.full[,c(2,4,14)]
plants <- TreatmentPSR.full[,c(2,4,14)]
interactions <- means.full[,c(5,1,4)]


# we need to change some of the column headers so the merge will happen properly, and add a label column to each table
colnames(insects) <- c("Treatment","ExactSeason","completeness")
colnames(plants) <- c("Treatment","ExactSeason","completeness")
colnames(interactions) <- c("Treatment","ExactSeason","completeness")

insects$Level <- factor("Insects")
plants$Level <- factor("Plants")
interactions$Level <- factor("Interactions")

# merge them together

all <- rbind(insects,plants,interactions)
summary(all)

mean(insects$completeness)
mean(plants$completeness)
mean(interactions$completeness)

# t tests for differences between treatments
t.test(completeness ~ Treatment, data=insects)
t.test(completeness ~ Treatment, data=plants)
t.test(completeness ~ Treatment, data=interactions)


# as well as the raw points we need a table with the mean and se and conf limits for each level+treatment combo
# we could do it from model predictions

model1 <- glm(completeness ~ Treatment * Level,     # including interaction ensures model pred = mean
              data = all)

summary(model1)
drop1(model1, test="F")

model1a <- glm(completeness ~ Treatment + Level,
               data = all)

summary(model1a)
drop1(model1a, test="F")

newdata <- expand.grid(Treatment=c("Fire","NoFire"),Level=c("Insects","Plants","Interactions"),completeness=1)
mm <- model.matrix(terms(model1),newdata)  

newdata$completeness <- predict(model1,newdata=newdata,type="response")

predict <- data.frame(predict(model1,newdata=newdata,type="response",se.fit=T))

newdata <- cbind(newdata,predict)  

newdata$ucl <- newdata$completeness+1.96*newdata$se.fit
newdata$lcl <- newdata$completeness-1.96*newdata$se.fit

summary(newdata)

all$Treatment <- relevel(all$Treatment,"Fire")
newdata$Treatment <- relevel(newdata$Treatment,"Fire")

newdata$Treatment <- revalue(newdata$Treatment, c("Fire"="Burned","NoFire"="Unburned"))
all$Treatment <- revalue(all$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

newdata$Level <- revalue(newdata$Level, c("Insects"="Moths"))
all$Level <- revalue(all$Level, c("Insects"="Moths"))



### figures

g1 <- ggplot(data = all, aes(x = Treatment, y = completeness, group = Level)) +
              geom_point(cex=0.5,position=position_dodge(width=0.5),color="gray45")+
              ylim(0,100)+
              geom_errorbar(data = newdata, aes(x = Treatment, y = fit, ymin = lcl,ymax = ucl),position=position_dodge(width=0.5),
                colour="black",size=0.6, width=0.5)+
              geom_point(data = newdata, aes(x = Treatment, y = fit, group = Level, shape = Treatment),cex = 1.4,position=position_dodge(width=0.5),
                         colour="black", fill="white",stroke=1.2)+
              scale_shape_manual(values=c(21,16),
                     labels=c("Burned","Unburned"))+
              xlab("Samples")+  ylab("Sampling completeness (%)")+
              facet_wrap(~Level, switch = "x", scales = "free_x")+
              theme(legend.text = element_text(size=15))+
              theme(legend.title = element_text(size=15))+
              theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"))+
              theme(panel.background=element_rect(fill="white"),
                    panel.margin = unit(0, "lines"),
                    strip.background = element_blank(),
                    panel.grid.major.x=element_blank(),
                    panel.grid.major.y=element_line(colour="gray70"),
                    panel.grid.minor=element_blank(),
                    panel.border=element_rect(color="black",fill=F,size=1),
                    text=element_text(size=15),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.title=element_blank(),
                    legend.key = element_blank())

g1



ggsave("Fig2.svg", plot = g1, device = "svg", path = "Results/UpdatedFigs/", width = 14.255, height = 10.515, units = "cm")





























