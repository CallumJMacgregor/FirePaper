###############################################################################################
####   Script for analysing and plotting NMDS of moth, plant and interaction communities   ####
###############################################################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","arm","MASS","scales","AICcmodavg","svglite","effects","plyr","gridExtra","reshape2","vegan")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### load up Callum's custom set of functions
k <- c("CheckResidsFunction.R","CheckConvergenceFunction.R","MultiplotFunction.R")
lapply(k,source)

#### first, let's do the moths

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

summary(dframe2)


# now introduce the season order data

# we also want this variable to be ordered so that we can look for any trend over time
# I've just crudely put this into a text file that we can import and merge in

dframe2 <- merge(SeasonOrder,dframe2)
summary(dframe2$ExactSeason)
summary(dframe2$Order)


# now we're ready to start analysing 
summary(dframe2)

# we want to collapse this one step further so that we have how many individuals of each species in each network (n = 18)

dframe3 <- ddply(dframe2, .(Family_Species, Treatment, Season, ExactSeason), numcolwise(sum))

summary(dframe3)
summary(dframe3$ExactSeason)

# now we need to cast this so that each sample is a single row, with a column for each insect species

dframe3.cast <- dcast(dframe3, Treatment + Season + ExactSeason ~ Family_Species,
                      value.var = "Count",
                      fun.aggregate = sum)

summary(dframe3.cast)
dframe3.cast[,1:3]


# and create a community version

community.m <- subset(dframe3.cast, select=-c(Treatment,Season,ExactSeason))

## now we're ready to do community stats on this

# NMDS

NMDS_moths <- metaMDS(community.m, k=2, trymax=100)

plot(NMDS_moths)

ordiplot(NMDS_moths, type = "n")
orditorp(NMDS_moths,display = "sites", cex=1.25,air=0.01)

ordiellipse(NMDS_moths, group=dframe3.cast$Treatment, show="NoFire")
ordiellipse(NMDS_moths, group=dframe3.cast$Treatment, show="Fire", col="red")


AdMoths <- adonis(community.m ~ Treatment*Season,
                  data = dframe3.cast, 
                  strata = dframe3.cast$ExactSeason,
                  perm=1e5)

AdMoths$aov.tab

# now prepare the same thing for plants


### Plants ####

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe4<-read.table("Data/PlantTransects.txt", header=TRUE)

dframe4$Year <- factor(dframe4$Year)
dframe4$Treatment <- relevel(dframe4$Treatment, "NoFire")
summary(dframe4)


# add in the exact season

dframe4$ExactSeason <- do.call(paste, c(dframe4[c("Season","Year")], sep = "_"))
dframe4$ExactSeason <- as.factor(dframe4$ExactSeason)
summary(dframe4$ExactSeason)


dframe4 <- merge(SeasonOrder,dframe4)
summary(dframe4$ExactSeason)
summary(dframe4$Order)

# We don't need the individual TransectID for this analysis, so get rid of it:

dframe4r <- dframe4[,c(1:2, 4:16)]
summary(dframe4r)

# each row currently contains 1 plant species, so add a column for "Count"; 
# this is 1 in every instance except the transect with zero flowers
dframe4r$SpeciesRichness <- ifelse(dframe4r$PlantSpecies=="none",0,1)


# now we're ready to start analysing 
summary(dframe4r)

# we want to collapse this one step further so that we have how the average cover of each species in each network (n = 18)

dframe5 <- ddply(dframe4r, .(PlantSpecies, Treatment, Season, ExactSeason), numcolwise(mean))

summary(dframe5)
summary(dframe5$ExactSeason)

# now we need to cast this so that each sample is a single row, with a column for each insect species

dframe5.cast <- dcast(dframe5, Treatment + Season + ExactSeason ~ PlantSpecies,
                      value.var = "PlantCoverage",
                      fun.aggregate = sum)

summary(dframe5.cast)
dframe5.cast[,1:3]


# and create a community version

community.p <- subset(dframe5.cast, select=-c(Treatment,Season,ExactSeason))

## now we're ready to do community stats on this

# NMDS

NMDS_plants <- metaMDS(community.p, k=2, trymax=100)

plot(NMDS_plants)

ordiplot(NMDS_plants, type = "n")
orditorp(NMDS_plants,display = "sites", cex=1.25,air=0.01)

ordiellipse(NMDS_plants, group=dframe5.cast$Treatment, show="NoFire")
ordiellipse(NMDS_plants, group=dframe5.cast$Treatment, show="Fire", col="red")



AdPlants <- adonis(community.p ~ Treatment*Season,
                  data = dframe5.cast, 
                  strata = dframe5.cast$ExactSeason,
                  perm=1e5)

AdPlants$aov.tab



### interactions


# finally let's try to do this with interaction data (slightly more complex)

# we need to create an edgelist of interactions for the data from each network

# let's start just by melting it

summary(dframe1)

dframe6 <- melt(dframe1,
                id=c("SampleID","Family_Species","Date","Site","SlideNumber","PollenCount","Treatment","SamplingDay","Sample","Season","Month","Year"))
  
  
# now bind together the two species-name columns to get an interaction-name column

dframe6$Interaction <- as.factor(paste(dframe6$Family_Species,dframe6$variable, sep = "-"))

# create an abundance for each interaction, then collapse the dataframe to count these abundances in each sample
dframe6$Count <- ifelse(dframe6$value == 0, 0, 1)  

dframe7 <- ddply(dframe6, .(Treatment, Sample, Season, Year, Interaction), numcolwise(sum))

summary(dframe7)

dframe7 <- dframe7[,c("Treatment","Sample","Season","Year","Interaction","Count")]

# add in the exact season

dframe7$ExactSeason <- do.call(paste, c(dframe7[c("Season","Year")], sep = "_"))
dframe7$ExactSeason <- as.factor(dframe7$ExactSeason)
summary(dframe7$ExactSeason)


dframe7 <- merge(SeasonOrder,dframe7)
summary(dframe7$ExactSeason)
summary(dframe7$Order)



# now we're ready to start analysing 
summary(dframe7)

# we want to collapse this one step further so that we have the frequency of each interaction in each network (n = 18)

dframe8 <- ddply(dframe7, .(Interaction, Treatment, Season, ExactSeason), numcolwise(sum))

summary(dframe8)

# now we need to cast this so that each sample is a single row, with a column for each insect species

dframe8.cast <- dcast(dframe8, Treatment + Season + ExactSeason ~ Interaction,
                      value.var = "Count",
                      fun.aggregate = sum)

summary(dframe8.cast)
dframe8.cast[,1:3]


# and create a community version

community.i <- subset(dframe8.cast, select=-c(Treatment,Season,ExactSeason))

## now we're ready to do community stats on this

# NMDS

NMDS_inters <- metaMDS(community.i, k=2, trymax=100)

plot(NMDS_inters)

ordiplot(NMDS_inters, type = "n")
orditorp(NMDS_inters,display = "sites", cex=1.25,air=0.01)

ordiellipse(NMDS_inters, group=dframe8.cast$Treatment, show="NoFire")
ordiellipse(NMDS_inters, group=dframe8.cast$Treatment, show="Fire", col="red")



AdInters <- adonis(community.i ~ Treatment*Season,
                   data = dframe8.cast, 
                   strata = dframe8.cast$ExactSeason,
                   perm=1e5)

AdInters$aov.tab


## finally, save all three plots as a panel figure

svg("Results/NMDSfig.svg", width = 4, height = 12, bg = "white")

par( mfrow = c(3,1), oma = c(0,0,2,0))


ordiplot(NMDS_moths, type = "n")
orditorp(NMDS_moths,display = "sites", cex=1.25,air=0.01)

ordiellipse(NMDS_moths, group=dframe3.cast$Treatment, show="NoFire")
ordiellipse(NMDS_moths, group=dframe3.cast$Treatment, show="Fire", col="red")



ordiplot(NMDS_plants, type = "n")
orditorp(NMDS_plants,display = "sites", cex=1.25,air=0.01)

ordiellipse(NMDS_plants, group=dframe5.cast$Treatment, show="NoFire")
ordiellipse(NMDS_plants, group=dframe5.cast$Treatment, show="Fire", col="red")




ordiplot(NMDS_inters, type = "n")
orditorp(NMDS_inters,display = "sites", cex=1.25,air=0.01)

ordiellipse(NMDS_inters, group=dframe8.cast$Treatment, show="NoFire")
ordiellipse(NMDS_inters, group=dframe8.cast$Treatment, show="Fire", col="red")



dev.off()






