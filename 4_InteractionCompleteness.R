################################################################
####   Script for network analysis by season and treatment  ####
################################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("plyr","data.table")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### load up Callum's custom set of functions
f <- c("NetworkFunction.R","SpecAccumFunctions.R")
lapply(f, source)



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

dframe1$ExactSeason <- do.call(paste, c(dframe1[c("Season","Year","Treatment")], sep = "_"))
dframe1$ExactSeason <- as.factor(dframe1$ExactSeason)


### prepare the data for network analysis

# trim off extra columns and reorder
dframe1r <- dframe1[,c(2,7,10,length(dframe1),12:(length(dframe1)-1))]


# break up dataframe per species and per Season-Year-Treatment
# now we have a massive number of dframes, each containing every individual sampled of a given species in a given season/treatment combo

dframes <- split(dframe1r, list(dframe1r$ExactSeason,dframe1r$Family_Species))  # this creates a list of smaller dframes, one for each level of sample
summary(dframes)


# now apply another custom function to extract interaction richness and extrapolated Chao richness for every species-season-treatment combo
interaction.completeness <- lapply(dframes, interaction.complete, cols=5, threshold=0)

# this function only gets applied more than the threshold number of individuals of a species are sampled, but here we set threshold=0 so it will work everywhere
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

# for the sake of argument, let's also do a weighted interaction completeness - more attention is paid to more generalist insects?

weighted.mean(ICgood$completeness, ICgood$nchao, na.rm=T)


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


### quick test for bias in sampling completeness
means.full$Treatment <- factor(ifelse(grepl("NoFire",means.full$ExactSeason),"NoFire","Fire"))


plot(WeightedMean ~ Treatment, data = means.full)


model1 <- glm(WeightedMean ~ Treatment,
              data = means.full)

summary(model1)
drop1(model1, test="Chi")

# no significant effect of Treatment on sampling completeness


##################### development ##############################
