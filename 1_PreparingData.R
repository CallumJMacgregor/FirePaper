####################################################
####   Script for preparing data for analysis   ####
####################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("reshape2","plyr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### Insects and pollen

### read in the data - this is the raw data as input by Paula. 
# I have manually added a column called 'SampleID' with a unique number for each row of the dataframe;
# I have renamed "Family/Species" to "Family_Species", 
# and added _1 to the first PollenType and PollenNumber column for consistency
# Finally, I've put a 0 in "PollenNumber_1" for insects with no pollen - this is important
dframe1<-read.csv("Data/BanzaDataNoct.csv", header=TRUE)
#dframe1<-read.csv(file.choose()) 


summary(dframe1) # Check it's imported correctly

### Remove redundant PollenTypes 14 & 15
dframe1 <- dframe1[,1:33]
summary(dframe1)

### Make binary "Lepidoptera" variable a factor
dframe1$Lepidoptera <- factor(dframe1$Lepidoptera)

### Melt the dataframe
names(dframe1) # check what columns you have and what they're called

# create a dframe with all the family/species names for each sample
melt1 <- melt(dframe1[,c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount",
                         "PollenType_1","PollenType_2","PollenType_3","PollenType_4",
                         "PollenType_5","PollenType_6","PollenType_7","PollenType_8",
                         "PollenType_9","PollenType_10","PollenType_11","PollenType_12","PollenType_13")],
                      id=c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount"))

#label each row of the new dframe and remove the old column headings
melt1$level.order <- seq(1, length(melt1$SampleID),1)
melt1 <- melt1[,-8]

# create a dframe with all the pollen counts for each sample
melt2 <- melt(dframe1[,c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount",
                         "PollenNumber_1","PollenNumber_2","PollenNumber_3","PollenNumber_4",
                         "PollenNumber_5","PollenNumber_6","PollenNumber_7","PollenNumber_8",
                         "PollenNumber_9","PollenNumber_10","PollenNumber_11","PollenNumber_12","PollenNumber_13")],
              id=c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount"))

# label each row of the new dframe and remove the old column headings
melt2$level.order <- seq(1, length(melt2$SampleID),1)
melt2 <- melt2[,-8]  #trim out the ID column - this is important for merging back together

# rename the 'value' column on the second dframe so it merges correctly
colnames(melt2) <- c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount", "value2","level.order")

# merge the dframes back together
dframe1.melt <- merge(melt1, melt2)

# remove the row labels as we are finished with them
dframe1.melt <- dframe1.melt[,-8]

# remove all the rows with NAs, keeping a single row for any moths with zero pollen (because of the zeroes in PollenNumber_1)
dframe1.melt <- dframe1.melt[complete.cases(dframe1.melt[,9]),]

# rename the remaining columns to original names
colnames(dframe1.melt) <- c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount", "PollenType","PollenNumber")

# view information about the columns
str(dframe1.melt)

# reformat the data into a matrix
matrix1 <- dcast(dframe1.melt, SampleID + Family_Species + Lepidoptera + Date +
                   Site + SlideNumber + PollenCount ~ PollenType,
                 value.var = "PollenNumber",
                 fun.aggregate = sum)



# this adds an extra blank column (I don't know why!) and a column for zero; let's delete them now
matrix1 <- matrix1[,c(1:7,10:length(matrix1))]


### create a variable for fire/no fire.
### This tells R to search for "NF" in each row of 'Site', to label rows containing "NF" as NoFire, and all other rows as Fire
matrix1$Treatment <- ifelse(grepl("NF",matrix1$Site),"NoFire","Fire")
matrix1$Treatment <- factor(matrix1$Treatment)

### create a variable for sample
# this will be important when we come to do network analysis
# first tell R to give each date a unique value
matrix1$SamplingDay <- factor(matrix1$Date,
                              levels = levels(matrix1$Date),
                              labels = 1:21)

# then create a new column including both Site and Date information, so each sample at each site has a unique level
matrix1$Sample <- do.call(paste, c(matrix1[c("Site","SamplingDay")], sep = "_"))


# create a variable for season and one for month
# four seasons of three months each:
# "winter" - Dec-Feb, "spring" - Mar-May, "Summer" - Jun-Aug, "Autumn" - Sep-Nov

# this code searches in the date column for the month and assigns each month to the appropriate categories

matrix1$Season <- ifelse(grepl("/11/",matrix1$Date),"Autumn",
                         ifelse(grepl("/12/",matrix1$Date),"Autumn",
                                ifelse(grepl("/01/",matrix1$Date),"Winter",
                                       ifelse(grepl("/02/",matrix1$Date),"Winter",
                  ifelse(grepl("/03/",matrix1$Date),"Winter",
                         ifelse(grepl("/04/",matrix1$Date),"Spring",
                                ifelse(grepl("/05/",matrix1$Date),"Spring",
                                       ifelse(grepl("/06/",matrix1$Date),"Spring",
                  ifelse(grepl("/07/",matrix1$Date),"Summer",
                         ifelse(grepl("/08/",matrix1$Date),"Summer",
                                ifelse(grepl("/09/",matrix1$Date),"Summer",
                                       ifelse(grepl("/10/",matrix1$Date),"Autumn",
                                              "Fail"))))))))))))
matrix1$Season <- factor(matrix1$Season)

matrix1$Month <- ifelse(grepl("/11/",matrix1$Date),"Nov",
                         ifelse(grepl("/12/",matrix1$Date),"Dec",
                                ifelse(grepl("/01/",matrix1$Date),"Jan",
                                       ifelse(grepl("/02/",matrix1$Date),"Feb",
                                              ifelse(grepl("/03/",matrix1$Date),"Mar",
                                                     ifelse(grepl("/04/",matrix1$Date),"Apr",
                                                            ifelse(grepl("/05/",matrix1$Date),"May",
                                                                   ifelse(grepl("/06/",matrix1$Date),"Jun",
                                                                          ifelse(grepl("/07/",matrix1$Date),"Jul",
                                                                                 ifelse(grepl("/08/",matrix1$Date),"Aug",
                                                                                        ifelse(grepl("/09/",matrix1$Date),"Sep",
                                                                                               ifelse(grepl("/10/",matrix1$Date),"Oct",
                                                                                                      "Fail"))))))))))))
matrix1$Month<-ordered(matrix1$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))


matrix1$Year <- ifelse(grepl("/04/2013",matrix1$Date),"1",
                       ifelse(grepl("/05/2013",matrix1$Date),"1",
                              ifelse(grepl("/06/2013",matrix1$Date),"1",
                                     ifelse(grepl("/07/2013",matrix1$Date),"1",
                                            ifelse(grepl("/08/2013",matrix1$Date),"1",
                                                   ifelse(grepl("/09/2013",matrix1$Date),"1",
                ifelse(grepl("/10/2014",matrix1$Date),"3",
                              ifelse(grepl("/11/2014",matrix1$Date),"3",
                                     ifelse(grepl("/12/2014",matrix1$Date),"3",
                                            ifelse(grepl("/01/2015",matrix1$Date),"3",
                                                   ifelse(grepl("/02/2015",matrix1$Date),"3",
                                                          ifelse(grepl("/03/2015",matrix1$Date),"3",
                                                                 ifelse(grepl("/04/2015",matrix1$Date),"3",
                                                                        ifelse(grepl("/05/2015",matrix1$Date),"3",
                                                    "2"))))))))))))))
matrix1$Year <- factor(matrix1$Year)


# check it's worked - there should be no values for Fail
summary(matrix1$Season)
summary(matrix1$Month)
summary(matrix1$Year)


# finally, reorder the columns to put pollen variables after the new variables
names(matrix1) # check what columns we have and what order
matrix1 <- matrix1[,c(1:7,77:82,8:76)] # tell R what order to put them in
names(matrix1) # check it's worked


# we want to check and remove "environmental contamination" - i.e. Pinus, Briza and Cupressus
# First, to check the quantity of it:

# create a reduced dataframe
matrix1c <- matrix1[,c(4:5,10:length(matrix1))]
summary(matrix1c)

# create a dummy variable that should be identical for every row of data
matrix1c$dummy <- ifelse(grepl("Fail",matrix1c$Season),"Fail","Pass")
matrix1c$dummy <- factor(matrix1c$dummy)
summary(matrix1c$dummy)


# create a dataframe with the total quantity of each pollen type
matrix1ce <- ddply(matrix1c, .(dummy), numcolwise(sum))
summary(matrix1ce)

# calculate the percentage of total pollen grains that are (a) Pinus and (b) Cupressus
matrix1ce$Total <- rowSums(matrix1ce[,2:length(matrix1ce)])

matrix1ce$percPinus <- matrix1ce$`Pinus sp.`*100/matrix1ce$Total
matrix1ce$percCupressus <- matrix1ce$`Cupressus sp.`*100/matrix1ce$Total
matrix1ce$percBriza <- matrix1ce$`Briza maxima`*100/matrix1ce$Total

# view percentages
str(matrix1ce$percPinus)
str(matrix1ce$percCupressus)
str(matrix1ce$percBriza)


# now remove these columns from the main working dataframe
matrix1$`Pinus sp.` <- NULL
matrix1$`Cupressus sp.` <- NULL
matrix1$`Briza maxima` <- NULL


# finally, we want to restrict the dataset to Lepidoptera only
matrix1L <- matrix1[matrix1$Lepidoptera==1,]
matrix1L <- matrix1L[,-3]   

# N.B. at this point any typing errors in species names will become obvious as they will be assigned an additional column.
# Therefore at this point, go back to the raw data, find these typing errors, correct them, and run this script again.


# there is an additional check for this later in this script, so we save the remaining steps for later







### plants

### read in the data - this is the raw data as input by Paula. 
# I have manually added a column called 'TransectID' with a unique number for each row of the dataframe;
# I have renamed some of the columns, 
# and added _1 to the first Plant and PlantCoverage column for consistency
# Finally, I've put a 0 in "PlantCoverage_1" for transects with no flowers - this is important
dframe2<-read.csv("Data/PlantTransects.csv", header=TRUE)

summary(dframe2) # Check it's imported correctly

### Remove redundant Plant columns 9 & 10
dframe2 <- dframe2[,1:22]
summary(dframe2)


### Melt the dataframe
names(dframe2) # check what columns you have and what they're called

# create a dframe with all the family/species names for each sample
melt3 <- melt(dframe2[,c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm",
                         "Plant_1","Plant_2","Plant_3","Plant_4",
                         "Plant_5","Plant_6","Plant_7","Plant_8")],
              id=c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm"))

#label each row of the new dframe and remove the old column headings
melt3$level.order <- seq(1, length(melt3$TransectID),1)
melt3 <- melt3[,-7]

# create a dframe with all the pollen counts for each sample
melt4 <- melt(dframe2[,c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm",
                         "PlantCoverage_1","PlantCoverage_2","PlantCoverage_3","PlantCoverage_4",
                         "PlantCoverage_5","PlantCoverage_6","PlantCoverage_7","PlantCoverage_8")],
              id=c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm"))

# label each row of the new dframe and remove the old column headings
melt4$level.order <- seq(1, length(melt4$TransectID),1)
melt4 <- melt4[,-7]  #trim out the ID column - this is important for merging back together

# rename the 'value' column on the second dframe so it merges correctly
colnames(melt4) <- c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm","value2","level.order")

# merge the dframes back together
dframe2.melt <- merge(melt3, melt4)

# remove the row labels as we are finished with them
dframe2.melt <- dframe2.melt[,-7]

# remove all the rows with NAs, keeping a single row for any moths with zero pollen (because of the zeroes in PollenNumber_1)
dframe2.melt <- dframe2.melt[complete.cases(dframe2.melt[,8]),]

# rename the remaining columns to original names
colnames(dframe2.melt) <- c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm", "PlantSpecies","PlantCoverage")

# view information about the columns
str(dframe2.melt)

# reformat the data into a matrix
matrix3 <- dframe2.melt



### create a variable for fire/no fire.
### This tells R to search for "NF" in each row of 'Site', to label rows containing "NF" as NoFire, and all other rows as Fire
matrix3$Treatment <- ifelse(grepl("NF",matrix3$Site),"NoFire","Fire")
matrix3$Treatment <- factor(matrix3$Treatment)

### create a variable for sample
# this will be important when we come to do network analysis
# first tell R to give each date a unique value
matrix3$SamplingDay <- factor(matrix3$Date,
                              levels = levels(matrix3$Date),
                              labels = 1:21)

# then create a new column including both Site and Date information, so each sample at each site has a unique level
matrix3$Sample <- do.call(paste, c(matrix3[c("Site","SamplingDay")], sep = "_"))


# create a variable for season and one for month
# four seasons of three months each:
# "winter" - Dec-Feb, "spring" - Mar-May, "Summer" - Jun-Aug, "Autumn" - Sep-Nov

# this code searches in the date column for the month and assigns each month to the appropriate categories

matrix3$Season <- ifelse(grepl("/11/",matrix3$Date),"Autumn",
                         ifelse(grepl("/12/",matrix3$Date),"Autumn",
                                ifelse(grepl("/01/",matrix3$Date),"Winter",
                                       ifelse(grepl("/02/",matrix3$Date),"Winter",
                                              ifelse(grepl("/03/",matrix3$Date),"Winter",
                                                     ifelse(grepl("/04/",matrix3$Date),"Spring",
                                                            ifelse(grepl("/05/",matrix3$Date),"Spring",
                                                                   ifelse(grepl("/06/",matrix3$Date),"Spring",
                                                                          ifelse(grepl("/07/",matrix3$Date),"Summer",
                                                                                 ifelse(grepl("/08/",matrix3$Date),"Summer",
                                                                                        ifelse(grepl("/09/",matrix3$Date),"Summer",
                                                                                               ifelse(grepl("/10/",matrix3$Date),"Autumn",
                                                                                                      "Fail"))))))))))))
matrix3$Season <- factor(matrix3$Season)

matrix3$Month <- ifelse(grepl("/11/",matrix3$Date),"Nov",
                        ifelse(grepl("/12/",matrix3$Date),"Dec",
                               ifelse(grepl("/01/",matrix3$Date),"Jan",
                                      ifelse(grepl("/02/",matrix3$Date),"Feb",
                                             ifelse(grepl("/03/",matrix3$Date),"Mar",
                                                    ifelse(grepl("/04/",matrix3$Date),"Apr",
                                                           ifelse(grepl("/05/",matrix3$Date),"May",
                                                                  ifelse(grepl("/06/",matrix3$Date),"Jun",
                                                                         ifelse(grepl("/07/",matrix3$Date),"Jul",
                                                                                ifelse(grepl("/08/",matrix3$Date),"Aug",
                                                                                       ifelse(grepl("/09/",matrix3$Date),"Sep",
                                                                                              ifelse(grepl("/10/",matrix3$Date),"Oct",
                                                                                                     "Fail"))))))))))))
matrix3$Month<-ordered(matrix3$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))


matrix3$Year <- ifelse(grepl("/04/2013",matrix3$Date),"1",
                       ifelse(grepl("/05/2013",matrix3$Date),"1",
                              ifelse(grepl("/06/2013",matrix3$Date),"1",
                                     ifelse(grepl("/07/2013",matrix3$Date),"1",
                                            ifelse(grepl("/08/2013",matrix3$Date),"1",
                                                   ifelse(grepl("/09/2013",matrix3$Date),"1",
                                                          ifelse(grepl("/10/2014",matrix3$Date),"3",
                                                                 ifelse(grepl("/11/2014",matrix3$Date),"3",
                                                                        ifelse(grepl("/12/2014",matrix3$Date),"3",
                                                                               ifelse(grepl("/01/2015",matrix3$Date),"3",
                                                                                      ifelse(grepl("/02/2015",matrix3$Date),"3",
                                                                                             ifelse(grepl("/03/2015",matrix3$Date),"3",
                                                                                                    ifelse(grepl("/04/2015",matrix3$Date),"3",
                                                                                                           ifelse(grepl("/05/2015",matrix3$Date),"3",
                                                                                                                  "2"))))))))))))))
matrix3$Year <- factor(matrix3$Year)

# check it's worked - there should be no values for Fail
summary(matrix3$Season)
summary(matrix3$Month)
summary(matrix3$Year)

matrix3$PlantSpecies <- ifelse(matrix3$PlantCoverage==0,"none",matrix3$PlantSpecies)


# finally, reorder the columns to put pollen variables after the new variables
names(matrix3) # check what columns we have and what order
matrix3 <- matrix3[,c(1:6,9:14,7:8)] # tell R what order to put them in
names(matrix3) # check it's worked



# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(matrix3, "Data\\PlantTransects.txt", sep="\t", row.names=FALSE)

# N.B. at this point any typing errors in species names will become obvious as they will be assigned an additional column.
# Therefore at this point, go back to the raw data, find these typing errors, correct them, and run this script again.







### basic summary statistics
# we want a selection of very basic statistics to cite at the start of the results section:

# total number of moths caught
nrow(matrix1L)

# number of moths identified to at least genus
matrix1s <- matrix1L[,c(1:2)]
matrix1s$Count <- 1
summary(matrix1s)

# save the table at this point for later use
matrix1all <- ddply(matrix1s, .(Family_Species), numcolwise(sum))
write.table(matrix1all, "Data\\AllSpecies.txt", sep="\t", row.names=FALSE)


# now remove anything not IDed to at least genus
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="none")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Unknown micro")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Unknown moth")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Lepidoptera")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Coleophoridae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Crambidae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Elachistidae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Gelechiidae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Geometridae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Gracillariidae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Noctuidae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Plutellidae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Pterophoridae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Pyralidae")
matrix1s <- subset(matrix1s, !matrix1s$Family_Species=="Tortricidae")

# number of rows is number of identified (to at least genus) individuals
nrow(matrix1s)

# number of (identified to at least genus) moth species caught
matrix1ns <- ddply(matrix1s, .(Family_Species), numcolwise(sum))

# number of rows is number of morphotypes IDed to genus or better
nrow(matrix1ns)


# number of moth families caught
# at this point I have manually added the family that each species belongs to in the file "Data\\AllSpecies.txt" 
# I have also recorded the taxonomic level to which each identification has been made
# and made a small number of corrections/updates to species names

# import the new table with this information
dframe5<-read.csv("Data/AllSpeciesEdited.csv", header = T)
summary(dframe5)

summary(dframe5$IdentifiedTo)
nrow(dframe5)

dframe5$NoIdentifiedSp <- ifelse(dframe5$IdentifiedTo=="Species", 1,0)
dframe5$NoMorphotypes <- 1

dframe5i <- ddply(dframe5, .(IdentifiedTo), numcolwise(sum))
dframe5f <- ddply(dframe5, .(Family), numcolwise(sum))

# we want both of these as summary tables
write.table(dframe5i, "Results\\IDlevel.txt", sep="\t", row.names=FALSE)
write.table(dframe5f, "Results\\Families.txt", sep="\t", row.names=FALSE)


### this is also the time to revisit our downstream scripts, as we have corrected the names

matrix1Lc <- merge(dframe5,matrix1L, by.x = "Family_Species", by.y="Family_Species")

# remove extra columns and rename corrected Family_Species column

matrix1Lca <- matrix1Lc[,c(8,2,9:length(matrix1Lc))]

colnames(matrix1Lca)[2] <- "Family_Species"


# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(matrix1Lca, "Data\\MatrixNoct.txt", sep="\t", row.names=FALSE)


# we want to duplicate this with a threshold applied, removing any interactions with <5 pollen grains counted
matrix1Lca.thresh <- matrix1Lca[,13:78]

matrix1Lca.thresh[matrix1Lca.thresh < 5] <- 0

matrix1Lca.thresh <- cbind(matrix1Lca[,1:12],matrix1Lca.thresh)

# output the thresholded table
write.table(matrix1Lca.thresh, "Data\\MatrixNoctThreshold.txt", sep="\t", row.names=FALSE)




# now we want to produce a separate file with the data organised by each sampling session at each site
# we will do this by aggregating samples from the same site and sampling session

matrix1r <- matrix1Lc[,c(5:6,9:length(matrix1Lc))]
summary(matrix1r)


matrix2 <- ddply(matrix1r, .(Site,Date,Treatment,SamplingDay,Sample,Season,Month,Year), numcolwise(sum))

# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(matrix2, "Data\\SamplesNoct.txt", sep="\t", row.names=FALSE)

# again, try adding the threshold
# first make sure all columns are numeric - some have disappeared
for (i in 13:78){
  matrix1Lca.thresh[,i] <- as.numeric(as.character(matrix1Lca.thresh[,i]))
}




matrix2.thresh <- ddply(matrix1Lca.thresh, .(Site,Date,Treatment,SamplingDay,Sample,Season,Month,Year), numcolwise(sum))

# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(matrix2.thresh, "Data\\SamplesNoctThreshold.txt", sep="\t", row.names=FALSE)




# number of plant species recorded
matrix3a <- matrix3
matrix3a$NoTransects <- 1
matrix3s <- ddply(matrix3a, .(PlantSpecies), numcolwise(sum))
matrix3s <- subset(matrix3s, !matrix3s$PlantSpecies=="none")
matrix3s <- subset(matrix3s, !matrix3s$PlantSpecies=="")

matrix3s <- matrix3s[,c(1,6)]

nrow(matrix3s)

# number of plant families recorded
write.table(matrix3s, "Data\\AllPlantSpecies.txt", sep="\t", row.names=FALSE)

# manually added family column and correct a couple of spelling errors
dframe6<-read.csv("Data/AllPlantSpecies.csv", header = T)
summary(dframe6)

dframe6$NoSpecies <- 1

dframe6f <- ddply(dframe6, .(Family), numcolwise(sum))

# we need this for a summary table
write.table(dframe6f, "Results\\PlantFamilies.txt", sep="\t", row.names=FALSE)



### some summary stats

# number and percentage of moths carrying pollen
matrix1Lc$PollenSum <- rowSums(matrix1Lc[,c(19:length(matrix1Lc))])
matrix1Lc$PollenYN <- ifelse(matrix1Lc$PollenSum==0,0,1)

sum(matrix1Lc$PollenYN)
nrow(matrix1Lc)

sum(matrix1Lc$PollenYN)*100/nrow(matrix1Lc)


# species and families of pollen-carriers
matrix1P <- subset(matrix1Lc, matrix1Lc$PollenYN==1)

matrix1PC <- matrix1P[,c(2,85,86)]

matrix1PC$Count <- 1

matrix1PCs <- ddply(matrix1PC, .(CorrectedName), numcolwise(sum))

dframe5s <- dframe5[,c(2,4,5)]

matrix1PCm <- merge(matrix1PCs,dframe5s, by="CorrectedName", all=T)

matrix1PCc <- matrix1PCm[!is.na(matrix1PCm$PollenYN),]
  
nrow(matrix1PCc)*100/nrow(matrix1PCm)

summary(matrix1PCc$IdentifiedTo)

matrix1PCf <- ddply(matrix1PCc, .(Family), numcolwise(sum))

nrow(matrix1PCf)-1




# number and percentage of pollen species
matrix1Pp <- matrix1P[,c(2,19:84)]

matrix1Ppc <- ddply(matrix1Pp, .(CorrectedName), numcolwise(sum))

rownames(matrix1Ppc) <- matrix1Ppc[,1]
matrix1Ppc <- matrix1Ppc[,-1]


pollen <- data.frame(t(matrix1Ppc))

pollen$TotalPollen <- rowSums(pollen)

pollen$PollenSpecies <- rownames(pollen)
pollen <- pollen[,c(length(pollen)-1,length(pollen))]


compare <- merge(pollen,dframe6,by.x = "PollenSpecies", by.y = "PlantSpecies", all=T)

compare[is.na(compare)] <- 0


compare$Class <- as.factor(ifelse(compare$TotalPollen==0,"NoPollen",
                                  ifelse(compare$NoTransects==0,"NoFlower","Both")))

summary(compare$Class)

58*100/69


# repeat with threshold

# number and percentage of moths carrying pollen
matrix1Lca.thresh$PollenSum <- rowSums(matrix1Lca.thresh[,c(13:length(matrix1Lca.thresh))])
matrix1Lca.thresh$PollenYN <- ifelse(matrix1Lca.thresh$PollenSum==0,0,1)

sum(matrix1Lca.thresh$PollenYN)
nrow(matrix1Lca.thresh)

sum(matrix1Lca.thresh$PollenYN)*100/nrow(matrix1Lca.thresh)


# species and families of pollen-carriers
matrix1P.thresh <- subset(matrix1Lca.thresh, matrix1Lca.thresh$PollenYN==1)

matrix1PC.thresh <- matrix1P.thresh[,c(2,79,80)]

matrix1PC.thresh$Count <- 1

matrix1PCs.thresh <- ddply(matrix1PC.thresh, .(Family_Species), numcolwise(sum))

dframe5s <- dframe5[,c(2,4,5)]

matrix1PCm.thresh <- merge(matrix1PCs.thresh,dframe5s, by.x="Family_Species", by.y="CorrectedName", all=T)

matrix1PCc.thresh <- matrix1PCm.thresh[!is.na(matrix1PCm.thresh$PollenYN),]

nrow(matrix1PCc.thresh)*100/nrow(matrix1PCm.thresh)

summary(matrix1PCc.thresh$IdentifiedTo)

matrix1PCf.thresh <- ddply(matrix1PCc.thresh, .(Family), numcolwise(sum))

nrow(matrix1PCf.thresh)-1




# number and percentage of pollen species
matrix1Pp.thresh <- matrix1P.thresh[,c(2,13:78)]

matrix1Ppc.thresh <- ddply(matrix1Pp.thresh, .(Family_Species), numcolwise(sum))

rownames(matrix1Ppc.thresh) <- matrix1Ppc.thresh[,1]
matrix1Ppc.thresh <- matrix1Ppc.thresh[,-1]


pollen.thresh <- data.frame(t(matrix1Ppc.thresh))

pollen.thresh$TotalPollen <- rowSums(pollen.thresh)

pollen.thresh$PollenSpecies <- rownames(pollen.thresh)
pollen.thresh <- pollen.thresh[,c(length(pollen.thresh)-1,length(pollen.thresh))]


compare.thresh <- merge(pollen.thresh,dframe6,by.x = "PollenSpecies", by.y = "PlantSpecies", all=T)

compare.thresh[is.na(compare.thresh)] <- 0


compare.thresh$Class <- as.factor(ifelse(compare.thresh$TotalPollen==0,"NoPollen",
                                  ifelse(compare.thresh$NoTransects==0,"NoFlower","Both")))

summary(compare.thresh$Class)

49*100/69


#################################### development #########################
