##################################################
####   Script for network analysis by sample  ####
##################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","bipartite","plyr","svglite","gridExtra")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### load up Callum's custom set of functions
f <- c("NetworkFunction.R")
lapply(f, source)



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

### prepare the data for network analysis

dframe1$SeasonTreatment <- do.call(paste, c(dframe1[c("Season","Year","Treatment")], sep = "_"))  

# select columns in specific order: Species, Network ID, Pollen
dframe1r <- dframe1[,c(2,9,13:78)]
dframe1s <- dframe1[,c(2,79,13:78)]


### pollen transport analysis

### per sample networks

# summarise pollen transport for each insect species within each sample
# this aggregates all pollen grains carried by an insect species for each plant species within each sample

dframe2 <- ddply(dframe1r, .(Family_Species,Sample), numcolwise(sum))

summary(dframe2)
# create a suitable matrix from the dframe
dframe2m <- ddply(dframe2, .(Family_Species), numcolwise(sum))
rownames(dframe2m) <- dframe2m[,1]
dframe2m <- dframe2m[,-1]

### plot a full-system network
plotweb(dframe2m)
# networklevel(dframe3m)
# trying to run networklevel() on this full network causes my computer to crash, so I've concealed it behind a # to stop it from running every time!


### that's fine, but we want networks (and descriptors) for each sample

# split each sample into a separate dataframe
dframesPT <- split(dframe2, list(dframe2$Sample))  # this creates a list of smaller dframes, one for each level of sample
summary(dframesPT)
# for example...
summary(dframesPT[1]) # the first dframe in the list is site F1 on sampling day 10
dframesPT[1]

### using a loop and a custom function called 'network', ...

### prepare the dataframes for network analysis and run it (this command will take a few minutes to run)
# the command will print "Fail" each time it detects a dataframe with too little data to produce a network
metricsPT <- lapply(dframesPT, network, index = c("linkage density","vulnerability","generality","niche overlap"))
summary(metricsPT)


###

### we now have a list of dataframes, each one containing the network metricsPT of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
metricsPT.merge <- do.call("cbind", metricsPT)    # merge the data with one sample per column
colnames(metricsPT.merge) <- names(metricsPT)     # assign the sample names to each column
metricsPT.merge <- data.frame(t(metricsPT.merge))
metricsPT.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
dframeP<-read.table("Data/SamplesNoct.txt", header=TRUE)
summary(dframeP)

# we'll need a bit more data later on, so save this in a separate dframe
dframePS <- dframeP[,1:8]

# we don't need any of the pollen data, so:
dframeP <- dframeP[,1:8]

# but we do need a Treatment variable, so:
dframeP$Treatment <- ifelse(grepl("NF",dframeP$Site),"NoFire","Fire")
dframeP$Treatment <- factor(dframeP$Treatment)

# for merge to work, we need to set the row names
rownames(dframeP) <- dframeP$Sample

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metricsPT.full <- merge(dframeP, metricsPT.merge, by=0)


## at this point we also want to generate robustness estimates
# we are not simply using networklevel to do this (like for the other metrics)
# because robustness is simulated, and networklevel doesn't appear to run enough simulations as default for the mean to stabilise

robustPT <- lapply(dframesPT, bootstrap_robustness, repeats = 5)
summary(robustPT)


### we now have a list of dataframes, each one containing the network metricsPT of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
robustPT.merge <- do.call("rbind", robustPT)    # merge the data with one sample per column
robustPT.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metricsPT.robust <- merge(metricsPT.full, robustPT.merge, by.x=1,by.y=0)




### before we can analyse this data, we need to remove any samples which have failed

# subset the columns to keep only those which do not have "Fail" values
metricsPT.good <- metricsPT.robust[ which( ! metricsPT.robust$linkage.density %in% "Network failed") , ]

# remove the duplicate 'row.names' column - this information is the same as 'Sample'
metricsPT.good <- metricsPT.good[,-1]


### hooray! we now have network metrics for each sample with sufficient data to produce a network


### output the dataframe to a .txt file to use in downstream analysis
write.table(metricsPT.good, "Data\\NetworkmetricsPTbySample.txt", sep="\t", row.names=FALSE)



### per season networks

# summarise pollen transport for each insect species within each season
# this aggregates all pollen grains carried by an insect species for each plant species within each sample



dframe3 <- ddply(dframe1s, .(Family_Species,SeasonTreatment), numcolwise(sum))

# create a suitable matrix from the dframe
dframe3m <- ddply(dframe3, .(Family_Species), numcolwise(sum))
rownames(dframe3m) <- dframe3m[,1]
dframe3m <- dframe3m[,-1]

### plot a full-system network
plotweb(dframe3m)
# networklevel(dframe3m)
# trying to run networklevel() on this full network causes my computer to crash, so I've concealed it behind a # to stop it from running every time!


### that's fine, but we want networks (and descriptors) for each sample

# split each sample into a separate dataframe
dframesPTS <- split(dframe3, list(dframe3$SeasonTreatment))  # this creates a list of smaller dframes, one for each level of sample
summary(dframesPTS)
# for example...
summary(dframesPTS[1]) # the first dframe in the list is Autumn of year 2
dframesPTS[1]


### using a loop and a custom function called 'network', ...

### prepare the dataframes for network analysis and run it (this command will take a few minutes to run)
# the command will print "Fail" each time it detects a dataframe with too little data to produce a network
metricsPTS <- lapply(dframesPTS, network, index = c("linkage density","vulnerability","generality","niche overlap"))
summary(metricsPTS)


###

### we now have a list of dataframes, each one containing the network metricsPT of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
metricsPTS.merge <- do.call("cbind", metricsPTS)    # merge the data with one sample per column
colnames(metricsPTS.merge) <- names(metricsPTS)     # assign the sample names to each column
metricsPTS.merge <- data.frame(t(metricsPTS.merge))
metricsPTS.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
summary(dframePS)

# we'll need to replicate the SeasonTreatment variable
dframePS$SeasonTreatment <- do.call(paste, c(dframePS[c("Season","Year","Treatment")], sep = "_")) 
dframePS$SeasonTreatment <- factor(dframePS$SeasonTreatment)

# we need to collapse this down into the same 18 rows as the metrics
dframePSc <- ddply(dframePS, .(Treatment,Season,Year,SeasonTreatment), numcolwise(sum))

# we don't need the SamplingDay variable any more
dframePSc <- dframePSc[,-5]


# for merge to work, we need to set the row names
rownames(dframePSc) <- dframePSc$SeasonTreatment


# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metricsPTS.full <- merge(dframePSc, metricsPTS.merge, by=0)


## at this point we also want to generate robustness estimates
# we are not simply using networklevel to do this (like for the other metrics)
# because robustness is simulated, and networklevel doesn't appear to run enough simulations as default for the mean to stabilise

robustPTS <- lapply(dframesPTS, bootstrap_robustness, repeats = 1000)
summary(robustPTS)


### we now have a list of dataframes, each one containing the network metricsPT of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
robustPTS.merge <- do.call("rbind", robustPTS)    # merge the data with one sample per column
robustPTS.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metricsPTS.robust <- merge(metricsPTS.full, robustPTS.merge, by.x=1,by.y=0)


### before we can analyse this data, we need to remove any samples which have failed

# subset the columns to keep only those which do not have "Fail" values
metricsPTS.good <- metricsPTS.robust[ which( ! metricsPTS.robust$linkage.density %in% "Fail") , ]

# remove the duplicate 'row.names' column - this information is the same as 'Sample'
metricsPTS.good <- metricsPTS.good[,-1]


### hooray! we now have network metrics for each sample with sufficient data to produce a network


### output the dataframe to a .txt file to use in downstream analysis
write.table(metricsPTS.good, "Data\\NetworkmetricsPTbySeason.txt", sep="\t", row.names=FALSE)



degPTS <- lapply(dframesPTS, degreedist)

degPTS.merge <- do.call("cbind", degPTS)    # merge the data with one sample per column
colnames(degPTS.merge) <- names(degPTS)     # assign the sample names to each column
degPTS.merge <- data.frame(t(degPTS.merge))
degPTS.merge



### flower-visitor analysis

# change each interaction to a 1 for flower-visitor analysis
dframe1f <- dframe1r
dframe1f[,3:68][dframe1f[,3:68] > 0] <- 1


# summarise plant-insect interactions for each insect species within each sample
# this produces semi-quantitative data - e.g. if two individuals of an insect species
# interact with a plant species in a sample,
# the entry for that interaction within that sample will be 2

dframe3 <- ddply(dframe1f, .(Family_Species,Sample), numcolwise(sum))
  
# create a suitable matrix from the dframe
dframe3m <- ddply(dframe3, .(Family_Species), numcolwise(sum))
rownames(dframe3m) <- dframe3m[,1]
dframe3m <- dframe3m[,-1]

### plot a full-system network
plotweb(dframe3m)
# networklevel(dframe3m)
# trying to run networklevel() on this full network causes my computer to crash, so I've concealed it behind a # to stop it from running every time!


### that's fine, but we want networks (and descriptors) for each sample

# split each sample into a separate dataframe
dframesFV <- split(dframe3, list(dframe3$Sample))  # this creates a list of smaller dframes, one for each level of sample
summary(dframesFV)
# for example...
summary(dframesFV[1]) # the first dframe in the list is site F1 on sampling day 10
dframesFV[1]

### using a loop and a custom function called 'network', ...

### prepare the dataframes for network analysis and run it (this command will take a few minutes to run)
# the command will print "Fail" each time it detects a dataframe with too little data to produce a network
metricsFV <- lapply(dframesFV, network, index = c("linkage density","vulnerability","generality","niche overlap"))
summary(metricsFV)

###

### we now have a list of dataframes, each one containing the network metricsFV of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
metricsFV.merge <- do.call("cbind", metricsFV)    # merge the data with one sample per column
colnames(metricsFV.merge) <- names(metricsFV)     # assign the sample names to each column
metricsFV.merge <- data.frame(t(metricsFV.merge))
metricsFV.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metricsFV.full <- merge(dframeP, metricsFV.merge, by=0)



## at this point we also want to generate robustness estimates
# we are not simply using networklevel to do this (like for the other metrics)
# because robustness is simulated, and networklevel doesn't appear to run enough simulations as default for the mean to stabilise

robustFV <- lapply(dframesFV, bootstrap_robustness, repeats = 5)
summary(robustFV)


### we now have a list of dataframes, each one containing the network metricsPT of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
robustFV.merge <- do.call("rbind", robustFV)    # merge the data with one sample per column
robustFV.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metricsFV.robust <- merge(metricsFV.full, robustFV.merge, by.x=1,by.y=0)



### before we can analyse this data, we need to remove any samples which have failed

# subset the columns to keep only those which do not have "Fail" values
metricsFV.good <- metricsFV.robust[ which( ! metricsFV.robust$linkage.density %in% "Network failed") , ]
  
# remove the duplicate 'row.names' column - this information is the same as 'Sample'
metricsFV.good <- metricsFV.good[,-1]



### output the dataframe to a .txt file to use in downstream analysis
write.table(metricsFV.good, "Data\\NetworkmetricsFVbySample.txt", sep="\t", row.names=FALSE)




### per season networks


# change each interaction to a 1 for flower-visitor analysis
dframe1fs <- dframe1s
dframe1fs[,3:length(dframe1fs)][dframe1fs[,3:length(dframe1fs)] > 0] <- 1

# summarise pollen transport for each insect species within each season
# this aggregates all pollen grains carried by an insect species for each plant species within each sample


dframe4 <- ddply(dframe1fs, .(Family_Species,SeasonTreatment), numcolwise(sum))

# create a suitable matrix from the dframe
dframe4m <- ddply(dframe4, .(Family_Species), numcolwise(sum))
rownames(dframe4m) <- dframe4m[,1]
dframe4m <- dframe4m[,-1]

### plot a full-system network
plotweb(dframe4m)
# networklevel(dframe3m)
# trying to run networklevel() on this full network causes my computer to crash, so I've concealed it behind a # to stop it from running every time!


### that's fine, but we want networks (and descriptors) for each sample

# split each sample into a separate dataframe
dframesFVS <- split(dframe4, list(dframe4$SeasonTreatment))  # this creates a list of smaller dframes, one for each level of sample
summary(dframesFVS)
# for example...
summary(dframesFVS[1]) # the first dframe in the list is site F1 on sampling day 10
dframesFVS[1]


### using a loop and a custom function called 'network', ...

### prepare the dataframes for network analysis and run it (this command will take a few minutes to run)
# the command will print "Fail" each time it detects a dataframe with too little data to produce a network
metricsFVS <- lapply(dframesFVS, network, index = c("linkage density","vulnerability","generality","niche overlap"))
summary(metricsFVS)

###

### we now have a list of dataframes, each one containing the network metricsPT of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
metricsFVS.merge <- do.call("cbind", metricsFVS)    # merge the data with one sample per column
colnames(metricsFVS.merge) <- names(metricsFVS)     # assign the sample names to each column
metricsFVS.merge <- data.frame(t(metricsFVS.merge))
metricsFVS.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used dataframe:
# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metricsFVS.full <- merge(dframePSc, metricsFVS.merge, by=0)


## at this point we also want to generate robustness estimates
# we are not simply using networklevel to do this (like for the other metrics)
# because robustness is simulated, and networklevel doesn't appear to run enough simulations as default for the mean to stabilise

robustFVS <- lapply(dframesFVS, bootstrap_robustness, repeats = 1000)
summary(robustFVS)


### we now have a list of dataframes, each one containing the network metricsPT of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
robustFVS.merge <- do.call("rbind", robustFVS)    # merge the data with one sample per column
robustFVS.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metricsFVS.robust <- merge(metricsFVS.full, robustFVS.merge, by.x=1,by.y=0)


### before we can analyse this data, we need to remove any samples which have failed

# subset the columns to keep only those which do not have "Fail" values
metricsFVS.good <- metricsFVS.robust[ which( ! metricsFVS.robust$linkage.density %in% "Fail") , ]

# remove the duplicate 'row.names' column - this information is the same as 'Sample'
metricsFVS.good <- metricsFVS.good[,-1]


### hooray! we now have network metrics for each sample with sufficient data to produce a network


### output the dataframe to a .txt file to use in downstream analysis
write.table(metricsFVS.good, "Data\\NetworkmetricsFVbySeason.txt", sep="\t", row.names=FALSE)



### degree distribution

# with all that done, we want to do a quick test on the raw interaction data
# looking at the frequency distribution of 'degree' (the number of interactions per species)

# let's return to that original data

# select columns in specific order: Species, Treatment, Pollen
dframe1dd <- dframe1[,c(2,7,13:78)]


# summarise pollen transport for each insect species within each treatment
# this aggregates all pollen grains carried by an insect species for each plant species within each treatment (across all samples)

dframeDD <- ddply(dframe1dd, .(Family_Species,Treatment), numcolwise(sum))

# we don't need the quantitative information here

dframeDD[,3:68][dframeDD[,3:68] > 0] <- 1

# now we calculate the degree for each insect species in each row

dframeDD$Degree <- rowSums(dframeDD[,3:68])

# now split it into two sub-frames

dframeDDb <- dframeDD[which(dframeDD$Treatment=="Fire"),]
dframeDDu <- dframeDD[which(dframeDD$Treatment=="NoFire"),]

hist(dframeDDb$Degree)
burned.insect.freq <- data.frame(table(dframeDDb$Degree))

hist(dframeDDu$Degree)
unburned.insect.freq <- data.frame(table(dframeDDu$Degree))

# now do a Kolgorov-Smirnov test on the two
# as our hypothesis is that unburned has more interactions than burned, we can make it one-sided
# if significant, unburned degree dist > burned degree dist ==> more interactions per pollinator species in unburned

ks.test(dframeDDb$Degree,dframeDDu$Degree,alternative="gr")


dframeDD$Treatment <- revalue(dframeDD$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

# now we want a histogram to depict this test
# first we need to calculate group means
mu <- ddply(dframeDD, "Treatment", summarise, grp.mean=mean(Degree), grp.se=se(Degree))
head(mu)





g1 <- ggplot(data=dframeDD, aes(x=Degree))+
            geom_histogram(binwidth=1,color="black",fill="white")+
            facet_grid(Treatment ~ .)+
            geom_vline(data=mu, aes(xintercept=grp.mean),color="black",
                       linetype="dashed",size=1)+
  ylab("Frequency (moths)")+
  theme(strip.text = element_text(size = 30),
        strip.background = element_rect(color="white",fill="white"),
        panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text=element_text(size=17,color="black"),
        plot.title=element_text(hjust=0.5))

g1


ggsave("FigS7a.svg", plot = g1, device = "svg", path = "Results/", width = 22, height = 15, units = "cm")


# now to repeat this for plants
# first extract the degrees for each plant species in each treatment

degburned <- data.frame(colSums(dframeDDb[,3:68]))
colnames(degburned) <- "Degree"
degburned$Treatment <- "Burned"

degunburned <- data.frame(colSums(dframeDDu[,3:68]))
colnames(degunburned) <- "Degree"
degunburned$Treatment <- "Unburned"

degplants <- rbind(degburned,degunburned)

# test to see if unburned is greater than burned

ks.test(degburned$Degree,degunburned$Degree,alternative="gr")

# now we want a histogram to depict this test
# first we need to calculate group means
mu2 <- ddply(degplants, "Treatment", summarise, grp.mean=mean(Degree), grp.se=se(Degree))
head(mu2)





g2 <- ggplot(data=degplants, aes(x=Degree))+
  geom_histogram(binwidth=1,color="black",fill="white")+
  facet_grid(Treatment ~ .)+
  geom_vline(data=mu2, aes(xintercept=grp.mean),color="black",
             linetype="dashed",size=1)+
  ylab("Frequency (plants)")+
  theme(strip.text = element_text(size = 30),
        strip.background = element_rect(color="white",fill="white"),
        panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=30),
        axis.text=element_text(size=17,color="black"),
        plot.title=element_text(hjust=0.5))

g2


ggsave("FigS7b.svg", plot = g2, device = "svg", path = "Results/", width = 22, height = 15, units = "cm")



# now let's stitch them into a single figure

m1 <- grid.arrange(g1,g2,ncol=1,nrow=2)

ggsave("FigS9full.svg",plot = m1, device = "svg", path = "Results/UpdatedFigs", width = 22, height = 30, units = "cm")


### by season
# this has worked well but it may be informative to do it separately for each season

# moths
# select columns in specific order: Species, Treatment, Pollen
dframe1dds <- dframe1[,c(2,7,10,13:78)]


# summarise pollen transport for each insect species within each treatment
# this aggregates all pollen grains carried by an insect species for each plant species within each treatment (across all samples)

dframeDDs <- ddply(dframe1dds, .(Family_Species,Season,Treatment), numcolwise(sum))

# we don't need the quantitative information here

dframeDDs[,4:69][dframeDDs[,4:69] > 0] <- 1

# now we calculate the degree for each insect species in each row

dframeDDs$Degree <- rowSums(dframeDDs[,4:69])

# now split it into two sub-frames

dframeDDsb <- dframeDDs[which(dframeDDs$Treatment=="Fire"),]
dframeDDsu <- dframeDDs[which(dframeDDs$Treatment=="NoFire"),]

# and these into four each by season
dframeDDsba <- dframeDDsb[which(dframeDDsb$Season=="Autumn"),]
dframeDDsbsp <- dframeDDsb[which(dframeDDsb$Season=="Spring"),]
dframeDDsbsu <- dframeDDsb[which(dframeDDsb$Season=="Summer"),]
dframeDDsbw <- dframeDDsb[which(dframeDDsb$Season=="Winter"),]

dframeDDsua <- dframeDDsu[which(dframeDDsu$Season=="Autumn"),]
dframeDDsusp <- dframeDDsu[which(dframeDDsu$Season=="Spring"),]
dframeDDsusu <- dframeDDsu[which(dframeDDsu$Season=="Summer"),]
dframeDDsuw <- dframeDDsu[which(dframeDDsu$Season=="Winter"),]


# now start testing
# autumn
hist(dframeDDsba$Degree)
hist(dframeDDsua$Degree)

ks.test(dframeDDsba$Degree,dframeDDsua$Degree,alternative="gr")

# spring
hist(dframeDDsbsp$Degree)
hist(dframeDDsusp$Degree)

ks.test(dframeDDsbsp$Degree,dframeDDsusp$Degree,alternative="gr")

# summer
hist(dframeDDsbsu$Degree)
hist(dframeDDsusu$Degree)

ks.test(dframeDDsbsu$Degree,dframeDDsusu$Degree,alternative="gr")

# winter
hist(dframeDDsbw$Degree)
hist(dframeDDsuw$Degree)

ks.test(dframeDDsbw$Degree,dframeDDsuw$Degree,alternative="gr")



dframeDDs$Treatment <- revalue(dframeDDs$Treatment, c("Fire"="Burned","NoFire"="Unburned"))

# we need to calculate group means
mus <- ddply(dframeDDs, c("Treatment","Season"), summarise, grp.mean=mean(Degree), grp.se=se(Degree))
print(mus)



# now to repeat this for plants
# first extract the degrees for each plant species in each treatment in each season

#autumn
degburnedsa <- data.frame(colSums(dframeDDsba[,4:69]))
colnames(degburnedsa) <- "Degree"
degburnedsa$Treatment <- "Burned"
degburnedsa$Season <- "Autumn"

degunburnedsa <- data.frame(colSums(dframeDDsua[,4:69]))
colnames(degunburnedsa) <- "Degree"
degunburnedsa$Treatment <- "Unburned"
degunburnedsa$Season <- "Autumn"


# spring
degburnedssp <- data.frame(colSums(dframeDDsbsp[,4:69]))
colnames(degburnedssp) <- "Degree"
degburnedssp$Treatment <- "Burned"
degburnedssp$Season <- "Spring"

degunburnedssp <- data.frame(colSums(dframeDDsusp[,4:69]))
colnames(degunburnedssp) <- "Degree"
degunburnedssp$Treatment <- "Unburned"
degunburnedssp$Season <- "Spring"

#summer
degburnedssu <- data.frame(colSums(dframeDDsbsu[,4:69]))
colnames(degburnedssu) <- "Degree"
degburnedssu$Treatment <- "Burned"
degburnedssu$Season <- "Summer"

degunburnedssu <- data.frame(colSums(dframeDDsusu[,4:69]))
colnames(degunburnedssu) <- "Degree"
degunburnedssu$Treatment <- "Unburned"
degunburnedssu$Season <- "Summer"

#winter
degburnedsw <- data.frame(colSums(dframeDDsbw[,4:69]))
colnames(degburnedsw) <- "Degree"
degburnedsw$Treatment <- "Burned"
degburnedsw$Season <- "Winter"

degunburnedsw <- data.frame(colSums(dframeDDsuw[,4:69]))
colnames(degunburnedsw) <- "Degree"
degunburnedsw$Treatment <- "Unburned"
degunburnedsw$Season <- "Winter"


degplantss <- rbind(degburnedsa,degunburnedsa, degburnedssp,degunburnedssp, degburnedssu,degunburnedssu, degburnedsw,degunburnedsw)


# now start testing
# autumn
hist(degburnedsa$Degree)
hist(degunburnedsa$Degree)

ks.test(degburnedsa$Degree,degunburnedsa$Degree,alternative="gr")

# spring
hist(degburnedssp$Degree)
hist(degunburnedssp$Degree)

ks.test(degburnedssp$Degree,degunburnedssp$Degree,alternative="gr")

# summer
hist(degburnedssu$Degree)
hist(degunburnedssu$Degree)

ks.test(degburnedssu$Degree,degunburnedssu$Degree,alternative="gr")

# winter
hist(degburnedsw$Degree)
hist(degunburnedsw$Degree)

ks.test(degburnedsw$Degree,degunburnedsw$Degree,alternative="gr")



# we need to calculate group means
muplantss <- ddply(degplantss, c("Treatment","Season"), summarise, grp.mean=mean(Degree), grp.se=se(Degree))
print(muplantss)


### median degree of 18 networks
# finally we want to follow the same steps to obtain a median degree for each network to take forward for analysis with other network metrics

# moths
# select columns in specific order: Species, Treatment, Pollen, Season, Year
dframe1ddsp <- dframe1[,c(2,7,10,12:78)]


# summarise pollen transport for each insect species within each treatment
# this aggregates all pollen grains carried by an insect species for each plant species within each treatment (across all samples)

dframeDDsp <- ddply(dframe1ddsp, .(Family_Species,Season,Treatment,Year), numcolwise(sum))

# we don't need the quantitative information here

dframeDDsp[,5:70][dframeDDsp[,5:70] > 0] <- 1

# now we calculate the degree for each insect species in each row

dframeDDsp$Degree <- rowSums(dframeDDsp[,5:70])


# now we want to calculate averages for each of the 18 networks
mumothnetworks <- ddply(dframeDDsp, c("Treatment","Season","Year"), summarise, grp.mean=mean(Degree), grp.median=median(Degree),grp.se=se(Degree))
print(mumothnetworks)
mumothnetworks$Level <- factor("Moths")


# plants

# the necessary information is already there in dframeDDsp cols 5:72, but it'll take some jiggery-pokery to extract it

# we need a variable to allow the networks to be split:

dframeDDsp$ExactSeason <- do.call(paste, c(dframeDDsp[c("Treatment","Season","Year")], sep = "_"))
dframeDDsp$ExactSeason <- as.factor(dframeDDsp$ExactSeason)
summary(dframeDDsp$ExactSeason)

dframesDD <- split(dframeDDsp, list(dframeDDsp$ExactSeason))

plantdegree <- function(x) {
  degree <- data.frame(colSums(x[,5:70]))
  colnames(degree) <- "Degree"
  degree$Treatment <- x[1,3]
  degree$Season <- x[1,2]
  degree$Year <- x[1,4]
  return(degree)
}

degree <- lapply(dframesDD, plantdegree)

degree.merge <- do.call("rbind", degree)

# now we want to calculate averages for each of the 18 networks
muplantsnetworks <- ddply(degree.merge, c("Treatment","Season","Year"), summarise, grp.mean=mean(Degree), grp.median=median(Degree),grp.se=se(Degree))
print(muplantsnetworks)
muplantsnetworks$Level <- factor("Plants")

# merge the two sets together
munetworks <- rbind(mumothnetworks,muplantsnetworks)

summary(munetworks)

# and write them to a txt file for downstream analysis
write.table(munetworks, "Data\\DegreeBySample.txt", sep="\t", row.names=FALSE)



