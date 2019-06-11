### Plants ####
# run this script immediately after script 2, with the same set of libraries and functions loaded

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe4<-read.table("Data/PlantTransects.txt", header=TRUE)
dframe4 <- read.csv(file.choose()) # PlantTransectsMelt ?

dframe4$Year <- factor(dframe4$Year)
dframe4$Treatment <- relevel(dframe4$Treatment, "NoFire")
summary(dframe4)


# We don't need the individual TransectID for this analysis, so get rid of it:

dframe4r <- dframe4[,c(3:15)]
summary(dframe4r)

# each row currently contains 1 plant species, so add a column for "Count"; 
# this is 1 in every instance except the transect with zero insects
dframe4r$SpeciesRichness <- ifelse(dframe4r$PlantSpecies=="none",0,1)

# summarise the dataframe to how many individuals of each species in each sample
dframe5 <- ddply(dframe4r, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year), numcolwise(mean))

dframe5SR <- ddply(dframe4r, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year), numcolwise(sum))


dframe5$Month <- ordered(dframe5$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

summary(dframe5)


# Floral abundance
plot(PlantCoverage ~ Treatment, dframe5)
plot(PlantCoverage ~ Season, dframe5)
plot(PlantCoverage ~ Month, dframe5)
plot(PlantCoverage ~ Year, dframe5)
plot(PlantCoverage ~ interaction(Treatment,Season), dframe5)

hist(dframe5$PlantCoverage)
hist(log(dframe5$PlantCoverage+1,10))
summary(dframe5$PlantCoverage)

# possibly overdispersion, possibly zero-inflation (1st quartile is 0)
# because it's average data it's non-integer so Poisson won't work
# but the log looks close to Gaussian so let's try that first
mean(dframe5$PlantCoverage)
var(dframe5$PlantCoverage)  # var > mean therefore data are overdispersed

# Gaussian

dframe5$tPlantCoverage <- log(dframe5$PlantCoverage + 1,10)
summary(dframe5$tPlantCoverage)
hist(dframe5$tPlantCoverage)


model5G <- lmer(tPlantCoverage ~ Treatment*Season
                + (1|Year) + (1|Site) + (1|Date),
                data = dframe5)

summary(model5G)
drop1(model5G, test = "Chi")

chkres(model5G,dframe5$Treatment,dframe5$Season) # these look ok-ish 
# but have a fairly clear negative trend and problems with the fit of season


# Poisson
# can't deal with non-integers so we'll have to use total rather than mean coverage

model5P <- glmer(PlantCoverage ~ Treatment*Season
                 + (1|Year) + (1|Site) + (1|Date),
                 family = poisson (link = "log"),
                 data = dframe5SR)

summary(model5P)
drop1(model5P, test = "Chi")

chkres(model5P,dframe5SR$Treatment,dframe5SR$Season) # these are pretty bad

# nbinom

model5NB <- glmer.nb(PlantCoverage ~ Treatment*Season
                     + (1|Year) + (1|Site) + (1|Date),
                     data = dframe5SR)

chkconv(model5NB)

summary(model5NB)
drop1(model5NB, test = "Chi")

chkres(model5NB,dframe5SR$Treatment,dframe5SR$Season) # these are pretty bad as well

# we need to bite the bullet and use zero-inflated models
# zero-inflated poisson

#model5ZIP <- glmmadmb(PlantCoverage ~  Treatment*Season
#                      + (1|Year) + (1|Site) + (1|Date), #Random effects
#                      zeroInflation=TRUE,
#                      family = "poisson",
#                      link = "log",
#                      data = dframe5SR)

#summary(model5ZIP)

#Anova(model5ZIP, type="III")
#drop1(model5ZIP, test="Chisq")

#chkres.zi(model5ZIP,dframe5SR$Treatment,dframe5SR$Season) # these are pretty horrible


# zero-inflated nbinom

#model5ZINB <- glmmadmb(PlantCoverage ~  Treatment*Season
#                       + (1|Year) + (1|Site) + (1|Date), #Random effects
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

summary(model5G)
drop1(model5G, test = "Chi")




# Floral Species Richness
plot(SpeciesRichness ~ Treatment, dframe5SR)
plot(SpeciesRichness ~ Season, dframe5SR)
plot(SpeciesRichness ~ Month, dframe5SR)
plot(SpeciesRichness ~ interaction(Treatment,Season), dframe5SR)

hist(dframe5SR$SpeciesRichness)

# possibly overdispersion (but not sure...)
mean(dframe5SR$SpeciesRichness)
var(dframe5SR$SpeciesRichness)  # var > mean therefore data are overdispersed


# Poisson model
model6P <- glmer(SpeciesRichness ~ Treatment*Season
                 + (1|Year) + (1|Site) + (1|Date),
                 family = poisson (link="log"),
                 data = dframe5SR)

summary(model6P)
drop1(model6P, test = "Chi")

chkres(model6P,dframe5SR$Treatment,dframe5SR$Season) # these are BAAAD - they look like the data needs a log transformation


# negative binomial
model6NB <- glmer.nb(SpeciesRichness ~ Treatment*Season
                     + (1|Year) + (1|Site) + (1|Date),
                     data = dframe5SR)

summary(model6NB)
drop1(model6NB, test="Chi")

chkres(model6NB,dframe5SR$Treatment,dframe5SR$Season)


# try a log-transformed gaussian as that worked previously
dframe5SR$tSpeciesRichness <- log(dframe5SR$SpeciesRichness+1,10)
hist(dframe5SR$tSpeciesRichness)

model6G <- lmer(tSpeciesRichness ~ Treatment*Season
                + (1|Year) + (1|Site) + (1|Date),
                data = dframe5SR)

summary(model6G)
drop1(model6G, test = "Chi")

chkres(model6G,dframe5SR$Treatment,dframe5SR$Season) # these are not terrible, though the zeroes in summer continue to have a visible effect

