############################################
####   Script for interaction turnover  ####
############################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("plyr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### load up Callum's custom set of functions
f <- c("NetworkFunction.R")
lapply(f, source)



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoctThreshold.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

### prepare the data for network analysis

# add a column to identify which network each row belongs to

dframe1$Network <- do.call(paste, c(dframe1[c("Season","Year","Treatment")], sep = "_"))
dframe1$Network <- as.factor(dframe1$Network)
summary(dframe1$Network)

# and another for which network pair each row belongs to

dframe1$Pair <- do.call(paste, c(dframe1[c("Season","Year")], sep = "_"))
dframe1$Pair <- as.factor(dframe1$Pair)
summary(dframe1$Pair)

# trim off extra columns and reorder
dframe1r <- dframe1[,c(2,7,(length(dframe1)-1),length(dframe1),13:(length(dframe1)-2))]

# change each interaction to a 1 for flower-visitation
dframe1r[,5:length(dframe1r)][dframe1r[,5:length(dframe1r)] > 0] <- 1

# summarise plant-insect interactions for each insect species within each sample
# this produces semi-quantitative data - e.g. if two individuals of an insect species
# interact with a plant species in a sample,
# the entry for that interaction within that sample will be 2

dframe2 <- ddply(dframe1r, .(Family_Species,Treatment,Pair,Network), numcolwise(sum))

# for later, we'll need a list of generic column names
colnames.gen <- colnames(dframe2)

# and a list of generic row names (= a list of insect species)
rownames.gen <- data.frame(levels(dframe2$Family_Species))


# create separate dframes depending on treatment

dframe.fire <- dframe2[ which(dframe2$Treatment=="Fire"), ]
dframe.nofire <- dframe2[ which(dframe2$Treatment=="NoFire"), ]


# split each sample into a separate dataframe
fires <- split(dframe.fire, list(dframe.fire$Pair))  # this creates a list of smaller dframes, one for each level of sample
summary(fires)

nofires <- split(dframe.nofire, list(dframe.nofire$Pair))
summary(nofires)

# we now have a list of networks for each treatment - the lists are in the same order (alphabetical, by season)
# therefore the pairs are fires[1] and nofires[1], to fires[n] and nofires[n]

output <- c("i","noposs","a","b","c","s0","sfp","sf","sp","Bfp","Bf2","Bp","B02","JBD")


# now run the loop - this is a lot of lines!
# scroll down until you find ###### end of loop ######, select from here to there, and run it all


for (i in 1:length(fires)){
  
  
  ## SECTION 1 - LOAD NETWORKS A AND B (CHANGING FILENAMES AS APPROPRIATE) AND 
  ## CHECK MATCH AND CONVERT TO ZEROS AND ONES
  
  
  
  # Read in files
  
  # instead, let's try reading in networks from a list of networks prepared in script 10
  netA <- data.frame(fires[i])
  netB <- data.frame(nofires[i])
  
  
  # set generic colnames
  colnames(netA) <- colnames.gen
  colnames(netB) <- colnames.gen
  
  # set generic rownames
  netA <- merge(rownames.gen, netA, by.x = 1, by.y = 1, all.x = TRUE)
  netB <- merge(rownames.gen, netB, by.x = 1, by.y = 1, all.x = TRUE)
  
  
  # set column 1 as rownames
  rownames(netA) <- netA[,1]
  rownames(netB) <- netB[,1]
  
  # throw out extra columns
  netA <- netA[,c(5:length(netA))]
  netB <- netB[,c(5:length(netB))]
  
  # set all NAs to zeroes
  netA[is.na(netA)] <- 0
  netB[is.na(netB)] <- 0
  
  
  # transpose the networks
  netA <- data.frame(t(netA))
  netB <- data.frame(t(netB))
  
  # add in new first column
  netA <- cbind(rownames(netA),netA)
  netB <- cbind(rownames(netB),netB)
  
  # rename the extra columns
  colnames(netA)[1] <- "Species"
  colnames(netB)[1] <- "Species"
  
  # extract the row and column names
  netArownames<-rownames(netA)
  netBrownames<-rownames(netB)
  netAcolnames<-colnames(netA[2:length(netA)])
  netBcolnames<-colnames(netB[2:length(netB)])
  
  # set up TRUE / FALSE matrix to check that matrix names match
  checkcolnames<-netAcolnames==netBcolnames
  checkrownames<-netArownames==netBrownames
  
  which(checkcolnames=="FALSE") # checks if columns and rows match
  which(checkrownames=="FALSE") # if output is integer (0) then matrix names match
  # if output is a number this identifies places in the list that do not match
  # e.g. an output of 3 and 7 in columns would suggest that matrix A and B 
  # differ in columns 3 and  7 (this could for example mean these two species have been 
  # swapped around in your data sheets.
  
  # Remove the plant name text column from the matrix leaving only the numbers
  onetl<-length(netA) # original column number before creating matrix
  netA<-as.matrix(netA[,2:onetl]) # removes first column so only numbers
  netB<-as.matrix(netB[,2:onetl]) # remain in the matrices
  
  # convert data to zero or one where necessary
  netA[netA > 0] <-1 # ensure matrix consists of 0s and 1s
  netB[netB > 0] <-1
  
  
  
  ## SECTION 2 - DETERMINE THE PRESENCE / ABSENCE OF FLOWERS AND POLLINATORS IN 
  ## (I.E. DETERMINE IF F OR P IS INVOLVED IN ANY INTERACTION WITH ANY OTHER F OR P)
  ## BY SUMMING THE COLUMNS AND ROWS
  
  netc<-length(netA[1,]) # count number of columns
  netr<-length(netA[,1]) # count number of rows
  netc
  netr
  
  # create coded vector to determine which pollinators are present in matrix A
  PA<-numeric(netc)
  for (n in 1:netc){
    if(sum(netA[,n])>0){PA[n]<-1}}
  PA
  
  # create coded vector to determine which flowers are present in matrix A
  FA<-numeric(netr)
  for (n in 1:netr){
    if(sum(netA[n,])>0){FA[n]<-1}}
  FA
  
  # create coded vector to determine which pollinators are present in matrix B
  PB<-numeric(netc)
  for (n in 1:netc){
    if(sum(netB[,n])>0){PB[n]<-1}}
  PB
  
  # create coded vector to determine which flowers are present in matrix B
  FB<-numeric(netr)
  for (n in 1:netr){
    if(sum(netB[n,])>0){FB[n]<-1}}
  FB
  
  
  ## SECTION 3 - CHANGE THE MATRIX INTO A LONG FORMAT TABLE FOR PROCESSING
  ## (WITH ROW AND COLUMN NAMES LAID OUT IN CORRECT ORDER)
  
  
  # turn matrix into a simple vector
  cnetA<-c(netA)
  cnetA 
  
  # expand the data on species presence / absence to span the new vector
  expandpresabA<-expand.grid(FA,PA)
  # create a list of the herbivore and plant names in order
  expandnamesA<-expand.grid(netArownames, netAcolnames)
  expandnamesA
  # List the names of pollinators and flowers in the correct order alongside 
  # data for interactions and presence / absence of flowers and pollinators
  datafA<-data.frame("flowersA"=expandnamesA[,1],"pollinatorsA"=expandnamesA[,2],
                     "flowerpres"=expandpresabA[,1],"pollpres"=expandpresabA[,2], cnetA)
  
  
  # visual check
  head(datafA)
  
  
  # repeat for network B
  # turn matrix into a simple vector
  cnetB<-c(netB)
  cnetB 
  
  # expand the data on species presence / absence to span the new vector
  expandpresabB<-expand.grid(FB,PB)
  # create a list of the herbivore and plant names in order
  expandnamesB<-expand.grid(netBrownames, netBcolnames)
  expandnamesB
  # List the names of pollinators and flowers in the correct order alongside 
  # data for interactions and presence / absence of flowers and pollinators
  datafB<-data.frame("flowersB"=expandnamesB[,1],"pollinatorsB"=expandnamesB[,2],
                     "flowerpres"=expandpresabB[,1],"pollpres"=expandpresabB[,2], cnetB)
  
  
  # visual check
  head(datafB)
  
  
  
  ### SECTION 4 - TEST FOR AND REMOVE CASES WHERE THE FLOWER AND POLLINATOR EXIST IN AT 
  ### LEAST ONE NETWORK BUT NO INTERACTION EXISTS IN EITHER NETWORK,
  ### AND TEST FOR AND REMOVE CASES WHERE INTERACTIONS ARE IMPOSSIBLE IN BOTH NETWORKS
  combineddat<-data.frame(datafA,datafB)
  for (n in 1:length(datafA[,1]))
  {combineddat[n,11]<-sum(combineddat[n,3],combineddat[n,4])
  combineddat[n,12]<-sum(combineddat[n,8],combineddat[n,9])
  combineddat[n,13]<-sum(combineddat[n,3],combineddat[n,8])
  combineddat[n,14]<-sum(combineddat[n,4],combineddat[n,9])
  ifelse ((combineddat[n,11]==2 | combineddat[n,12]==2) & (combineddat[n,5]==0 & combineddat[n,10]==0), 
          combineddat[n,15]<-1, combineddat[n,15]<-0)}
  networknulls<-which(combineddat[,15]==1) # identify which rows contain these  
  networknulls                             # 'nointeraction' conditions
  
  
  
  ### remove cases where the flower is absent from both networks, and do the same for the insect
  flowernulls <- which(combineddat[,13]==0)
  insectnulls <- which(combineddat[,14]==0)
  
  # combine the lists of rows for removal
  allnulls <- c(networknulls,flowernulls,insectnulls)
  
  
  # remove all those rows
  if(length(allnulls!=0))
  {
    datafA<-datafA[-allnulls,]
    datafB<-datafB[-allnulls,]
  }
  
  
  
  ## SECTION 5 - HERE THE VECTOR - newcnetA or B - IS CODED TO 
  ## DETERMINE THE POTENTIAL FOR INTERACTION AND / OR THE PRESENCE OF 
  ## POLLINATORS AND FLOWERSS
  
  # fp indicates presence of both flower and pollinator, f indicates only the flower is
  # present, p only the pollinator, ni indicates neither are present.
  
  # for network A
  
  # create the new vector newcnetA
  newcnetA<-numeric(length(datafA[,1])) 
  
  for(n in 1:length(datafA[,1]))
  {
    ifelse (datafA[n,5]==1, newcnetA[n]<-"fp",                                                      # if interaction exists, assumed both present
            ifelse(datafA[n,5]==0 & datafA[n,3]==1 & datafA[n,4]==1, newcnetA[n]<-"fp",                   # if no interaction but both present
                   ifelse (datafA[n,5]==0 & datafA[n,3]==1 & datafA[n,4]==0, newcnetA[n]<-"f",                 # flower only present 
                           ifelse(datafA[n,5]==0 & datafA[n,3]==0 & datafA[n,4]==1, newcnetA[n]<-"p",                # pollinator only present
                                  newcnetA[n]<-"ni"))))                                                                     # if all above fail, assumed neither present
  }
  
  head(newcnetA)
  
  # same process for network B
  
  # create the new vector newcnetB
  newcnetB<-numeric(length(datafB[,1])) 
  
  for(n in 1:length(datafB[,1]))
  {
    ifelse (datafB[n,5]==1, newcnetB[n]<-"fp",
            ifelse(datafB[n,5]==0 & datafB[n,3]==1 & datafB[n,4]==1, newcnetB[n]<-"fp",
                   ifelse (datafB[n,5]==0 & datafB[n,3]==1 & datafB[n,4]==0, newcnetB[n]<-"f",  
                           ifelse(datafB[n,5]==0 & datafB[n,3]==0 & datafB[n,4]==1, newcnetB[n]<-"p",
                                  newcnetB[n]<-"ni"))))
  }
  head(newcnetB)
  
  
  
  ## SECTION 6 - CREATE THE FINAL VECTOR WHICH CODES ALL THE POTENTIAL INTERACTIONS
  ## ETC. PARTITIONED ACCORDING TO (NOVOTNY ET AL. 2009 - INSECT CONSERVATION AND
  ## DIVERSITY - WHERE "I" EQUATES TO A POTENTIAL INTERACTION (I.E. FLOWER AND POLLINATOR ARE
  ## PRESENT IN BOTH NETWORKS) - THIS INCLUDES THE SHADED AND THE ZERO VALUE IN NOVOTNY ET AL.,
  ## TO DETERMINE THE SHADED (I.E. TRUE SHARED INTERACTIONS) AN ADDITIONAL STEP IS TAKEN WHERE THE 
  ## ORIGINAL NETWORKS ARE SUMMED AND ALL VALUES OF 2 INDICATE AN INTERACTION IN BOTH NETWORKS AND EQUATES 
  ## TO THE a VALUE IN THE FINAL CALCULATIONS, THE POTENTIAL BUT UNFULFILLED INTERACTIONS (ZERO IN
  ## NOVOTNY'S ARTICLE), ARE THEN CALCULATED AS THE COUNT OF I's-a. NB. THE CELLS LEFT BLANK IN NOVOTNY ET AL.
  ## ARE CODED ZERO HERE.
  
  finnet<-numeric(length(newcnetA)) # create a vector of appropriate length
  
  # code the final table (i.e. figure 1 table C in Novotny et al)
  for (n in 1:length(newcnetA))
  {
    ifelse(newcnetA[n]=="fp" & newcnetB[n]=="ni", finnet[n]<-"fp",                      # if both missing from one network and present in other
           ifelse(newcnetA[n]=="fp" & newcnetB[n]=="f", finnet[n]<-"p",                      # if insect missing from one and both present in other
                  ifelse(newcnetA[n]=="fp" & newcnetB[n]=="p", finnet[n]<-"f",                    # if flower missing from one and both present in other
                         ifelse(newcnetA[n]=="fp" & newcnetB[n]=="fp", finnet[n]<-"I",                 # if both present in both but no interaction
                                ifelse(newcnetA[n]=="f" & newcnetB[n]=="p", finnet[n]<-"0",                 # if insect missing from one and flower from other
                                       ifelse(newcnetA[n]=="f" & newcnetB[n]=="fp", finnet[n]<-"p",         # if insect missing from one and both present in other
                                              ifelse(newcnetA[n]=="p" & newcnetB[n]=="f", finnet[n]<-"0",             # if insect missing from one and flower from other
                                                     ifelse(newcnetA[n]=="p" & newcnetB[n]=="fp", finnet[n]<-"f",          # if flower missing from one and both present in other
                                                            ifelse(newcnetA[n]=="ni" & newcnetB[n]=="fp",finnet[n]<-"fp",       # if both missing from one network and present in other
                                                                   finnet[n]<-"error")))))))))                                  # all other pairs (i.e. ni + ni)
  }
  
  # visual check of finnet
  finnet
  
  # check for errors (if length is zero then all is well)
  length(finnet[finnet=="error"])
  
  noposs<-length(finnet[finnet=="0"]) # i.e. where one network has only the pollinator and the other
  # has only the flower - thus no interaction is possible - blanks
  # in Novotny's table
  paste ("No. of cases where no interaction is possible ", noposs)
  
  
  # sum the two networks - where 2s occur this is an interaction found in both networks
  
  sumnet<-netA+netB
  sumnet
  a<-length(sumnet[sumnet==2])
  paste ("No. of interactions in both networks (a) ", a) # number of interactions present in both grids 
  
  b<- length(netA[netA=="1"])-a
  paste ("No. of interactions in network A only (b) ", b) 
  
  c<-length(netB[netB=="1"])-a
  paste ("No. of interactions in network B only (c) ", c)
  
  
  paste("sum of (a+b+c) ", (a+b+c))
  
  s0<-length(finnet[finnet=="I"])-a  # b0 + c0 i.e. where both flower and pollinator exist in both 
  s0                                 # networks but interaction is only in one
  
  sfp<-length(finnet[finnet=="fp"]) # bph + cph i.e. where both flower and pollinator are absent 
  sfp                               # from one of the networks
  
  sf<-length(finnet[finnet=="f"]) # bp + cp i.e. where the flower is absent 
  sf                               # from one of the networks
  
  
  sp<-length(finnet[finnet=="p"]) # bh + ch i.e. where the pollinator is absent 
  sp                              # from one of the networks
  
  
  #calculate Jacard Beta Diversity
  Bfp<-sfp/(a+b+c) # 
  paste ("Beta fp ", Bfp) # change in both flower and pollinator
  
  Bf<-sf/(a+b+c)
  paste ("Beta f ", Bf)  # change in flower only
  
  Bp<-sp/(a+b+c)
  paste ("Beta p ", Bp)  # change in pollinator only
  
  B0<-s0/(a+b+c)
  paste ("Beta 0 ", B0)  # all necessary species present - change in interaction 
  
  
  JBD<-Bfp+Bf+Bp+B0
  paste ("Jaccard Beta Diversity ", JBD)
  
  result <- c(i,noposs,a,b,c,s0,sfp,sf,sp,Bfp,Bf,Bp,B0,JBD)
  
  output <- rbind(output,result)
  
}                                               

######  end of loop ######

# now tidy up the output table and add in information about which season is which

networks <- c("Season",names(fires))
final <- data.frame(cbind(networks,output))

rownames(final) <- final[,1]
colnames(final) <- as.character(unlist(final[1,]))

final <- final[-1,]


# finally, we want to allow the final dataframe to be sorted chronologically, not alphabetically

order <- c(1:9)
Season <- c("Spring_1","Summer_1","Autumn_2","Winter_2","Spring_2","Summer_2","Autumn_3","Winter_3","Spring_3")

sort <- data.frame(cbind(order,Season))
final <- merge(sort,final)


# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(final, "Results\\BetaDiversityTurnoverThreshold.txt", sep="\t", row.names=FALSE)







########################### development #######################
