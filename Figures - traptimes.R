########################################
### Plot trap start and finish times ###
########################################

### Clear the workspace
rm(list=ls())

# set the working directory as the script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### install if necessary and then load the libraries you need

j <- c("ggplot2","ggmap","plyr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### load up Callum's custom set of functions
f <- c("MultiplotFunction.R")
lapply(f, source)


### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.csv("Data/TrapTimes.csv", header=TRUE)

summary(dframe1)

dframe1$Month <- ifelse(grepl("-11-",dframe1$Date),"Nov",
                        ifelse(grepl("-12-",dframe1$Date),"Dec",
                               ifelse(grepl("-01-",dframe1$Date),"Jan",
                                      ifelse(grepl("-02-",dframe1$Date),"Feb",
                                             ifelse(grepl("-03-",dframe1$Date),"Mar",
                                                    ifelse(grepl("-04-",dframe1$Date),"Apr",
                                                           ifelse(grepl("-05-",dframe1$Date),"May",
                                                                  ifelse(grepl("-06-",dframe1$Date),"Jun",
                                                                         ifelse(grepl("-07-",dframe1$Date),"Jul",
                                                                                ifelse(grepl("-08-",dframe1$Date),"Aug",
                                                                                       ifelse(grepl("-09-",dframe1$Date),"Sep",
                                                                                              ifelse(grepl("-10-",dframe1$Date),"Oct",
                                                                                                     "Fail"))))))))))))
dframe1$Month<-ordered(dframe1$Month, levels=c("Dec","Nov","Oct","Sep","Aug","Jul","Jun","May","Apr","Mar","Feb","Jan"))

summary(dframe1)



means <- ddply(dframe1, .(Month,Job), numcolwise(mean))
max <- ddply(dframe1, .(Month,Job), numcolwise(max))
min <- ddply(dframe1, .(Month,Job), numcolwise(min))

colnames(means) <- c("Month","Job","Mean")
colnames(max) <- c("Month","Job","Max")
colnames(min) <- c("Month","Job","Min")



dframe1.5 <- merge(means,max)
dframe2 <- merge (dframe1.5,min)


scaleFUN <- function(x) sprintf("%.2f", x)


g1 <- ggplot(dframe2, aes(x=Mean,y=Month,group=Job))+
  geom_point(size=1.5,colour="black",fill="black",stroke=2,alpha=1,
             aes(shape=Job))+
  scale_shape_manual(values=c(25,24),
                     labels=c("Collection","Set-up"))+
  geom_errorbarh(xmin=dframe2$Min,xmax=dframe2$Max,
                 colour="black",size=0.8, height=0.5)+
  scale_x_continuous(limits=c(0,24),breaks = c(0,4,8,12,16,20,24),labels=scaleFUN)+
  xlab("Time") + ylab("Month")+
  theme(panel.background=element_rect(fill="white"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major.x=element_line(colour="gray70"),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=15),
        axis.text=element_text(color="black"),
        legend.title=element_blank())


g1
