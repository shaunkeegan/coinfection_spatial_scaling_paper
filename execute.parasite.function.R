##  Generates dataframes with columns for neighbourhood nematode and Eimeria variables 
##      at each specified neighbourhood size corresponding to each focal capture
##


## CLEAR MEMORY
rm(list = ls())



## LOAD DATA 
data <- read.csv("data input/raw data.csv", header = T)



## LOAD FUNCTIONS AND LIBRARIES
source("exppoints.R")
source("prop.on.grid.R")
source("make.dist.vector.R")
source("parasite.function.vary.time.R")

library(gtools) #for permutation, used to construct location combinations for grid coordinates



## SET NUMERIC CONSTANTS 
min.point <- 0    # Minimum value of lat/long on a given grid
max.point <- 50   # Maximum value of lat/long on a given grid
trap.interval <- 10    #Distance between adjacent traps



# DISTANCES FOR NEIGHBOURHOOD SIZES TO EXAMINE
dist <- make.dist.vector(max.point)
dist <- sort(dist)
dist <- dist[c(2:10,14)]
dist_trunc <- round(dist, digits = 1)



## SELECT DATA VARIABLES FOR neighbour AND focal PARASITE COMBINATION TO PROCESS
focal.variables <- c("EhungInt", "EhungINF") 
focal.NH <- "EhungINF"  #the focal Eimeria in the neighbourhood [needs to match same sp as focal.variables]
#focal.variables <- c("EappInt", "EappINF") 
#focal.NH <- "EappINF"  #the focal Eimeria in the neighbourhood [needs to match same sp as focal.variables]

neighbour.variable.INF <- "HpolINF"     #Infection status variable for NH parasite
neighbour.variable.EPG <- "H.poly.EPG"  #Infection burden (EPG) variable for NH parasite



## SELECT neighbour AND focal DATA SUBSETS 
neighbour <- data
focal <- neighbour[which (neighbour$treated.Y.N == 0),] #Focals should be untreated

keep.foc <- c("ID", "ID.CapDate", "X", "Y", "Capture.date","Grid", "year", "Sex", "Age", focal.variables, neighbour.variable.INF, neighbour.variable.EPG)
keep.nei <- c("ID", "ID.CapDate", "X", "Y", "Capture.date","Grid", "year", neighbour.variable.INF, neighbour.variable.EPG, focal.NH) 
focal <- na.omit(focal[,keep.foc])
dim(focal)	
neighbour <- na.omit(neighbour[,keep.nei])
dim(neighbour)	



## CREATE CUMULATIVE NEIGHBOURHOOD VARIABLES
# Proportion of grid 
PropOnGrid <- rep("PropOnGrid", length(dist))
PropOnGrid <- paste(PropOnGrid, dist_trunc, sep = "_")

# Neighbourhood nematode infection counts
NemNoInf <- rep("NemNoInf", length(dist))
NemNoInf <- paste(NemNoInf, dist_trunc, sep = "_")
NemNoUninf <- rep("NemNoUninf", length(dist))
NemNoUninf <- paste(NemNoUninf, dist_trunc, sep = "_")
NemPrev <- rep("NemPrev", length(dist))
NemPrev <- paste(NemPrev, dist_trunc, sep = "_")

NemNoTotal <- rep("NemNoTotal", length(dist))
NemNoTotal <- paste(NemNoTotal, dist_trunc, sep = "_")

NemTotBurden <- rep("NemTotBurden", length(dist))
NemTotBurden <- paste(NemTotBurden, dist_trunc, sep = "_")
NemMeanBurden <- rep("NemMeanBurden", length(dist))
NemMeanBurden <- paste(NemMeanBurden, dist_trunc, sep = "_")

#Neighbourhood Eimeria infection counts
EimNoInf <- rep("EimNoInf", length(dist))
EimNoInf <- paste(EimNoInf, dist_trunc, sep = "_")
EimNoUninf <- rep("EimNoUninf", length(dist))
EimNoUninf <- paste(EimNoUninf, dist_trunc, sep = "_")
EimPrev <- rep("EimPrev", length(dist))
EimPrev <- paste(EimPrev, dist_trunc, sep = "_")
EimNoTotal <- rep("EimNoTotal", length(dist))   
EimNoTotal <- paste(EimNoTotal, dist_trunc, sep = "_")



## GENERATE NEIGHBOURHOOD DATA FOR ANALYSIS
daysbefore <- 17
days17 <- parasite.function.vary.time(neighbour, focal, neighbour.variable.INF, neighbour.variable.EPG, focal.NH, dist, daysbefore)
daysbefore <- 34
days34 <- parasite.function.vary.time(neighbour, focal, neighbour.variable.INF, neighbour.variable.EPG, focal.NH, dist, daysbefore)
daysbefore <- 51
days51 <- parasite.function.vary.time(neighbour, focal, neighbour.variable.INF, neighbour.variable.EPG, focal.NH, dist, daysbefore)
daysbefore <- 365
daysall <- parasite.function.vary.time(neighbour, focal, neighbour.variable.INF, neighbour.variable.EPG, focal.NH, dist, daysbefore)



## OUTPUT TO csv FILES FOR STORAGE AND SUBSEQUENT PROCESSING
write.csv(days17, "data output/DATA - Ehung_17.csv") 
write.csv(days34, "data output/DATA - Ehung_34.csv")
write.csv(days51, "data output/DATA - Ehung_51.csv")
write.csv(daysall, "data output/DATA - Ehung_all.csv")

#write.csv(days17, "data output/DATA - Eapp_17.csv")
#write.csv(days34, "data output/DATA - Eapp_34.csv")
#write.csv(days51, "data output/DATA - Eapp_51.csv")
#write.csv(daysall, "data output/DATA - Eapp_all.csv")



