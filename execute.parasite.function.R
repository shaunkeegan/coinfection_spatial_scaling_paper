## --------------------- EXECUTE: parasite.function
## Reads in raw data file
## Calculates unique neighbourhood sizes
## Calls 'parasite.function.vary.time' to create dataset of nematode prevalences around each focal at each specified neighbourhood size
##





## ----------------------------------------------------------- LOAD DATA 

setwd ("~/Documents/Projects/Coinfection Spatial Scaling/
       Project Files/2_neighbourhood_calculation")

data <- read.csv("input/whole_dataset_woodmouse.csv", header = T)


## ------------------------------------------------------ LOAD FUNCTIONS


source("funcs/make.dist.vector.R")
source("funcs/parasite.function.vary.time.R")
source("funcs/exppoints.R")
source("funcs/prop.on.grid.R")
library(gtools) #for 'permutations' function to create list of all possible neighbourhood sizes


## ----------------------------------------------------- INPUT VARIABLES

# ---------- Numeric Variables 

min.point <- 0    # Minimum value of lat/long
max.point <- 50  # Maximum value of lat/long on a given grid
trap.interval <- 10


# ---------- Create distance vector of all possible neighbourhood sizes

dist <- make.dist.vector(max.point, trap.interval)
dist <- sort(dist)
dist <- dist[c(2:10,14)]
dist_trunc <- round(dist, digits = 1)


# ---------- Name Variables 

PropOnGrid <- rep("PropOnGrid", length(dist))
PropOnGrid <- paste(PropOnGrid, dist_trunc, sep = "_")
NoInf <- rep("NoInf", length(dist))
NoInf <- paste(NoInf, dist_trunc, sep = "_")
NoUninf <- rep("NoUninf", length(dist))
NoUninf <- paste(NoUninf, dist_trunc, sep = "_")
Prev <- rep("Prev", length(dist))
Prev <- paste(Prev, dist_trunc, sep = "_")
NoTotal <- rep("NoTotal", length(dist))
NoTotal <- paste(NoTotal, dist_trunc, sep = "_")


## ------------------------------------------------ SUBSET DATA


focal.variable <- "EhungInt" #Specify Ehung or Eapp as needed
neighbour.variable <- "HpolINF"


neighbour <- data
focal <- neighbour[which (neighbour$treated.Y.N == 0),]


keep1 <- c("ID", "ID.CapDate", "X", "Y", "Capture.date","Grid", "year", "Sex", "Age", focal.variable, neighbour.variable)
keep2 <- c("ID", "ID.CapDate", "X", "Y", "Capture.date","Grid", "year", neighbour.variable)
focal <- na.omit(focal[,keep1])
dim(focal)	
neighbour <- na.omit(neighbour[,keep2])
dim(neighbour)	

## ------------------------------------------------ GENERATE DATA FOR ANALYSIS


daysbefore <- 17
days17 <- parasite.function.vary.time(neighbour, focal, neighbour.variable, dist, daysbefore)
daysbefore <- 34
days34 <- parasite.function.vary.time(neighbour, focal, neighbour.variable, dist, daysbefore)
daysbefore <- 51
days51 <- parasite.function.vary.time(neighbour, focal, neighbour.variable, dist, daysbefore)
daysbefore <- 365
daysall <- parasite.function.vary.time(neighbour, focal, neighbour.variable, dist, daysbefore)

#Output data files: change names to Ehung or Eapp as needed
write.csv(days17, "input/Ehung_hpol17.csv")
write.csv(days34, "input/Ehung_hpol34.csv")
write.csv(days51, "input/Ehung_hpol51.csv")
write.csv(daysall, "input/Ehung_hpolall.csv")



