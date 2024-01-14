##  Creates vector of all possible neighbourhood sizes 
##  from 0 to specified 'max.point' on a trapping grid with trap locations 'trap.interval'. apart
##

##----------------------------------------------------------------------------##
##                              ARGUMENTS                                     ##
##----------------------------------------------------------------------------##
##
## maxdist: (numeric) - maximum specified neighbourhood size
##
##----------------------------------------------------------------------------##


##  USES trap.interval specified in 'execute.parasite.function' for grid structure 


make.dist.vector <- function(maxdist){ 

  #Create a sequence of all 1-d coordinates up to 'maxdist' from 0
  vals <- seq(0, maxdist, by=trap.interval)	
  
  #create a 2-d mini-grid anchored at hypothetical focal location {0,0}, with all coordinates from that point
  distmatrix <- as.data.frame(permutations(n = length(vals), r=2, v=vals, repeats.allowed = T)) 
  colnames(distmatrix) <- c("long", "lat") 
  
  #Calculate Euclidean distances of each point in distmatrix from {0,0} 
  distmatrix$distfromfoc <- sqrt((0 - distmatrix$long)^2 + (0 - distmatrix$lat)^2) 
  
  #Extract the unique distances across the mini-grid, and sort by distance
  dist <- unique(distmatrix$distfromfoc)
  dist <- sort(dist)

  #Extract only those distances that lie within 'maxdist' of the focal
  dist <- dist[which (dist <= maxdist)]
  return(dist)
}
