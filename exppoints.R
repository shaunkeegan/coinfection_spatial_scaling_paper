##  Returns the number of trap locations within specified distance of the central point
##

##----------------------------------------------------------------------------##
##                              ARGUMENTS                                     ##
##----------------------------------------------------------------------------##
##
##  distance: (numeric) - radius (in m) of the desired neighbourhood 
##
##----------------------------------------------------------------------------##


##  USES max.point and trap.interval specified in 'execute.parasite.function' for grid structure 


exppoints <- function(distance){
  
  #Create a sequence of all 1-d coordinates up to 'max.point' in both directions of central point
  vals <- seq(-max.point, max.point, by=trap.interval)
  
  #create a 2-d mini-grid centred on {0,0}, with coordinates in all directions of that central point
  distmatrix <- as.data.frame(permutations(n=length(vals), r=2, v=vals, repeats.allowed=T))
  colnames(distmatrix) <- c("long", "lat")  
  
  #Calculate Euclidean distances of each point in distmatrix from the central point 
  distmatrix$dist <- sqrt((0 - distmatrix$long)^2 + (0 - distmatrix$lat)^2)
  
  #Return the number of points lying within specified distance of the central point
  return(nrow(distmatrix[distmatrix$dist<=distance,]))

}
