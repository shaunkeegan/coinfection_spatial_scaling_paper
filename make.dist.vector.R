## Creates vector of all possible neighbourhood sizes 
##  from 0 to specified 'max.point' on a trapping grid with trap locations 'trap.interval'. apart
##


make.dist.vector <- function(max.point, trap.interval){ # Define function
  
  ##----------------------------------------------------------------------------##
  ##                                  ARGS                                      ##
  ##----------------------------------------------------------------------------##
  ##  max.point - Maximum radius (in m) around focal animal                     ##
  ##  trap.interval - Distances (in m) between trapping points on grid          ##
  ##----------------------------------------------------------------------------##
  
 
  
  
  vals3 <- seq(0, max.point, by=trap.interval)	
  distmatrix3 <- as.data.frame(permutations(n = length(vals3), r=2, v=vals3, repeats.allowed = T)) # n=vec size; r=number to choose each time
  colnames(distmatrix3) <- c("lat", "long") # Give them sensible names
  
  tempmat <- distmatrix3
  tempmat$distfromfoc <- sqrt((0 - tempmat$lat)^2 + (0 - tempmat$long)^2) 
  
  dist <- unique(tempmat$distfromfoc)
  dist <- sort(dist)
  dist <- dist[which (dist <= max.point)]
  return(dist)
}
