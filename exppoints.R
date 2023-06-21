## ---------------------- exppoints()
## Takes value 'distance' 
##  returns the number of trap locations within that distance of the central point
##



exppoints <- function(distance){
  
  #create a minigrid centred on {0,0}, with coordinates within x metres in all directions of that central point
  vals2 <- seq(-max.point,max.point,by=trap.interval)
  distmatrix2 <- as.data.frame(permutations(n=length(vals2),r=2,v=vals2,repeats.allowed=T))            #n=vec size; r=number to choose each time
  colnames(distmatrix2) <- c("lat", "long")                 #give them sensible names
  distmatrix2$dist <- sqrt((0 - distmatrix2$lat)^2 + (0 - distmatrix2$long)^2)
  
  return(nrow(distmatrix2[distmatrix2$dist<=distance,]))

}
