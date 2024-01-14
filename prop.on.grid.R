##  Calculates the proportion of trap locations on grid within distance 'd' of a focal's specified location
##

##----------------------------------------------------------------------------##
##                              ARGUMENTS                                     ##
##----------------------------------------------------------------------------##
##   
##  focal.long: (numeric) - The x coordinate of the focal's location
##  focal.lat: (numeric) - The y coordinate of the focal's location
##  d: (numeric) - The specified neighbourhood size (in m) around the focal
##  dnumpoint: (dataframe) - comprising all possible neighbourhood sizes around the focal (column 1)
##      and the number of trap locations lying within those neighbourhood sizes (column 2)
##
##----------------------------------------------------------------------------##


##  USES max.point and trap.interval specified in 'execute.parasite.function' for grid structure 


prop.on.grid <- function(focal.long, focal.lat, d, dnumpoint){
  
  #Create a sequence of all 1-d coordinates arbitrarily up to twice 'max.point' from point 0
  vals <- seq(0, 2*max.point, by=trap.interval)	
  
  #create a 2-d mini-grid anchored at hypothetical focal location {0,0}, with all coordinates from that point
  distmatrix <- as.data.frame(permutations(n = length(vals), r=2, v=vals, repeats.allowed = T)) 
  colnames(distmatrix) <- c("long", "lat") 
  
  #Calculate Euclidean distances of each point in distmatrix from the focal's location
  distmatrix$distfromfoc <- sqrt((focal.long - distmatrix$long)^2 + (focal.lat - distmatrix$lat)^2)
  
  #subset out just those locations that are within distance 'd' of the focal
  temp <- distmatrix[which(distmatrix$distfromfoc <= d),]	
  
  ##select from dnumpoint all the points that would lie within distance d of the focal on an infinite grid
  val <- dnumpoint[dnumpoint$dist_vec == round(d,1), "numpoints"]
  
  #return the proportion of possible points that actually lie on the grid
  ifelse(length(val)==0, return(0), return(nrow(temp)/val))	

}