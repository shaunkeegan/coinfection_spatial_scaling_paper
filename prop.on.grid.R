## ---------------------- prop.on.grid()
## takes a focal's lat and long, and a specified distance (d)
##  returns the proportion of trap locations on grid within distance d of the focal's location
##



prop.on.grid <- function(focal.lat, focal.long, d, dnumpoint){
  
  #create a dataframe with all coordinates in a (min.point to max.point) x (min.point to max.point) grid
  vals3 <- seq(min.point, max.point, by=10)	
  distmatrix3 <- as.data.frame(permutations(n = length(vals3), r=2, v=vals3, repeats.allowed = T)) # n=vec size; r=number to choose each time
  colnames(distmatrix3) <- c("lat", "long") # Give them sensible names
  

  
  tempmat <- distmatrix3
  tempmat$distfromfoc <- sqrt((focal.lat - tempmat$lat)^2 + (focal.long - tempmat$long)^2)	#distance of each gridpoint from focal
  temp <- tempmat[which(tempmat$distfromfoc <= d),]		#subset out just those locations that are within 'd' of the focal
  val <- dnumpoint[dnumpoint$dist_vec == d, "numpoints"]
  ifelse(length(val)==0, return(0), return(nrow(temp)/val))	#

}