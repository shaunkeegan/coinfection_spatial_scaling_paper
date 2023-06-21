## ---------------------- parasite.function.vary.time()
##
##
##
##

parasite.function.vary.time <- function (dat, dat2, neighbour.variable, dist, time.window){
  ## -------------------- inputs ---
  ## dat            - dataset of neighbours
  ## dat2           - dataset of focals
  ## neighbour.variable  - the column in 'dat' for neighbour infection status by the effecting parasite  
  ## dist           - the list of neighbourhood sizes to work through
  ## time.window    - the desired window to consider neighbour effects prior to the focal's capture date
  ## -------------------------------
  
  for (NH.dist in 1:length(dist)){ # cycle through all possible neighbourhood sizes
    
    no <- length(unique(dat2$ID.CapDate)) # All ID.CapDate entries
    
    for (i in 1:no){ # loop through all ID.CapDate entries
      
      infsum <- 0
      uninfsum <- 0
      
      #Get X,Y,date for focal individual
      focal.X <- dat2[i,"X"]
      focal.Y <- dat2[i,"Y"]
      focal.date <- dat2[i,"Capture.date"]
      
      #Create vector of all possible neighbourhood sizes 
      dist_vec <- make.dist.vector(dist[NH.dist], trap.interval)
      dist_vec <- sort(dist_vec)
      
      #generate list of number of points for each value of 'dist' from central point
      numpoints <- apply(as.array(dist_vec),1,exppoints)	#gives a list 1, 5, 9, 13, 21, ... for all values in 'dist_vec'
      #bind that list to dist_vec
      dnumpoint <- as.data.frame(cbind (dist_vec, numpoints))
      dnumpoint$dist_vec <- round(dist_vec, digits = 1)
      dat2[i,PropOnGrid[NH.dist]] <- prop.on.grid(focal.X, focal.Y, dist_vec[NH.dist], dnumpoint["dist_vec"])
      
      
      # Subset animals on same grid in same year 
      GRID <- as.character(dat2[i, "Grid"])
      YEAR <- dat2[i,"year"]
      dat.subset <- dat[which (dat$Grid == GRID & dat$year == YEAR),]
      
      #Go through each non-focal, work out distance from focal, 
      # if within specified neighbourhood and specified capture window of focal,
      # increment number of animals and infected animals, depending on infection status
      
      if(time.window == "ALL"){
        for (j in 1:length(dat.subset$ID.CapDate)){ # Loop through all other animals
          neighbour.Y <- dat.subset[j,"Y"]
          neighbour.X <- dat.subset[j,"X"]
          temp.date <- dat.subset[j,"Capture.date"]
          
          neighbourdist <- sqrt((focal.X - neighbour.X)^2 + (focal.Y - neighbour.Y)^2)
          
          if(neighbourdist <= dist[NH.dist] &
             temp.date < focal.date ){
            
            ifelse(dat.subset[j,neighbour.variable]==1, infsum <- infsum + 1, 
                   ifelse(dat.subset[j,neighbour.variable]==0, uninfsum <- uninfsum + 1, NA))	
            
          }	#Close if
        }	#Close j loop
      }else{
        for (j in 1:length(dat.subset$ID.CapDate)){ # Loop through all other animals
          neighbour.Y <- dat.subset[j,"Y"]
          neighbour.X <- dat.subset[j,"X"]
          temp.date <- dat.subset[j,"Capture.date"]
          window.date <- as.Date(focal.date) - time.window
          
          neighbourdist <- sqrt((focal.X - neighbour.X)^2 + (focal.Y - neighbour.Y)^2)
          
          if(neighbourdist <= dist[NH.dist] &
             temp.date < focal.date &
             temp.date >= window.date){
            
            ifelse(dat.subset[j,neighbour.variable]==1, infsum <- infsum + 1, 
                   ifelse(dat.subset[j,neighbour.variable]==0, uninfsum <- uninfsum + 1, NA))	
            
          }	#Close if
        }	#Close j loop
      }
      
      
      #Update relevant neighbourhood count variables in main dataset
      dat2[i,NoInf[NH.dist]] <- infsum
      dat2[i,NoUninf[NH.dist]] <- uninfsum
      dat2[i,Prev[NH.dist]] <- infsum/(infsum + uninfsum)
      dat2[i, NoTotal[NH.dist]] <- infsum + uninfsum 
      
      
    }
    
    
  }
  
  return(dat2)
  
}