##  Collates NH variables: NemInf/Uninf/NemPrev, Nemtotbuden/Nemmeanburden; Eiminf/Eimuninf/EimPrev for each focal entry at each NH size
##

##----------------------------------------------------------------------------##
##                              ARGUMENTS                                     ##
##----------------------------------------------------------------------------##
##   
##  neighbourdat: (dataframe) - neighbour dataset for characterising neighbourhood infection levels around each focal
##  focaldat: (dataframe) - focal dataset, including individual-level variables for each focal
##  neighbour.variable.INF: (character) - which variable in neighbourdat relates to neighbour nematode infection status (presence/absence) 
##  neighbour.variable.EPG: (character) - which variable in neighbourdat relates to neighbour nematode EPG
##  focalsp.NH.INF: (character) - which variable in neighbourdat relates to neighbour Eimeria infection status (presence/absence)
##  dist (numeric vector) - all possible  neighbourhood sizes (in m) around focal
##  time.window (character or numeric) - the time window to consider neighbours within ("ALL" indicates any time)
##
##----------------------------------------------------------------------------##



parasite.function.vary.time <- function (neighbourdat, focaldat, neighbour.variable.INF, neighbour.variable.EPG, focalsp.NH.INF, dist, time.window){

  for (NH.dist in 1:length(dist)){ # cycle through all possible neighbourhood sizes within 'dist'
    
    no <- length(unique(focaldat$ID.CapDate)) # Number of unique ID.CapDate entries to process
    
    for (i in 1:no){ # loop through all unique ID.CapDate entries
      
      #Set NH parasite cumulative variables to 0 for current focal
      Neminfsum <- 0
      Nemuninfsum <- 0
      Nemburdsum <- 0  
      
      Eiminfsum <- 0   
      Eimuninfsum <- 0

      
      #Get X,Y,date for current focal individual
      focal.X <- focaldat[i,"X"]
      focal.Y <- focaldat[i,"Y"]
      focal.date <- focaldat[i,"Capture.date"]
      
      
      # Calculate Percentage of the NH in the grid
      #1) All possible distances to trap locations within current NH size
      dist_vec <- make.dist.vector(dist[NH.dist])
      dist_vec <- sort(dist_vec)
      
      #2) generate list of number of points for each value of 'dist' from central point
      numpoints <- apply(as.array(dist_vec), 1, exppoints)	#gives a list 1, 5, 9, 13, 21, ... for all values in 'dist'

      #3) bind that list to dist_vec and round to 1 dp
      dnumpoint <- as.data.frame(cbind (dist_vec, numpoints))
      dnumpoint$dist_vec <- round(dist_vec, digits = 1)
      
      #4) calculate proportion of neighbourhood on the grid around the focal, and add to focaldat df
      focaldat[i,PropOnGrid[NH.dist]] <- prop.on.grid(focal.X, focal.Y, dist_vec[NH.dist], dnumpoint)
      
      
      #Subset neighbourhood data to the focal's Grid and Year of capture
      GRID <- as.character(focaldat[i, "Grid"])
      YEAR <- focaldat[i,"year"]
      neighbourdat.subset <- neighbourdat[which (neighbourdat$Grid == GRID & neighbourdat$year == YEAR),] 
      
      
      #Check what time lag to consider, subset neighbourdat accordingly, and process
      
      if(time.window == "ALL"){
        for (j in 1:length(neighbourdat.subset$ID.CapDate)){ # Loop through each neighbour in turn
          neighbour.Y <- neighbourdat.subset[j,"Y"]
          neighbour.X <- neighbourdat.subset[j,"X"]
          temp.date <- neighbourdat.subset[j,"Capture.date"]
          
          #Calculate neighbour's capture distance from focal
          neighbourdist <- sqrt((focal.X - neighbour.X)^2 + (focal.Y - neighbour.Y)^2)
          
          #If neighbour lies within current NH size and was captured before focal, process...
          if(neighbourdist <= dist[NH.dist] &
             temp.date < focal.date ){
            
            #Accumulate nem infection status (number of infecteds) in NH
            ifelse(neighbourdat.subset[j,neighbour.variable.INF]==1, Neminfsum <- Neminfsum + 1, 
                   ifelse(neighbourdat.subset[j,neighbour.variable.INF]==0, Nemuninfsum <- Nemuninfsum + 1, NA))	
            #Accumulate nem infection burden (mean EPG) in NH
            ifelse(is.na(neighbourdat.subset[j,neighbour.variable.EPG]), NA, 
                   Nemburdsum <- Nemburdsum + neighbourdat.subset[j,neighbour.variable.EPG])
            #Accumulate focal infection status (number of infecteds) in NH
            ifelse(neighbourdat.subset[j,focalsp.NH.INF]==1, Eiminfsum <- Eiminfsum + 1, 
                   ifelse(neighbourdat.subset[j,focalsp.NH.INF]==0, Eimuninfsum <- Eimuninfsum + 1, NA))	
            
          }	#Close if neighbour is to be processed
        }	#Close j (neighbour) loop 
        
      }else{  #time.window!="ALL"
        
        for (j in 1:length(neighbourdat.subset$ID.CapDate)){ # Loop through each neighbour
          neighbour.Y <- neighbourdat.subset[j,"Y"]
          neighbour.X <- neighbourdat.subset[j,"X"]
          temp.date <- neighbourdat.subset[j,"Capture.date"]
          window.date <- as.Date(focal.date) - time.window
          
          #Calculate neighbour's capture distance from focal
          neighbourdist <- sqrt((focal.X - neighbour.X)^2 + (focal.Y - neighbour.Y)^2)
          
          #If neighbour lies within current NH size and was captured before focal but within window.date, process...
          if(neighbourdist <= dist[NH.dist] &
             temp.date < focal.date &
             temp.date >= window.date){
            
            #Accumulate nem infection status (number of infecteds) in NH
            ifelse(neighbourdat.subset[j,neighbour.variable.INF]==1, Neminfsum <- Neminfsum + 1, 
                   ifelse(neighbourdat.subset[j,neighbour.variable.INF]==0, Nemuninfsum <- Nemuninfsum + 1, NA))	
            #Accumulate nem infection burden (mean EPG) in NH
            ifelse(is.na(neighbourdat.subset[j,neighbour.variable.EPG]), NA, 
                   Nemburdsum <- Nemburdsum + neighbourdat.subset[j,neighbour.variable.EPG])
            #Accumulate focal infection status (number of infecteds) in NH
            ifelse(neighbourdat.subset[j,focalsp.NH.INF]==1, Eiminfsum <- Eiminfsum + 1, 
                   ifelse(neighbourdat.subset[j,focalsp.NH.INF]==0, Eimuninfsum <- Eimuninfsum + 1, NA))
            
          }	#Close if neighbour is to be processed 
        }	#Close j (neighbour) loop 
      }
      
      
      #Collate all output variables
      Nemtotnum <- Neminfsum + Nemuninfsum
      Eimtotnum <- Eiminfsum + Eimuninfsum
      focaldat[i,NemNoInf[NH.dist]] <- Neminfsum
      focaldat[i,NemNoUninf[NH.dist]] <- Nemuninfsum
      focaldat[i,NemPrev[NH.dist]] <- ifelse(Nemtotnum > 0, Neminfsum/Nemtotnum, 0)
      focaldat[i, NemNoTotal[NH.dist]] <- Nemtotnum / focaldat[i, PropOnGrid[NH.dist]]
      
      #Note, these are logged + 1 burdens
      focaldat[i,NemTotBurden[NH.dist]] <- log(Nemburdsum+1)
      focaldat[i,NemMeanBurden[NH.dist]] <- ifelse(Nemtotnum > 0, log(Nemburdsum/Nemtotnum+1), 0)
      
      focaldat[i,EimNoInf[NH.dist]] <- Eiminfsum
      focaldat[i,EimNoUninf[NH.dist]] <- Eimuninfsum
      focaldat[i,EimPrev[NH.dist]] <- ifelse(Eimtotnum > 0, Eiminfsum/Eimtotnum, 0)
      focaldat[i, EimNoTotal[NH.dist]] <- Eimtotnum / focaldat[i, PropOnGrid[NH.dist]] 
      
    }
    
    
  }
  
  return(focaldat)
  
}