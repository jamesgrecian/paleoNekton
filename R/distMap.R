###############
### distMap ###
###############

# Climate change velocity using multivariate analogs,
# method described in:
# Ordonez & Williams (2013) DOI: 10.1007/s10584-013-0752-1

# function taken from Ryan Reisinger github repo
#https://github.com/ryanreisinger/PEIfuture/Working/Scripts/old/10-Compare climates.R


distMap <- function(historical,
                    future,
                    thr = 0.5,
                    which.distance = "SQeuclidean") {
  
  require(analogue)
  require(geosphere)
  
  pb <- txtProgressBar(min = 0, max = nrow(historical), style = 3)
  
  frm <- data.frame()
  
  # Add an ID column for matching
  historical$ID <- 1:nrow(historical)
  future$ID <- 1:nrow(future)
  
  # Get the future conditions
  fut <- future[, 3:ncol(future)]
  fut <- fut[complete.cases(fut), ]
  futIDs <- fut$ID
  fut$ID <- NULL
  
  for (i in 1:nrow(historical)) {
    
    setTxtProgressBar(pb, i)
    
    his <- historical[i, 3:ncol(historical)]
    his <- his[complete.cases(his), ]
    hisIDs <- his$ID
    his$ID <- NULL
    
    # Run only if the cell in question has no missing values
    if(nrow(his) > 0) {
      
      d <- analogue::distance(his, fut, method = which.distance)
      
      nnval <- min(d)
      if(nnval < thr) {
        nndx <- which(d == min(d))
        nnid <- futIDs[nndx]
        
        # Get the data for the nearest neighbour
        this.future.cell <- future[future$ID == nnid, ]
        
        # Calculate distance and bearing to the nearest neighbour
        this.historical.cell <- historical[i, c("x", "y")]
        
        d <- distGeo(p1 = this.historical.cell[ , c("x", "y")], p2 = this.future.cell[ , c("x", "y")])/1000
        b <- bearing(p1 = this.historical.cell[ , c("x", "y")], p2 = this.future.cell[ , c("x", "y")])
        target.lon <- this.future.cell$x
        target.lat <- this.future.cell$y
      } else {
        d <- NA
        b <- NA
        target.lon <- NA
        target.lat <- NA
      }
      
    } else {
      d <- NA
      b <- NA
      target.lon <- NA
      target.lat <- NA
    }
    giveback <- historical[i, ]
    giveback$ID <- NULL
    giveback$distance <- d
    giveback$bearing <- b
    giveback$target.lon <- target.lon
    giveback$target.lat <- target.lat
    
    frm <- rbind(frm, giveback)
    rm(giveback)
    rm(his)
  }
  
  return(frm)
  
  close(pb)
}