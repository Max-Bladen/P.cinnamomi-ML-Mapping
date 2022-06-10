### ======================================================================== ###
### ------------------------------- LIBRARIES ------------------------------ ###
### ======================================================================== ###

library(DescTools)
library(raster)
library(readxl)
library(sp)

### ======================================================================== ###
### ------------------------------- FUNCTIONS ------------------------------ ###
### ======================================================================== ###


### ======================================================================== ###
## slope: df of three columns. long, lat and val which is the slope value at that pos
## aspect: df of three columns. long, lat and val which is the aspect value at that pos
## sunPos: df of elevation and azimuth for each hour averaged across each month
#
# Iterates over all slope and aspect value pairs. For each, finds average HS 
# by iterating over all sun positions for the year. 
#
## return: df of three columns. long, lat and val which is the PRR value at that pos
CalculatePRR <- function(slope, 
                         aspect, 
                         sunPos) 
  {
  
  if (dim(slope)[1] != dim(aspect)[1]){ # ensure same nrow of 'slope' and 'aspect'
    stop("Differing dimensions of aspect and slope")
  }
  
  # set output df, same format as 'slope' and 'aspect'
  output <- data.frame(matrix(NA, nrow = nrow(slope), ncol = 3))
  colnames(output) <- c("long", "lat", "val")
  output$long <- slope$long
  output$lat <- slope$lat
  
  # for each row ...
  for (r in 1:nrow(slope)) {
    slp <- slope[r,]$val
    asp <- aspect[r,]$val
    
    hS <- 0
    
    for (s in 1:nrow(sunPos)) { # and for each sun measurement ...
      zenith <- (90 - sunPos[s,"E"]) * pi/180
      azimuth <- sunPos[s,"A"] * pi/180
      
      hS <- hS + HillShade(slp, asp, azimuth, zenith) # calculate cumulative hS
    }
    output[r, "val"] <- hS/nrow(sunPos) # return average for that cell
    
    p <- floor(r*100/nrow(slope))
    cat("\rCalculating PRR Progress: ", p, "%,  ")
  }
  
  return(output)
}


### ======================================================================== ###
## slope: df of three columns. long, lat and val which is the slope value at that pos
## aspect: df of three columns. long, lat and val which is the aspect value at that pos
#
# Iterates over all slope and aspect value pairs. For each, calculates Sun Index
# using formula.
#
## return: df of three columns. long, lat and val which is the Sun.Index value at that pos
CalculateSunIndex <- function(slope, 
                              aspect) 
  {
  
  if (dim(slope)[1] != dim(aspect)[1]){ # ensure same nrow of 'slope' and 'aspect'
    stop("Differing dimensions of aspect and slope")
  }
  
  # set output df, same format as 'slope' and 'aspect'
  output <- data.frame(matrix(NA, nrow = nrow(slope), ncol = 3))
  colnames(output) <- c("long", "lat", "val")
  output$long <- slope$long
  output$lat <- slope$lat
  
  # calculate Sun Index
  output$val <- cos(aspect$val) * tan(slope$val) * 100
  return(output)
}


### ======================================================================== ###
## dfs: list of 3-column dfs (long, lat, val), each of which correspond to a feature. 
#
# Iterates over all individual feature dataframes and combines them in a column-
# wise fashion.
# 
## return: df containing all input feature values for the long and lat coordinates.
##         ncol(df) = length(dfs), nrow(df) = nrow(dfs[[1]])
CondenseRasterDFs <- function(dfs) 
  {
  
  df <- dfs[[1]][, c("long", "lat")]
  n <- names(dfs)
  for (i in 1:length(dfs)) {
    df[, n[i]] <- dfs[[i]][, "val"]
  }
  return(df)
}


### ======================================================================== ###
## DN1: 3-column df (long, lat, val) representing a reflectance ratio feature
## DN2: 3-column df (long, lat, val) representing a different reflectance ratio feature
#
# Uses NDRI formula for the two input DNs. Does this down the whole DF
#
## return: 3-column df (long, lat, val) representing a NDRI measure
DF.NDRI <- function(DN1, 
                    DN2) 
  {
  if (length(which(DN1$long != DN2$long))!=0) {
    stop("long columns are not equal")
  }
  if (length(which(DN1$lat != DN2$lat))!=0) {
    stop("lat columns are not equal")
  }
  
  new.df <- cbind(DN1[, c(1,2)], (DN1[,3]-DN2[,3])/(DN1[,3]+DN2[,3]))
  colnames(new.df) <- c("long", "lat", "val")
  return(new.df)
}


### ======================================================================== ###
## r: 'raster' type object
#
# Converts raster to a "long" df, such that each row represents a pair of 
# present long and lat coords with the associated raster value at that point.
#
## return: 3-column df (long, lat, val) representing values of input raster 
##         at each cell
Extract.RasterValues.To.LongDF <- function(r) 
  {
  
  r.sptl.pnts <- as.data.frame(r,xy=TRUE); coordinates(r.sptl.pnts)=~x+y
  df <- as.data.frame(r.sptl.pnts)
  colnames(df) <-  c("long", "lat", "val")
  
  return(df)
}


### ======================================================================== ###
## slp: slope of landscape, val or vector
## asp: aspect of landscape, val or vector
## azi: azimuth of sun, val or vector
## zen: zenith of sun, val or vector
#
#
#
## return: 
HillShade <- function(slp, 
                      asp, 
                      azi, 
                      zen) 
  {
  return(255 * (cos(slp) * cos(zen) + sin(slp) * sin(zen) * cos(azi - asp)))
}



### ======================================================================== ###
## inputDF: sun data in the form of average elevation (zenith) and azimuth
##          values at each hour for a set of days.
#
# Converts the input df in a new df. Contains the average zenith and azimuth
# values for each hour of each month. Cleans the dataframe up by removed "--"
# which represent null in the 'inputDF'
#
## return: Long form of input data with invalid rows removed
ReduceSunData <- function(inputDF) 
  {
  
  ### ------- PRE-PROCESS DATA ------- ###
  colnames(inputDF)[1] <- "date" # change first column name to "date"
  inputDF <- inputDF[-nrow(inputDF), ] # remove last row (from 2023)
  
  # extract date into "day" and "month" column
  inputDF[, "month"] <- substr(inputDF[,"date"], 6, 7)
  inputDF[, "day"] <- substr(inputDF[,"date"], 9, 10)
  
  # replace all "--" with NA
  is.na(inputDF) <- inputDF == "--" 
  
  # remove columns with all NAs
  NAcols <- unname(which(sapply(inputDF, function(x)all(is.na(x)))))
  inputDF <- inputDF[, -NAcols]
  
  
  ### ------- REDUCE MONTHS ------- ###
  
  # set months list
  monthsList <- unique(inputDF$month)
  
  # produce dataframe which takes mean for each time period across all days
  # in each month. note the -3 in ncol, removes the date, month and day feature
  averagedDF <- data.frame(matrix(NA, nrow = 12, ncol = ncol(inputDF)-3)) 
  
  # set dimension names accordingly
  rownames(averagedDF) <- monthsList
  colnames(averagedDF) <- colnames(inputDF)[!colnames(inputDF) 
                                            %in% c("date", "month", "day")]
  
  # iterate over each month, extracting data for that month then find mean for each
  # time (both elevation and azimuth). pass this mean value to the averagedDF df
  m <- 1
  for (month in monthsList) {
    monthData <- inputDF[inputDF$month == month, ]
    for (c in colnames(averagedDF)) {
      averagedDF[m, c] <- mean(as.numeric(monthData[,c]), na.rm = TRUE)
    }
    m<-m+1
  }
  
  # remove extraneous information from column names 
  colnames(averagedDF) <- substr(colnames(averagedDF), 1, 4)
  
  ### ------- REFORMAT REDUCED MONTHS ------- ###
  
  # extract the columns which correspond to elevation and azimuth
  measurement <- substr(colnames(averagedDF), 1, 1)
  # extract the columns which correspond to time
  time <- substr(colnames(averagedDF), 3, 4)
  
  # set empty dataframe. /2 is due to the alternating "E" (elevation) and "A" (azimuth)
  outputDF <- data.frame(matrix(NA, 
                                nrow = length(monthsList)*length(time)/2, 
                                ncol = 4))
  colnames(outputDF) <- c("month", "time", "E", "A")
  
  # set the month and time features of the dataframe
  outputDF$month <- unlist(lapply(monthsList, rep, 
                                  times = length(unique(time))))
  outputDF$time <- rep(unique(time), length(monthsList))
  
  # for both elevation and azimuth, extract the relevant columns, set a list with
  # all values for that measurement and pass to the final dataframe
  for (measure in unique(measurement)) {
    cols <- unname(which(stri_detect_fixed(colnames(averagedDF), measure)))
    measureCol <- c()
    # for each month
    for (m in rownames(averagedDF)) {
      measureCol <- c(measureCol, unlist(averagedDF[m, cols]))
    }
    outputDF[, measure] <- measureCol
  }
  
  finalRowsToRemove <- which(is.na(outputDF$E))
  outputDF <- outputDF[-finalRowsToRemove, ]
  
  finalRowsToRemove <- which(outputDF$E < 0)
  outputDF <- outputDF[-finalRowsToRemove, ]
  
  rownames(outputDF) <- 1:length(rownames(outputDF))
  
  return(outputDF)
}


### ======================================================================== ###
## mat: matrix of values
## lowCol: color to use for the minimum 'mat' value
## highCol: color to use for the mmximum 'mat' value
#
# Plots a heatmap of a spatial value. Needs an matrix not a long dataframe
#
## return: 'pheatmap' object (see Value of ?pheatmap)
Spatial.HM <- function(mat, 
                       lowCol = "blue", 
                       highCol = "red") 
  {
  pheatmap(mat, cluster_rows = F, cluster_cols = F, 
           colorRampPalette(c(lowCol, "white", highCol))(50),
           show_rownames = F, show_colnames = F)
}