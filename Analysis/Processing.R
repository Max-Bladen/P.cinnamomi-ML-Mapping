
### ======================================================================== ###
### ------------------------------- DIRECTORY ------------------------------ ###
### ======================================================================== ###

wd = "C:/...P. cinnamomi ML Mapping/Analysis" # SET YOUR WORKING DIRECTORY 
setwd(wd)
source("Internals/Internals.Processing.Feat.Eng.R")

### ======================================================================== ###
### -------------------------- READ IN LABEL DATA -------------------------- ###
### ======================================================================== ###

### ------------------------------------------------------------------------ ###

# extract site 1 data sheets from xlsx
d.sheet1 <- data.frame(as.matrix(read_excel("Data/Labelled.Data/Raw.Site.1.Data.xlsx")))
d.sheet2 <- data.frame(as.matrix(read_excel("Data/Labelled.Data/Raw.Site.1.Data.xlsx", sheet = 2)))

# extract site 2 data sheets from xlsx
d.sheet3 <- data.frame(as.matrix(read_excel("Data/Labelled.Data/Raw.Site.2.Data.xlsx")))


### ------------------------------ d.sheet1 -------------------------------- ###

# remove NAs
which(is.na(d.sheet1), arr.ind = T) # showed rows 151 & 152 have NAs
d.sheet1 <- d.sheet1[1:150,]

site1.in.polygon <- d.sheet1[, c(5, 4, 2)] # extract long, lat and disease
colnames(site1.in.polygon) <- c("long", "lat", "disease")


### ------------------------------ d.sheet2 -------------------------------- ###

# rows corresponding to the desired polygon of points
# -1 due to conversion from excel to dataframe
polygon.rows <- 805:1060 - 1 

# all points of this polygon are post infested, hence extract long & lat
site1.of.polygon <- cbind(d.sheet2[polygon.rows, c(13, 12)], rep("P", length(polygon.rows)))
colnames(site1.of.polygon) <- c("long", "lat", "disease")


### ------------------------------ d.sheet3 -------------------------------- ###

 # extract long, lat and disease
site2.polygon <- d.sheet3[, c(5, 4, 1)]
colnames(site2.polygon) <- c("long", "lat", "disease")

# bring class labels in line with site 1 data
site2.polygon$disease[site2.polygon$disease == "Active Dieback front"] <- "A"
site2.polygon$disease[site2.polygon$disease == "Post Infested"] <- "P"
site2.polygon$disease[site2.polygon$disease == "Dieback free"] <- "N"

which(is.na(site2.polygon)) # check for NAs (there are none)


### ======================================================================== ###
### -------------------------- SAVE LABELLED DATA -------------------------- ###
### ======================================================================== ###

### ---------------------------- Save per site ----------------------------- ###
# write.csv(site1.in.polygon, file="Data/site1.in.polygon.csv")
# write.csv(site1.of.polygon, file="Labelled.Data/site1.of.polygon.csv")
# write.csv(site2.polygon, file="Labelled.Data/site2.polygon.csv")

### --------------------------- Save per class ----------------------------- ###

# active.points <- all.points[all.points$disease=="A", ]
# post.points <- all.points[all.points$disease=="P", ]
# free.points <- all.points[all.points$disease=="N", ]
# 
# write.csv(active.points, file="Labelled.Data/Subsets/active.points.csv")
# write.csv(post.points, file="Labelled.Data/Subsets/post.points.csv")
# write.csv(free.points, file="Labelled.Data/Subsets/free.points.csv")

### --------------------------- Save all ------------------------------ ###

all.points <- rbind(site1.in.polygon, site1.of.polygon, site2.polygon)
write.csv(all.points, file="Data/Labelled.Data/Long.Lat.Class.csv")















### ======================================================================== ###
### -------------------- CONDENSE LABELS AND RASTER DATA ------------------- ###
### ======================================================================== ###

Desired.CRS <- crs(raster("Data/Spatial.Data/Clipped.DEM.tif"))

IVs.DF <- CondenseRasterDFs(readRDS("Data/All.Raster.DFs.rda"))
DV.DF <- read.csv("Data/Labelled.Data/Long.Lat.Class.csv")[, 2:4]

### ======================================================================== ###

# convert labels to integer
DV.DF$disease[DV.DF$disease=="A"] <- 1
DV.DF$disease[DV.DF$disease=="P"] <- 2
DV.DF$disease[DV.DF$disease=="N"] <- 3

# clean row names
DV.DF <- DV.DF[-which(duplicated(DV.DF)),]
rownames(DV.DF) <- 1:nrow(DV.DF)

### ======================================================================== ###

# adjust coordinate reference system
pts <- SpatialPointsDataFrame(coords = DV.DF[,c(1,2)], data = DV.DF,
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
pts.lyr <- spTransform(pts, Desired.CRS)

### ======================================================================== ###

# combine IV and DV dataframes
IV.DV.DF <- DV.DF
IV.DV.DF$long <- pts.lyr$long
IV.DV.DF$lat <- pts.lyr$lat

# convert long and lat to x and y
IV.DV.DF$x <- floor(pts.lyr$long - min(IVs.DF$long))
IV.DV.DF$y <- floor(pts.lyr$lat - min(IVs.DF$lat))

IV.DV.DF$disease <- as.numeric(IV.DV.DF$disease) # ensure not factor

IV.DV.DF <- IV.DV.DF[, c(1,2,4,5,3)]

### ======================================================================== ###

for(r in 1:nrow(IV.DV.DF)) { # 1:446
  
  for (var in colnames(IVs.DF)[3:21]) {
  
    point <- c(IV.DV.DF$long[[r]], IV.DV.DF$lat[[r]])
    
    rows.x <- which(round(abs(point[1]-IVs.DF$long)) == min(round(abs(point[1]-IVs.DF$long))))
    rows.y <- which(round(abs(point[2]-IVs.DF$lat)) == min(round(abs(point[2]-IVs.DF$lat))))
    
    IV.DV.DF[r, var] <- mean(IVs.DF[intersect(rows.x, rows.y), var])
    
  }
  cat("\r", r)
}

saveRDS(IV.DV.DF, file = "Data/Labelled.Predictors.rda")
