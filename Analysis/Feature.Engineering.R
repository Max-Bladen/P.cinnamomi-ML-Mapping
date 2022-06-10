
### ======================================================================== ###
### ------------------------------- DIRECTORY ------------------------------ ###
### ======================================================================== ###

wd = "C:/...P. cinnamomi ML Mapping/Analysis" # SET YOUR WORKING DIRECTORY 
setwd(wd)
source("Internals/Internals.Processing.FeatEng.R")

### ======================================================================== ###
### -------------------------- READ IN INPUT DATA -------------------------- ###
### ======================================================================== ###

### ------------------------------ Reflectance ----------------------------- ###

Blue.Raster <- raster("Data/Spatial.Data/Clipped.Blue.tif")
Green.Raster <- raster("Data/Spatial.Data/Clipped.Green.tif")
Red.Raster <- raster("Data/Spatial.Data/Clipped.Red.tif")

NIR.Raster <- raster("Data/Spatial.Data/Clipped.NIR.tif")

SWIR1.Raster <- raster("Data/Spatial.Data/Clipped.SWIR1.tif")
SWIR2.Raster <- raster("Data/Spatial.Data/Clipped.SWIR2.tif")

### ---------------------------- Brightness Temp. -------------------------- ###

TIRS1.Raster <- raster("Data/Spatial.Data/Clipped.TIRS1.tif")
TIRS2.Raster <- raster("Data/Spatial.Data/Clipped.TIRS2.tif")

### -------------------------- Terrain Properties -------------------------- ###

Aspect.Raster <- raster("Data/Spatial.Data/Clipped.Aspect.tif")
Slope.Raster <- raster("Data/Spatial.Data/Clipped.Slope.tif")

### ------------------------------- TPI & DEM ------------------------------ ###

DEM.Raster <- raster("Data/Spatial.Data/Clipped.DEM.tif")
TPI100.Raster <- raster("Data/Spatial.Data/Clipped.TPI100.tif")
TPI1000.Raster <- raster("Data/Spatial.Data/Clipped.TPI1000.tif")

### -------------------------- SUN POSITION DATA --------------------------- ###

sunData <- data.frame(as.matrix(read_excel("Data/Sun.Data/Sun2022.xlsx")))

sunData <- ReduceSunData(sunData)


### ======================================================================== ###
### ------------------------------- PROJECTION ----------------------------- ###
### ======================================================================== ###

All.Rasters <- list(Blue.Raster=Blue.Raster,
                    Green.Raster=Green.Raster,
                    Red.Raster=Red.Raster,
                    NIR.Raster=NIR.Raster,
                    SWIR1.Raster=SWIR1.Raster,
                    SWIR2.Raster=SWIR2.Raster,
                    TIRS1.Raster=TIRS1.Raster,
                    TIRS2.Raster=TIRS2.Raster,
                    Aspect.Raster=Aspect.Raster,
                    Slope.Raster=Slope.Raster,
                    DEM.Raster=DEM.Raster,
                    TPI100.Raster=TPI100.Raster,
                    TPI1000.Raster=TPI1000.Raster)

Reflectance.Rasters <- list(Blue.Raster=Blue.Raster,
                            Green.Raster=Green.Raster,
                            Red.Raster=Red.Raster,
                            NIR.Raster=NIR.Raster,
                            SWIR1.Raster=SWIR1.Raster,
                            SWIR2.Raster=SWIR2.Raster,
                            TIRS1.Raster=TIRS1.Raster,
                            TIRS2.Raster=TIRS2.Raster)

Geo.Rasters <- list(Aspect.Raster=Aspect.Raster,
                    Slope.Raster=Slope.Raster,
                    DEM.Raster=DEM.Raster,
                    TPI100.Raster=TPI100.Raster,
                    TPI1000.Raster=TPI1000.Raster)

# reproject all reflectance rasters to bring in line with geo rasters
Blue.Raster <- projectRaster(Blue.Raster, DEM.Raster)
Green.Raster <- projectRaster(Green.Raster, DEM.Raster)
Red.Raster <- projectRaster(Red.Raster, DEM.Raster)
NIR.Raster <- projectRaster(NIR.Raster, DEM.Raster)
SWIR1.Raster <- projectRaster(SWIR1.Raster, DEM.Raster)
SWIR2.Raster <- projectRaster(SWIR2.Raster, DEM.Raster)
TIRS1.Raster <- projectRaster(TIRS1.Raster, DEM.Raster)
TIRS2.Raster <- projectRaster(TIRS2.Raster, DEM.Raster)

All.Rasters <- list(Blue.Raster=Blue.Raster,
                    Green.Raster=Green.Raster,
                    Red.Raster=Red.Raster,
                    NIR.Raster=NIR.Raster,
                    SWIR1.Raster=SWIR1.Raster,
                    SWIR2.Raster=SWIR2.Raster,
                    TIRS1.Raster=TIRS1.Raster,
                    TIRS2.Raster=TIRS2.Raster,
                    Aspect.Raster=Aspect.Raster,
                    Slope.Raster=Slope.Raster,
                    DEM.Raster=DEM.Raster,
                    TPI100.Raster=TPI100.Raster,
                    TPI1000.Raster=TPI1000.Raster)

Reflectance.Rasters <- list(Blue.Raster=Blue.Raster,
                            Green.Raster=Green.Raster,
                            Red.Raster=Red.Raster,
                            NIR.Raster=NIR.Raster,
                            SWIR1.Raster=SWIR1.Raster,
                            SWIR2.Raster=SWIR2.Raster,
                            TIRS1.Raster=TIRS1.Raster,
                            TIRS2.Raster=TIRS2.Raster)



### ======================================================================== ###
### ------------------------------- EXTRACT DF ----------------------------- ###
### ======================================================================== ###

### ------------------------------ Reflectance ----------------------------- ###

Blue.DF <- Extract.RasterValues.To.LongDF(Blue.Raster)
Green.DF <- Extract.RasterValues.To.LongDF(Green.Raster)
Red.DF <- Extract.RasterValues.To.LongDF(Red.Raster)

NIR.DF <- Extract.RasterValues.To.LongDF(NIR.Raster)

SWIR1.DF <- Extract.RasterValues.To.LongDF(SWIR1.Raster)
SWIR2.DF <- Extract.RasterValues.To.LongDF(SWIR2.Raster)

### ---------------------------- Brightness Temp. -------------------------- ###

TIRS1.DF <- Extract.RasterValues.To.LongDF(TIRS1.Raster)
TIRS2.DF <- Extract.RasterValues.To.LongDF(TIRS2.Raster)

### -------------------------- Terrain Properties -------------------------- ###

Aspect.DF <- Extract.RasterValues.To.LongDF(Aspect.Raster)
Slope.DF <- Extract.RasterValues.To.LongDF(Slope.Raster)

### ------------------------------- TPI & DEM ------------------------------ ###

DEM.DF <- Extract.RasterValues.To.LongDF(DEM.Raster)
TPI100.DF <- Extract.RasterValues.To.LongDF(TPI100.Raster)
TPI1000.DF <- Extract.RasterValues.To.LongDF(TPI1000.Raster)



### ======================================================================== ###
### -------------------------- FEATURE ENGINEERING ------------------------- ###
### ======================================================================== ###

### -------------------------- Vegetative Indicies ------------------------- ###

# -- Plant Pigment Ratio -- #
PPR.DF <- DF.NDRI(Green.DF, Blue.DF)
#Spatial.HM(PPR.DF, "blue", "chartreuse")

# -- Plant Vigour Ratio -- #
PVR.DF <- DF.NDRI(Green.DF, Red.DF)
#Spatial.HM(PVR, "red", "chartreuse")

# -- Water Pigment Ratio -- #
WPR.DF <- DF.NDRI(NIR.DF, SWIR2.DF)
#Spatial.HM(WPR, "red", "blue")

# -- Normalised Difference Vegetation Index -- #
NDVI.DF <- DF.NDRI(NIR.DF, Red.DF)
#Spatial.HM(NDVI, "red", "darkorchid1")

### ------------------------- Sun-Terrian Indicies ------------------------- ###

# -- Convert to Radians  -- #
Aspect.DF$val <- DegToRad(Aspect.DF$val)
Slope.DF$val <- DegToRad(Slope.DF$val)

# -- Potential Relative Radiation (mean HillShade) -- #
PRR.DF <- CalculatePRR(Slope.DF, Aspect.DF, sunData)
#Spatial.HM(PRR)

# -- Sun Index -- #
Sun.Index.DF <- CalculateSunIndex(Slope.DF, Aspect.DF)
#Spatial.HM(Sun.Index)




### ======================================================================== ###
### ------------------- CONVERT ALL DATA INTO SINGLE DF -------------------- ###
### ======================================================================== ###

All.Raster.DFs <- list(Blue = Blue.DF,
                       Green = Green.DF,
                       Red = Red.DF,
                       NIR = NIR.DF,
                       SWIR1 = SWIR1.DF,
                       SWIR2 = SWIR2.DF,
                       TIRS1 = TIRS1.DF,
                       TIRS2 = TIRS2.DF,
                       DEM = DEM.DF,
                       TPI100 = TPI100.DF,
                       TPI1000 = TPI1000.DF,
                       PPR = PPR.DF,
                       PVR = PVR.DF,
                       WPR = WPR.DF,
                       NDVI = NDVI.DF,
                       Aspect = Aspect.DF,
                       Slope = Slope.DF,
                       PRR = PRR.DF,
                       Sun.Index = Sun.Index.DF
)

#saveRDS(All.Raster.DFs, "Data/All.Raster.DFs.rda")

### ======================================================================== ###
### ------------------------ PROCESSING ON ALL DATA ------------------------ ###
### ======================================================================== ###


All.Raster.DFs <- readRDS("Data/All.Raster.DFs.rda")

IVs.DF <- CondenseRasterDFs(All.Raster.DFs)
IVs.DF[100:105,] # check it looks all good

# check which rows contain NA's
NAs <- which(is.na(IVs.DF), arr.ind = T)
NA.rows <- unique(NAs[,1])
# then remove these rows:
IVs.DF <- IVs.DF[-NA.rows, ]

rownames(IVs.DF) <- 1:nrow(IVs.DF)

saveRDS(IVs.DF, "Data/Predictors.DF.rda")
