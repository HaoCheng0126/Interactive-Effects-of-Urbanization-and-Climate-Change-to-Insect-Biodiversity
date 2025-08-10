##%######################################################%##
#                                                          #
####      Organise Green Space info for sites       ####
#                                                          #
##%######################################################%##

# This script calculates the proportion of natural habitats surrounding each site 
# in PREDICTS subset at multiple spatial scales: 1, 3, 5, and 10 km.

# directories
dataDir <- "0_data/"
predictsDir <- "4_PREDICTSMatchClimateIndex/"
outDir <- "5_PREDICTSMatchPropNatHab/"

if(!dir.exists(outDir)) dir.create(outDir)

# sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)

# load libraries
library(raster)
library(sp)
library(predictsFunctions)
library(snow)

# required crs
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# read in the predicts data with climate info
predicts_sp <- readRDS(paste0(predictsDir,"PREDICTSSitesWithClimateData_update.rds"))

# Load the GLC10 raster layer
glc10 <- raster(paste0(dataDir, "FROM_GLC10.tif"), crs=wgs84)

# Create binary rasters for different land use types
forest <- (glc10 == 20) * 20    
grassland <- (glc10 == 30) * 30   
shrubland <- (glc10 == 40) * 40     
wetland <- (glc10 == 60) * 60        

# Set values to NA for other land types
forest[forest == 0] <- NA
grassland[grassland == 0] <- NA
shrubland[shrubland == 0] <- NA
wetland[wetland == 0] <- NA

# Save each land cover type as a TIF file
writeRaster(forest, filename = paste0(outDir, "forest_glc10.tif"), format="GTiff", overwrite=TRUE)
writeRaster(grassland, filename = paste0(outDir, "grassland_glc10.tif"), format="GTiff", overwrite=TRUE)
writeRaster(shrubland, filename = paste0(outDir, "shrubland_glc10.tif"), format="GTiff", overwrite=TRUE)
writeRaster(wetland, filename = paste0(outDir, "wetland_glc10.tif"), format="GTiff", overwrite=TRUE)

# load the GLC10 raster layers for forest, grassland, shrubland, and wetland
forest <- raster(paste0(dataDir, "forest_glc10.tif"), crs=wgs84)
grassland <- raster(paste0(dataDir, "grassland_glc10.tif"), crs=wgs84)
shrubland <- raster(paste0(dataDir, "shrubland_glc10.tif"), crs=wgs84)
wetland <- raster(paste0(dataDir, "wetland_glc10.tif"), crs=wgs84)

# run in parallel
nCores <- parallel::detectCores()

st1 <- Sys.time()

cl <- snow::makeCluster(nCores-1)

# export data to cores
snow::clusterExport(
  cl = cl,
  list = c('predicts_sp', 'forest', 'grassland', 'shrubland', 'wetland', 'buffer', 'crop', 'cellStats'),
  envir = environment())

# for each site, get the proportion of natural habitats surrounding the site
GreSpace <- data.frame(t(parSapply(cl = cl, X = (1:nrow(predicts_sp)), FUN = function(i){
  cat(paste0("Processing site ", i, " of ", nrow(predicts_sp), "\r"))
  
  # Define buffers
  buffers <- c(1000, 3000, 5000, 10000)
  results <- numeric(length(buffers) * 4)  # 4 habitat types
  
  for (j in seq_along(buffers)) {
    buff <- buffer(predicts_sp[i, ], width=buffers[j])
    
    # Crop the habitat rasters
    forest_crop <- crop(forest, buff)
    grassland_crop <- crop(grassland, buff)
    shrubland_crop <- crop(shrubland, buff)
    wetland_crop <- crop(wetland, buff)
    
    # Calculate the mean proportion of each habitat type
    results[(j-1)*4 + 1] <- cellStats(forest_crop, stat="mean", na.rm=TRUE)
    results[(j-1)*4 + 2] <- cellStats(grassland_crop, stat="mean", na.rm=TRUE)
    results[(j-1)*4 + 3] <- cellStats(shrubland_crop, stat="mean", na.rm=TRUE)
    results[(j-1)*4 + 4] <- cellStats(wetland_crop, stat="mean", na.rm=TRUE)
  }
  
  return(results)
})))

snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # Time difference

# Assign results back to the predicts_sp data frame
predicts_sp$Forest_1000 <- GreSpace$V1
predicts_sp$Grassland_1000 <- GreSpace$V2
predicts_sp$Shrubland_1000 <- GreSpace$V3
predicts_sp$Wetland_1000 <- GreSpace$V4

predicts_sp$Forest_3000 <- GreSpace$V5
predicts_sp$Grassland_3000 <- GreSpace$V6
predicts_sp$Shrubland_3000 <- GreSpace$V7
predicts_sp$Wetland_3000 <- GreSpace$V8

predicts_sp$Forest_5000 <- GreSpace$V9
predicts_sp$Grassland_5000 <- GreSpace$V10
predicts_sp$Shrubland_5000 <- GreSpace$V11
predicts_sp$Wetland_5000 <- GreSpace$V12

predicts_sp$Forest_10000 <- GreSpace$V13
predicts_sp$Grassland_10000 <- GreSpace$V14
predicts_sp$Shrubland_10000 <- GreSpace$V15
predicts_sp$Wetland_10000 <- GreSpace$V16

saveRDS(object = predicts_sp, file = paste0(outDir,"PREDICTSSitesWithClimateAndNatHab.rds"))

t.end <- Sys.time()

print(round(t.end - t.start, 0))

sink()
