defineNativeGridVIC <- function(nativeResKm = 1. ){
  #
  # Define a basic native grid at some relevant resolution for the semi-processed underlying data 
  
  # RAW -> Semi-processed Native Grid -> HUC specific input files
  #
  options(stringsAsFactors = F)
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source("/Volumes/MiniPro/Climatics/Code/Master/utilities.R")
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,rgdal,raster,ggplot2)

  # intial parameters
  nativeResDeg <- round(nativeResKm/kmToDeg,3)
  
  # define the native grid based on the US watersheds spatial layer
  #
  # https://viewer.nationalmap.gov/basic/?basemap=b1&category=nhd&title=NHD%20View#startUp
  #
  huc17 <- spTransform(readOGR(dsn=paste0(climaticsDataDir,"Watersheds_Shapefile/HUC17ShapeWithCanda/WBDHU2.shp"), projVals)
  #huc17Simp <- rgeos::gSimplify(huc17,tol=.002,topologyPreserve = T)
  hucSpdf = spTransform(readOGR(dsn="/Volumes/Minipro/Climatics/Data/HDB/wbdhu2_a_us_september2016.gdb"), projVals)
  #hucSimp <- rgeos::gSimplify(hucSpdf,tol=.002,topologyPreserve = T)
  sEdge <- 23; nEdge <- 52; eEdge <- -66; wEdge <- -125
  hucSpdf2 <- crop(hucSpdf, extent(wEdge,eEdge,sEdge,nEdge))
  combined <- union(huc17,hucSpdf2)
  combinedSimp <- rgeos::gSimplify(combined,tol=.002,topologyPreserve = T)
  plot(combinedSimp)
  nativeExt <- extent( combinedSimp)
  westEdge <- (floor(extent(nativeExt)[1] * 1e3)) / 1e3
  southEdge <- (floor(extent(nativeExt)[3] * 1e3)) / 1e3
  xDim <- (round((extent(nativeExt)[2] - westEdge) / nativeResDeg) * 10) / 10
  yDim <- (round((extent(nativeExt)[4] - southEdge) / nativeResDeg) * 10) / 10
  eastEdge <- westEdge + nativeResDeg * xDim
  northEdge <- southEdge + nativeResDeg * yDim
  nCells <- xDim * yDim
  nativeGrid <- raster(xmn=westEdge,xmx=eastEdge, ymn=southEdge, ymx=northEdge, nrows=yDim, ncols= xDim, crs = projVals)
  nativeGrid <- setValues(nativeGrid,rep(1,ncell(nativeGrid)))
  nativeArea <- area(nativeGrid)
  nativeLon <- coordinates(nativeArea)[,1]
  nativeLat <- coordinates(nativeArea)[,2] 
  nativeMask <- mask(nativeGrid,combinedSimp)
  combinedSimp <- crop(combinedSimp,nativeExt)
  plot(nativeMask, col="blue")
  plot(combinedSimp,add=T)
  outfile <- paste0(climaticsDataDir,"RDS/VICNativeGridUSA.raster.rds")
  saveRDS(nativeMask, file = outfile)
  outfile <- paste0(climaticsDataDir,"RDS/VICNativeExtentUSA.shape.rds")
  saveRDS(combinedSimp, file = outfile)
}
defineNativeGridVIC()