createVICForcingRVic <- function(vicGrid = 0, yr = "1980",
                             thisFilename = "Test"){
  
  # Read in a VIC domain netcdf file
  # Use Merra2 data to create the met forcing file
  #
  # k. westrick 1/18/2019
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICForcingMerra2.R"))
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,raster,VICmodel)
  metList <- createVicForcingMerra2(vicGrid = vicGrid,yr = yr)
  for(x in 1:9){
    print(metList[x])
    print(range(getValues(metList[[x]])))
  }
  precTotBrickMmHr <- metList[[1]]
  precBiasCorrBrickMmHr <- metList[[2]]
  tasBrickC <- metList[[3]]
  pressureBrickKPa <- metList[[4]]
  windBrick <- metList[[5]]
  vpBrickKPa <- metList[[6]]
  swBrick <- metList[[7]]
  lwBrick <- metList[[8]]
  range(getValues(lwBrick))
  vioplot::vioplot(as.vector(getValues(lwBrick)))
  hist(getValues(lwBrick))
  merraTerHt <- metList[[9]]
  #
  # read in pararameters file to get the "true" elevation (for lapse rate correction)
  
  merraElevRaster <- resample(merraTerHt,vicGrid) * vicGrid

  # transfer information over to this vic grid
  #
  tasRaster <- resample(tasBrickC,vicGrid) #[[1:240]]
  precTotRaster <- resample(precTotBrickMmHr,vicGrid) #[[1:240]]
  precTotBCRaster <- resample(precBiasCorrBrickMmHr,vicGrid) #[[1:240]]
  pressureRaster <- resample(pressureBrickKPa,vicGrid) #[[1:240]]
  windRaster <- resample(windBrick,vicGrid) #[[1:240]]
  vpRaster <- resample(vpBrickKPa,vicGrid) #[[1:240]]
  swRaster <- resample(swBrick,vicGrid) #[[1:240]]
  lwRaster <- resample(lwBrick,vicGrid) #[[1:240]]
  #colList <- names(STEHE$forcing)
  forcing <- list(
    PREC = matrix(precTotBCRaster,nc=ncell(precTotBCRaster),byrow = T),
    TEMP = matrix(tasRaster,nc=ncell(precTotBCRaster),byrow = T),
    SW = matrix(swRaster,nc=ncell(precTotBCRaster),byrow = T),
    LW = matrix(lwRaster,nc=ncell(precTotBCRaster),byrow = T),
    PRESS = matrix(pressureRaster,nc=ncell(precTotBCRaster),byrow = T),
    VP = matrix(vpRaster,nc=ncell(precTotBCRaster),byrow = T),
    WIND = matrix(windRaster,nc=ncell(precTotBCRaster),byrow = T)
  )
  return(forcing)
}