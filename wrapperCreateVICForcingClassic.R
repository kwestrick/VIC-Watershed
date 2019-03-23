wrapperCreateVICForcingClassic <- function(yrStart = 0, yrEnd = 0,
                                           dataDir = ""
                                           ){
  
  # A wrapper to quickly call the routine to create the forcing files 
  # for the classic version of VIC. The program assumes that the program
  # "driverVICInputsClassic.R" has been run and that the various 
  # geographic files exist in the appropriate directories
  #
  #
  options(stringsAsFactors = F)
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,rgdal,raster,openxlsx,ggmap,rgeos)
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  stackFile <- list.files(dataDir,pattern = "stack_allRasters_",full.names = T)[1]
  stackData <- stack(gsub(".grd","",stackFile))
  names(stackData)
  watershedFile <- list.files(dataDir,pattern = "watershed_shape_",full.names = T)
  thisWatershedBdry <- readRDS(watershedFile)
  yrs <- as.character(seq(yrStart,yrEnd,1))
  for(yr in yrs){
    createVICForcingClassic(vicGrid = vicGrid, 
                            yr = yr,
                            saveDir = saveDir)
    gc()
  }
}
wrapperCreateVICForcingClassic(yrStart = 2014, yrEnd = 2019,
  dataDir = "/Volumes/MiniPro/Climatics/Operational/Data/KeenleysideDam_16km/" 
 
  