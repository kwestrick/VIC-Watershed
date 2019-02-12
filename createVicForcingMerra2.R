createVicForcingMerra2 <- function(vicGrid=0, yr = "1980"){
  # 
  # Creates the met forcings netcdf file for VIC
  #
  # Input: vicGrid: a raster mask of the area to create the raster bricks
  # used to create the netcdf. Note that the brick contains the entire extent,
  # and is not trimmed (masked) to the watershed
  #
  # Output: A list of nine fields that can be used to create the netcdf
  # forcing input for VIC
  #
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,raster,humidity)
  # for the details on the creation of the merra2DescriptorsTable.csv see the Vic Evernotes for 19 Sep 2018
  # or use the PDF: https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf to define the fields required
  merra2Descriptors <- read.table(paste0(climaticsCodeDir,"/Master/merra2DescriptorsTable.csv"), sep=",",header = T)
  subsetListDir <- "/Volumes/MiniPro/Climatics/Data/MERRA2/"
  yrDir <- paste0(subsetListDir,"/",yr)
  fileList <- list.files(yrDir)
  dtg <- matrix(unlist(strsplit(fileList,"\\.")),nc=5,byrow = T)[,3]
  uniqueCollections <- unique(merra2Descriptors$merra2_download_fileExt)
  lapply(uniqueCollections, print)
  #
  # expand the VIC grid to include the gridpoints just outside the extent
  #
  infile <- paste0(climaticsDataDir,"MERRA2/MERRA2_101.const_2d_asm_Nx.00000000.nc4.nc")
  foo <- coordinates(raster(infile, varname = "PHIS",levels=c(1)))
  lons <- unique(foo[,1])
  lats <- unique(foo[,2])
  extent(vicGrid)
  expXmin <- max(lons[lons < extent(vicGrid)[1]]) - .625
  expXmax <- min(lons[lons > extent(vicGrid)[2]]) + .625
  expYmin <- max(lats[lats < extent(vicGrid)[3]]) - .5
  expYmax <- min(lats[lats > extent(vicGrid)[4]]) + .5
  trimExt <- extent(expXmin,expXmax,expYmin,expYmax)
  # 
  # get terrain height to apply various lapse information to

  cat("Processing model surface height",fill=T)
  infile <- paste0(climaticsDataDir,"MERRA2/MERRA2_101.const_2d_asm_Nx.00000000.nc4.nc")
  varNames <- names(nc_open(infile)$var)
  terHt <- (crop(raster(infile, varname = "PHIS",levels=c(1)),trimExt)) / 9.80665
  plot(terHt)
  plot(vicGrid, add=T)
  cat("Range of values:",range(getValues(terHt)), fill =T)
  cat("----------------------------------",fill=T)
  mos <- stringr::str_sub(paste0('0',seq(1,12)),-2)
# for (month in mos){
  #cat(paste0("Creating VIC forcing for ",monthList[as.numeric(month)]," ",yr, " based on Merra2 data"),fill=T)
  #subList <- fileList[substr(dtg,5,6) == month]
  #outTitle <- paste0("VIC forcing for ",monthList[as.numeric(month)]," ",yr, " based on Merra2 data")
  cat(paste0("Creating VIC forcing for ",yr, " based on Merra2 data"),fill=T)
  subList <- fileList
  outTitle <- paste0("VIC forcing for ",yr, " based on Merra2 data")
  outHistory <- paste0("Created: ",Sys.time())
  outAuthor <- "Kenneth Westrick, Climatics"
  #
  # ------------- Temperature, surface pressure, wind, and vapor pressure
  #
  thisColl <- uniqueCollections[2] 
  thisList <- subList[grepl(thisColl,subList)]
  thisList <- paste0(yrDir,"/",thisList)
  varNames <- names(nc_open(thisList[1])$var)
  #
  # air temperature (VIC: AIR_TEMP; Average air temperature; C)
  cat("Processing air temperature",fill=T)
  tasBrickK <- crop(do.call(addLayer, lapply(thisList,function(x){
    return(brick(x, varname = "T2M",levels=c(1:24)))
  })),trimExt)
  tasBrickC <- tasBrickK - 273.15
  cat("Range of values:",range(getValues(tasBrickC)), fill =T)
  cat("----------------------------------",fill=T)
  #
  # surface pressure (VIC: PRESSURE; Atmospheric pressure; kPa)
  cat("Processing pressure",fill=T)
  pressureBrickPa <- crop(do.call(addLayer, lapply(thisList,function(x){
    return(brick(x, varname = "PS",levels=c(1:24)))
  })),trimExt)
  pressureBrickKPa <- pressureBrickPa / 1000. # unit conversion from pascals kilo-pascals
  cat("Range of values:",range(getValues(pressureBrickKPa)), fill =T)
  cat("----------------------------------",fill=T)
  #
  # windspeed based on u and v 10m winds (VIC: WIND; wind speed; m/s)
  cat("Processing wind",fill=T)
  u10mBrick <- crop(do.call(addLayer, lapply(thisList,function(x){
    return(brick(x, varname = "U10M",levels=c(1:24)))
  })),trimExt)
  v10mBrick <- crop(do.call(addLayer, lapply(thisList,function(x){
    return(brick(x, varname = "V10M",levels=c(1:24)))
  })),trimExt)
  windBrick <- sqrt(u10mBrick^2 + v10mBrick^2)
  cat("Range of values:",range(getValues(windBrick)), fill =T)
  cat("----------------------------------",fill=T)
  #
  # vapor pressure (VIC: WIND; wind speed; m/s)
  cat("Processing vapor pressure",fill=T)
  qvBrick <- crop(do.call(addLayer, lapply(thisList,function(x){
    return(brick(x, varname = "QV2M",levels=c(1:24)))
  })), trimExt)
  #cat("calculating vapor pressure from specific humidity...", fill=T)
  rhBrick <- overlay(qvBrick,tasBrickK,pressureBrickPa,fun=SH2RH)
  rhBrick[rhBrick > 100.] <-100.
  esBrick <- calc(tasBrickK,SVP)
  vpBrickKPa <- overlay(rhBrick,esBrick,fun=WVP2) / 1000.
  cat("Range of values:",range(getValues(qvBrick)), fill =T)
  cat("----------------------------------",fill=T)
  #
  # ------------- Shortwave and Longwave-------------------------------
  #
  thisColl <- uniqueCollections[6] 
  thisList <- subList[grepl(thisColl,subList)]
  thisList <- paste0(yrDir,"/",thisList)
  varNames <- names(nc_open(thisList[1])$var)
  #
  # shortwave (VIC: SWDOWN; Incoming shortwave radiation; W/m^2)
  cat("Processing shortwave radiation",fill=T)
  swBrick <- crop(do.call(addLayer, lapply(thisList,function(x){
    return(brick(x, varname = "SWGDN",levels=c(1:24)))
  })),trimExt)
  cat("Range of values:",range(getValues(swBrick)), fill =T)
  cat("----------------------------------",fill=T)
  #
  # longwave (VIC: LWDOWN; Incoming longwave radiation; W/m^2)
  cat("Processing longwave radiation",fill=T)
  lwBrick <- crop(do.call(addLayer, lapply(thisList,function(x){
    return(brick(x, varname = "LWGAB",levels=c(1:24)))
  })),trimExt)
  cat("Range of values:",range(getValues(lwBrick)), fill =T)
  cat("----------------------------------",fill=T)
  # ------------- Precipitation----------------------------------
  #
  thisColl <- uniqueCollections[5] 
  thisList <- subList[grepl(thisColl,subList)]
  thisList <- paste0(yrDir,"/",thisList)
  varNames <- names(nc_open(thisList[1])$var)
  #
  # totalPrecip (VIC: PREC;	Total precipitation; (rain and snow);	mm)
  cat("Processing precipitation",fill=T)
  precTotBrickKgM2S <- crop(do.call(addLayer, lapply(thisList,function(x){
    return(brick(x, varname = "PRECTOT",levels=c(1:24)))
  })),trimExt)
  precTotBrickMmHr <- precTotBrickKgM2S * 60 * 60
  cat("Range of values:",range(getValues(precTotBrickMmHr)), fill =T)
  cat("----------------------------------",fill=T)
  #
  # totalPrecipCorr (VIC: PREC;	Total precipitation; (rain and snow);	mm)
  cat("Processing bias corrected precipitation",fill=T)
  precBiasCorrBrickKgM2S <- crop(do.call(addLayer, lapply(thisList,function(x){
    return(brick(x, varname = "PRECTOTCORR",levels=c(1:24)))
  })),trimExt)
  precBiasCorrBrickMmHr <- precBiasCorrBrickKgM2S * 60 * 60 # (1 kg / m^2) = 1 mm
  cat("Range of values:",range(getValues(precBiasCorrBrickMmHr)), fill =T)
  cat("----------------------------------",fill=T)
  #
  outList <- list(precTotBrickMmHr,precBiasCorrBrickMmHr,
                  tasBrickC,pressureBrickKPa,
                  windBrick,vpBrickKPa,
                  swBrick,lwBrick, terHt)
  names(outList) <- c("PREC","PRECBIASCORR","AIR_TEMP","PRESSURE",
                      "WIND","VP","SWDOWN","LWDOWN","TERHT")
  return(outList)
}
