createMetForcing <- function(maskRaster = 0, yr = "1980"){
  
  # Read in a VIC domain netcdf file
  # Use Merra2 data to create the met forcing file
  #
  # k. westrick 1/18/2019
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  rotate <- function(x) t(apply(x, 2, rev))
  missingValue <- -999.99
  myFlip <- function(s=0){
    tmp <- flip(t(s),direction='x')
    tmp[is.na(tmp)] <- missingValue
    return(tmp)
  }
  
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,raster)
  metList <- createVicForcingMerra2(vicGrid = maskRaster,yr = yr)
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
  
  merraElevRaster <- resample(merraTerHt,maskRaster) * maskRaster

  # transfer information over to this vic grid
  #
  tasRaster <- (resample(tasBrickC,maskRaster) * maskRaster)#[[1:240]]
  precTotRaster <- (resample(precTotBrickMmHr,maskRaster) * maskRaster)#[[1:240]]
  precTotBCRaster <- (resample(precBiasCorrBrickMmHr,maskRaster) * maskRaster)#[[1:240]]
  pressureRaster <- (resample(pressureBrickKPa,maskRaster) * maskRaster)#[[1:240]]
  windRaster <- (resample(windBrick,maskRaster) * maskRaster)#[[1:240]]
  vpRaster <- (resample(vpBrickKPa,maskRaster) * maskRaster)#[[1:240]]
  swRaster <- (resample(swBrick,maskRaster) * maskRaster)#[[1:240]]
  lwRaster <- (resample(lwBrick,maskRaster) * maskRaster)#[[1:240]]

  #################
  #
  # CREATE THE NETCDF FILES
  #
  #################
  
  #
  vicLon <- coordinates(maskRaster)[,1]
  vicLat <- coordinates(maskRaster)[,2]
  #
  # dimensions
  #
  timeUnitsDescriptor <- paste0("hours since ",yr,"-01-01")
  #timeUnitsDescriptor <- paste0("hours since 1949-01-01")
  times <- as.integer(seq(0,(dim(tasRaster)[3])-1,1))
  dimTime <- ncdim_def("time",timeUnitsDescriptor,times,calendar = "proleptic_gregorian")
  dimLon <- ncdim_def("lon","degrees_east",as.double(unique(vicLon)),longname = "longitude of grid cell center") 
  dimLat <- ncdim_def("lat","degrees_north",as.double(sort(unique(vicLat))),longname = "latitude of grid cell center")
  #
  # thisDTG <- substr(as.character(gsub("-","",Sys.time())),1,8)
  thisDTG <- "test"
  var_time <- ncvar_def("time","",list(dimTime),missingValue,"",prec="integer")
  var_lat <- ncvar_def("lat","degrees_north",dimLat,missingValue,"latitude of grid cell center",prec="single")
  var_lon <- ncvar_def("lon","degrees_east",dimLon,missingValue,"longitude of grid cell center",prec="single")
  var_mask <- ncvar_def("mask","",list(dimLon,dimLat),missingValue,longname = "domain mask", prec = "double")
  var_prcp <- ncvar_def("prcp","mm/step",list(dimLon,dimLat,dimTime),missingValue,longname = "", prec = "float")
  var_tas <- ncvar_def("tas","",list(dimLon,dimLat,dimTime),missingValue,longname = "", prec = "float")
  var_dswrf <- ncvar_def("dswrf","",list(dimLon,dimLat,dimTime),missingValue,longname = "", prec = "float")
  var_dlwrf <- ncvar_def("dlwrf","",list(dimLon,dimLat,dimTime),missingValue,longname = "", prec = "float")
  var_pres <- ncvar_def("pres","",list(dimLon,dimLat,dimTime),missingValue,longname = "", prec = "float")
  var_vp <- ncvar_def("vp","",list(dimLon,dimLat,dimTime),missingValue,longname = "", prec = "float")
  var_wind <- ncvar_def("wind","",list(dimLon,dimLat,dimTime),missingValue,longname = "", prec = "float")

 
  metNcdfFile <- paste0("/Volumes/MiniPro/Climatics/Code/VIC/vic/drivers/image/forcings/myStehekinForcings.1980.nc")
  metNcdf <- nc_create(metNcdfFile,list(var_mask,
                                         var_tas,var_prcp,var_dswrf,
                                         var_dlwrf,var_pres,var_vp,
                                         var_wind),verbose = F,force_v4 = TRUE)
  ncvar_put(metNcdf,var_tas,as.array(myFlip(tasRaster)))
  ncvar_put(metNcdf,var_mask,as.matrix(myFlip(maskRaster)))
  ncvar_put(metNcdf,var_prcp,as.array(myFlip(precTotRaster)))
  ncvar_put(metNcdf,var_dswrf,as.array(myFlip(swRaster)))
  ncvar_put(metNcdf,var_dlwrf,as.array(myFlip(lwRaster)))
  ncvar_put(metNcdf,var_pres,as.array(myFlip(pressureRaster)))
  ncvar_put(metNcdf,var_vp,as.array(myFlip(vpRaster)))
  ncvar_put(metNcdf,var_wind,as.array(myFlip(windRaster)))
  nc_close(metNcdf)
  
}