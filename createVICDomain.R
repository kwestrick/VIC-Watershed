createVICDomain <- function(vicBdry = 0, 
                            resolutionKm = 4.,
                            thisFilename = "Test"){
  
  missingValue <- -999.99
  myFlip <- function(s=0,missingValue = NULL){
    tmp <- flip(t(s),direction='x')
    tmp[is.na(tmp)] <- missingValue
    return(tmp)
  }
  # define VIC grid
  nativeMask <- readRDS(file = paste0(climaticsDataDir,"RDS/VICNativeGridUSA.raster.rds"))
  initialExtent <- extent(vicBdry)
  t1 <- cellsFromExtent(nativeMask, initialExtent, expand = T)
  initialExtentSnap <- extentFromCells(nativeMask,t1)
  origMask <- crop(mask(nativeMask, vicBdry),initialExtentSnap, snap = "out")
  if (resolutionKm > 1 && resolutionKm %% 1 == 0) vicGrid <- aggregate(origMask,fact = resolutionKm)
  origMask[is.nan(origMask)] <- NA

  plot(vicGrid,colNA="blue")
  plot(vicBdry, add=T)
  
  # create run_cell mask
  runCell <- vicGrid
  runCell[!is.na(runCell)] <- 1
  runCell[is.na(runCell)] <- 0
  plot(runCell)
  
  # gridcell; sequence of numbers for the grid cells
  gridcell <- runCell
  gridcell[gridcell == 1] <- seq(1,ncell(gridcell[gridcell == 1]))
  plot(gridcell)
  plot(vicBdry, add=T, border="grey40")
  
  # area
  gridArea <- area(vicGrid) * 1e6
  totalVicArea <- area(origMask) * 1e6
  vicArea <- aggregate(totalVicArea * origMask, fun= sum, fact = resolutionKm)
  vicArea[is.nan(vicArea)] <- NA
  totalVicArea <- aggregate(totalVicArea, fun= sum, fact = resolutionKm)
  totalVicArea[is.nan(totalVicArea)] <- NA
  as.matrix(vicArea) / 1e6
  as.matrix(totalVicArea) / 1e6
  plot(vicArea,colNA="black")
  plot(vicBdry, add=T, border="grey40")
  vicArea <- round(vicArea,0)
  
  # fraction of grid cell active
  vicFrac <- signif(vicArea / totalVicArea,5)
  as.matrix(vicFrac)
  plot(vicFrac,colNA="black")
  plot(vicBdry, add=T, border="grey40")
  
# CREATE THE NETCDF FILES

  vicLon <- coordinates(vicArea)[,1]
  lonRaster <- rasterFromXYZ(cbind(coordinates(vicArea),coordinates(vicArea)[,1]))
  vicLat <- coordinates(vicArea)[,2]
  latRaster <- rasterFromXYZ(cbind(coordinates(vicArea),coordinates(vicArea)[,2]))
  
  #
  # dimensions
  #
  dimLon <- ncdim_def("lon","degrees_east",as.double(unique(vicLon)),longname = "longitude of grid cell center") 
  dimLat <- ncdim_def("lat","degrees_north",as.double(sort(unique(vicLat))),longname = "latitude of grid cell center")
  dimLayers <- ncdim_def("nlayer","",1:3,longname = "soil layer")
  
  ##########
  #
  # create the domain file
  # 
  ##########
  
  var_lat <- ncvar_def("lat","degrees_north",dimLat,missingValue,"latitude of grid cell center",prec="single")
  var_lon <- ncvar_def("lon","degrees_east",dimLon,missingValue,"longitude of grid cell center",prec="single")
  var_mask <- ncvar_def("mask","",list(dimLon,dimLat),missingValue,longname = "domain mask", prec = "integer")
  var_frac <- ncvar_def("frac","unitless",list(dimLon,dimLat),missingValue,longname = "fraction of grid cell that is active", prec = "double")
  var_area <- ncvar_def("area","m2",list(dimLon,dimLat),missingValue,longname = "area of grid cell", prec = "double")
  
  domainsNcdfFile <- paste0(climaticsCodeDir,"/VIC/vic/drivers/image/parameters/domain.",thisFilename,".nc")
  cat("writing ",domainsNcdfFile,fill=T)
  domainsNcdf <- nc_create(domainsNcdfFile,list(var_mask,var_frac,var_area),verbose = F)
  ncatt_put(domainsNcdf,0,"_missingValue",missingValue)
  as.matrix(runCell)
  runCell[runCell == 0] <- missingValue
  ncvar_put(domainsNcdf,var_mask,as.matrix(myFlip(runCell)))
  as.matrix(vicFrac)
  vicFrac[is.na(vicFrac)] <- missingValue
  ncvar_put(domainsNcdf,var_frac,as.matrix(myFlip(vicFrac)))
  gridArea[is.na(gridArea)] <- missingValue
  ncvar_put(domainsNcdf,var_area,as.matrix(myFlip(gridArea)))
  nc_close(domainsNcdf)
  return(vicGrid)
}

