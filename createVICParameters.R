createVICParameters <- function(vicBdry = 0, vicGrid = 0, thisFilename = "Test",
                                resolutionKm = NULL){
  # k. westrick 2/5/2019
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICForcingMerra2.R"))
  source(paste0(climaticsCodeDir,"Master/createLandCoverImage.R"))
  if(!requireNamespace("devtools")) install.packages("devtools")
  if (!require("pacman")) install.packages("pacman")
  
  pacman::p_load(ncdf4,raster)
  missingValue <- -999.99
  myFlip <- function(s=0){
    tmp <- flip(t(s),direction='x')
    tmp[is.na(tmp)] <- missingValue
    return(tmp)
  }
  vicGrid <- readRDS("./stehekin14km.vicGrid.rds")
  vicGrid[is.nan(vicGrid)] <- NA
  thisFilename <- "stehekin.20190203"
  # define and check for the availability of the required input files
  
  defineCheckInput <- function(){
    soilParams <- data.frame(
      # run_cell is set in 
      infilt = T,
      Ds = T,
      Dsmax = T,
      Ws = T,
      c = T,
      expt = T,
      Ksat = T,
      phi_s = T,
      init_moist = T,
      elev = T,
      depth = T,
      avg_T = T,
      dp = T,
      bubble = T,
      quartz = T,
      bulk_density = T,
      soil_density = T,
      organic = F,
      bulk_dens_org = F,
      soil_dens_org = F,
      off_gmt = T,
      Wcr_FRACT = T,
      Wpwp_FRACT = T,
      rough = T,
      snow_rough = T,
      annual_prec= T,
      resid_moist = T,
      fs_active = T,
      frost_slope = F,
      max_snow_distrib_slope = F,
      July_Tavg = T
    )
    soilParamsExluded <- names(soilParams[,which(soilParams == F)])
    soilParamsIncluded <- soilParams[,which(soilParams == T)]
    
    
    l1 <- list.files(paste0(climaticsDataDir,"VicInputGrids/"),pattern = "^Vic.soil.*.grd")
    soilVicFileList <- matrix(unlist(strsplit(l1,"\\.")),ncol = 4, byrow = T)[,3]
    missingMask <- !names(soilParamsIncluded) %in% soilVicFileList 
    missingSoilParams <- names(soilParamsIncluded[missingMask])
    if(length(missingSoilParams) == 0) cat("All base files are available\n") else {
      for(i in 1:length(missingSoilParams)){
        cat(paste0("missing ",missingSoilParams[i],"\n",fill=T))
      }
      stop("run processSoilForVic.R")
    }
    return(soilParams)
  }
  soilParams <- defineCheckInput()
  
  # we need 3 input files from the domain file
  domainsNcdfFile <- paste0(climaticsCodeDir,"/VIC/vic/drivers/image/parameters/domain.",thisFilename,".nc")
  domainPtr <- nc_open(domainsNcdfFile)
  as.matrix(names(domainPtr$var))
  runCell <- raster(domainsNcdfFile,varname = "mask")
  runCell[!is.na(runCell)] <- 1
  runCell[is.na(runCell)] <- 0
  plot(runCell)

  ############
  #
  # name and open the parameters file
  #
  ############
  vicLon <- coordinates(vicGrid)[,1]
  lonRaster <- rasterFromXYZ(cbind(coordinates(vicGrid),coordinates(vicGrid)[,1]))
  vicLat <- coordinates(vicGrid)[,2]
  latRaster <- rasterFromXYZ(cbind(coordinates(vicGrid),coordinates(vicGrid)[,2]))
  #
  # dimensions
  #
  dimLon <- ncdim_def("lon","degrees_east",as.double(unique(vicLon)),longname = "longitude of grid cell center") 
  dimLat <- ncdim_def("lat","degrees_north",as.double(sort(unique(vicLat))),longname = "latitude of grid cell center")
  dimLayers <- ncdim_def("nlayer","",1:3,longname = "soil layer")
  
  
  # SOIL
  var_layer <- ncvar_def("layer","",list(dimLayers),missingValue,longname = "soil_layer",prec = "integer")
  var_run_cell <- ncvar_def("run_cell","unitless",list(dimLon,dimLat), missingValue,"run_cell",prec = "integer")
  var_gridCell <- ncvar_def("gridCell","N/A",list(dimLon,dimLat),missingValue,prec = "integer", longname = "grid_cell")
  var_lats <- ncvar_def("lats","degrees",list(dimLon,dimLat),missingValue,prec = "double", longname = "Latitude of grid cell")
  var_lons <- ncvar_def("lons","degrees",list(dimLon,dimLat),missingValue,prec = "double", longname = "Longitude of grid cell")
  var_infilt <- ncvar_def("infilt", "mm/day",list(dimLon,dimLat),missingValue,prec = "double", longname = "infilt")
  var_Ds <- ncvar_def("Ds","fraction",list(dimLon,dimLat),missingValue,prec = "double", longname = "Ds")
  var_Dsmax <- ncvar_def("Dsmax","mm/day",list(dimLon,dimLat),missingValue,prec = "double", longname = "Dsmax")
  var_Ws <- ncvar_def("Ws","fraction",list(dimLon,dimLat),missingValue,prec = "double", longname = "Ws ")
  var_c <- ncvar_def("c","N/A",list(dimLon,dimLat),missingValue,prec = "double", longname = "c")
  var_expt <- ncvar_def("expt","N/A",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "expt")
  var_Ksat <- ncvar_def("Ksat","mm/day",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "kSat")
  var_phi_s <- ncvar_def("phi_s","mm/mm",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "phi_s")
  var_init_moist <- ncvar_def("init_moist","mm",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "init_moist")
  var_elev <- ncvar_def("elev","",list(dimLon,dimLat),missingValue,prec = "double", longname = "elev")
  var_depth <- ncvar_def("depth","m",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "depth")
  var_avg_T <- ncvar_def("avg_T","C",list(dimLon,dimLat),missingValue,prec = "double", longname = "avg_T")
  var_dp <- ncvar_def("dp","m",list(dimLon,dimLat),missingValue,prec = "double", longname = "dp")
  var_bubble <- ncvar_def("bubble","cm",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "bubble")
  var_quartz <- ncvar_def("quartz","fraction",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "quartz")
  var_bulk_density <- ncvar_def("bulk_density","kg/m3",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "bulk_density")
  var_soil_density <- ncvar_def("soil_density","kg/m3",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "soil_density")
  var_off_gmt <- ncvar_def("off_gmt","hours",list(dimLon,dimLat),missingValue,prec = "double", longname = "off_gmt")
  var_Wcr_FRACT <- ncvar_def("Wcr_FRACT","fraction",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "Wcr_FRACT")
  var_Wpwp_FRACT <- ncvar_def("Wpwp_FRACT","fraction",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "Wpwp_FRACT")
  var_rough <- ncvar_def("rough","m",list(dimLon,dimLat),missingValue,prec = "double", longname = "rough")
  var_snow_rough <- ncvar_def("snow_rough","m",list(dimLon,dimLat),missingValue,prec = "double", longname = "snow_rough")
  var_annual_prec <- ncvar_def("annual_prec","mm",list(dimLon,dimLat),missingValue,prec = "double", longname = "annual_prec")
  var_resid_moist <- ncvar_def("resid_moist","",list(dimLon,dimLat,dimLayers),missingValue,prec = "double", longname = "resid_moist")
  var_fs_active <- ncvar_def("fs_active","binary ",list(dimLon,dimLat),missingValue,prec = "integer", longname = "fs_active")
  var_July_Tavg <- ncvar_def("July_Tavg","C ",list(dimLon,dimLat),missingValue,prec = "double", longname = "July_Tavg")

  paramsNcdfFile <- paste0(climaticsCodeDir,"/VIC/vic/drivers/image/parameters/parameters.",thisFilename,".nc")
  varsForList <- c("var_run_cell","var_gridCell","var_lats","var_lons", paste0("var_",names(soilParams[which(soilParams == T)])))
  varList <- lapply(varsForList, get)

  writeSoil <- T
  if(writeSoil){
    paramsNcdf <- nc_create(paramsNcdfFile,varList,verbose = F)
    
    # process and place the basic fields in the netcdf
    # outCell (0,1 mask); retain runCell with NA's
    
    ncvar_put(paramsNcdf,var_run_cell,as.matrix(myFlip(runCell)))
    ncatt_put(paramsNcdf,var_run_cell,"description","1 = Run Grid Cell, 0 = Do Not Run")
    
    # gridCell (sequental numbering)
    gridCell <- vicGrid
    gridCell[!is.na(gridCell)] <- seq(1,ncell(gridCell[gridCell == 1]))
    plot(gridCell,colNA="black")
    plot(vicBdry,add=T)
    
    ncvar_put(paramsNcdf,var_gridCell,as.matrix(myFlip(gridCell)))
    ncatt_put(paramsNcdf,var_gridCell,"description","Grid cell number")
    
    # lats and lons
    ncvar_put(paramsNcdf,var_lats,as.matrix(myFlip(latRaster)))
    ncvar_put(paramsNcdf,var_lons,as.matrix(myFlip(lonRaster)))
    
    # process and place all the soil fields in the netcdf
    
    myAggregate <- function(inParam="",inFun="mean",inRes=0){
      inFile <- paste0(climaticsDataDir,"VicInputGrids/Vic.soil.",inParam,".grd")
      tmp <- crop(stack(inFile),vicBdry)
      if(inFun == "mean") outGrid <- aggregate(tmp,fact = inRes)
      if(inFun == "max") outGrid <- aggregate(tmp,fun = max, fact = inRes)
      outGrid[is.nan(outGrid)] <- NA
      return(outGrid)
    }
    soilVars <- names(soilParams[which(soilParams == T)])
    for(inParam in soilVars){
      cat(inParam,fill = T)
      varIn <- paste0("var_",inParam)
      assign(inParam, myAggregate(inParam,inRes=resolutionKm))
      tmpGrid <- mask(get(inParam), runCell, maskvalue = 0) * vicGrid
      if(dim(tmpGrid)[3] == 1) plot(tmpGrid,colNA="black",main=inParam) else plot(tmpGrid[[1]],colNA="black",main=inParam)
      plot(vicBdry, add=T, border="grey40")
      tmpGrid[is.na(tmpGrid)] <- missingValue
      if(dim(tmpGrid)[3] == 1) ncvar_put(paramsNcdf,get(varIn),as.matrix(myFlip(tmpGrid))) else {
        ncvar_put(paramsNcdf,get(varIn),as.array(myFlip(tmpGrid))) 
      }
      ncatt_put(paramsNcdf,get(varIn),"description","")
    }
    nc_close(paramsNcdf)
  }
  ##########
  #
  # Snow Parameters
  #
  ##########
  
  # Logic: if the resolution is greater than thresholdResolutionMakeSnowBands create snow bands
  # snow bands are limited to maxSnowBands
  # if the elevation difference is smaller than the elevDiffThreshold go with something smaller
  
  # NOTE: right now the routine does NOT have the ability to subset
  
  thresholdResolutionMakeSnowBands <- 5
  if(resolutionKm > thresholdResolutionMakeSnowBands) {
    # SNOW
    maxSnowBands <- 6
    ElevDiffThreshold <- 280 # this is the elevation difference based on the binning on whether to use maxSnowBands as the number of valid bins, or use something smaller if the 
    # nmbrSnowBands <- min(resolutionKm / 5, maxSnowBands)
    #
    inParam <- "elev"
    inFile <- paste0(climaticsDataDir,"VicInputGrids/Vic.soil.",inParam,".grd")
    elevHiRes <- mask(crop(stack(inFile),vicGrid),vicBdry)
    plot(elevHiRes)
    plot(vicBdry,add=T)
    
    centerCoords <- coordinates(vicGrid)
    maskVals <- getValues(vicGrid)
    coarseDF <- data.frame(lon=centerCoords[,1],lat=centerCoords[,2],val=maskVals)
    coarseDF <- coarseDF[!is.na(coarseDF$val),]
    
    # make the grid reference: https://gis.stackexchange.com/questions/113963/create-a-rectangular-buffer-from-point-shapefile-using-qgis-grass-saga-gis-r
    # keywords: gridbox
    x <- coarseDF[,1]
    y <- coarseDF[,2]
    sp <- SpatialPoints(data.frame(x,y),proj4string = crs(vicGrid))
    topleftCorner <- bbox(sp)[,1]
    columns <- length(unique(x))
    rows <- length(unique(y))
    cellWidth <- res(vicGrid)[1]
    cellHeight <- res(vicGrid)[2]
    grd <- GridTopology(topleftCorner, c(cellWidth,cellHeight), c(columns,rows))
    grdPolys <- as.SpatialPolygons.GridTopology(grd) #rectangles
    crs(grdPolys) <- crs(vicGrid)
    # remove the polys which don't have a point
    grdPolys <- grdPolys[over(sp, grdPolys),]
    nmbrCells <- length(grdPolys)
    retrieveSnowBand <- function(cellNmbr){
      
      # This needs to be fixed, right now when the elevation difference is small it still defaults to nmbrSnowBands
      # this should actually go to something smaller, as we don't want to process the snow bands if they are not needed
      # see the commented out lines in the below "else" statement.
      print(cellNmbr)
      thisCell <- grdPolys[cellNmbr,]
      plot(elevHiRes)
      plot(grdPolys, border="blue",add = T)
      plot(sp, col="black", add=T, pch="+",cex=2)
      plot(vicBdry, add=T, border="grey40")
      color <- c(2,2,3,4,5) 
      color_transparent <- adjustcolor(color, alpha.f = 0.5) 
      plot(thisCell,col=color_transparent,add=T)
      
      thisElev <- mask(crop(elevHiRes,thisCell),thisCell)
      # plot(thisElev)
      # plot(grdPolys, border="blue",add = T)
      # plot(vicBdry, add=T, border="grey40")
      allVals <- getValues(thisElev)
      validVals <- allVals[!is.na(allVals)]
      fractionInCell <- length(validVals) / length(allVals)
      minElev <- min(validVals)
      maxElev <- max(validVals)
      binSize <- (maxElev - minElev) / maxSnowBands
      elevDiffThreshold <- 100
      if(binSize < elevDiffThreshold){
        nmbrBins <- ceiling((maxElev - minElev) / elevDiffThreshold) + 1
        nmbrBins <- min(nmbrBins,maxSnowBands+1)
        histOut <- hist(validVals,breaks=seq(minElev, by = elevDiffThreshold, length.out = nmbrBins), 
                        freq=FALSE,col="orange",
                        main="Histogram",xlab="x",ylab="f(x)",yaxs="i",xaxs="i")
        elevBands <- round(histOut$mids,0)
        areaFrac <- round(histOut$counts / sum(histOut$counts),3)
        residual <- 1 - sum(areaFrac)
        if(residual != 0) areaFrac[which.max(areaFrac)] <- areaFrac[which.max(areaFrac)] + residual # true up for fractions to equal one
        pFactor <- areaFrac
        if(length(areaFrac) < maxSnowBands){ # we have to pad with zeros so the length of the snow bands == maxSnowBands
          areaFrac <- c(areaFrac,rep(0,maxSnowBands-length(areaFrac)))
          elevBands <- c(elevBands,rep(0,maxSnowBands-length(elevBands)))
          pFactor <- c(pFactor,rep(0,maxSnowBands-length(pFactor)))
        }
      } else {
        histOut <- hist(validVals,breaks=seq(min(validVals),max(validVals),l=maxSnowBands+1), 
                        freq=FALSE,col="orange",
                        main="Histogram",xlab="x",ylab="f(x)",yaxs="i",xaxs="i")
        binBracket <- (histOut$breaks[length(histOut$breaks)] - histOut$breaks[1]) / (maxSnowBands - 1)
        elevBands <- round(histOut$mids,0)
        areaFrac <- round(histOut$counts / sum(histOut$counts),3)
        residual <- 1 - sum(areaFrac)
        if(residual != 0) areaFrac[which.max(areaFrac)] <- areaFrac[which.max(areaFrac)] + residual # true up for fractions to equal one
        pFactor <- areaFrac
      }
      cellMeanElev <- mean(getValues(thisElev),na.rm=T)
      outList <- list(cellNmbr, elevBands,areaFrac,pFactor,cellMeanElev)
    }
    snowParamList <- do.call(list, lapply(1:nmbrCells, retrieveSnowBand))
    nmbrCols <- length(unlist(snowParamList[1]))
    snowParamMatrix <- matrix(unlist(snowParamList),nc=nmbrCols,byrow = T)
    elevBands <- snowParamMatrix[,2:(maxSnowBands + 1)]
    areaFract <- snowParamMatrix[,(2 + maxSnowBands):(maxSnowBands*2 + 1)]
    Pfactor <- snowParamMatrix[,(2 + maxSnowBands*2):(maxSnowBands*3 + 1)]
    grdPolysCoords <- gCentroid(grdPolys,byid = T)
    elevation <-rasterFromXYZ(cbind(coordinates(grdPolysCoords),elevBands))
    plot(elevation)
    areaFract <- rasterFromXYZ(cbind(coordinates(grdPolysCoords),areaFract))
    plot(areaFract)
    Pfactor <- rasterFromXYZ(cbind(coordinates(grdPolysCoords),Pfactor))
    
    dimSnowBand <- ncdim_def("snow_band","",1:maxSnowBands, longname = "snow band")
    var_cellnum <- ncvar_def("cellnum","N/A",list(dimLon,dimLat),missingValue,prec = "double", longname = "")
    var_areaFract <- ncvar_def("AreaFract","fraction",list(dimLon,dimLat,dimSnowBand),missingValue,prec = "double")
    var_elevation <- ncvar_def("elevation","m",list(dimLon,dimLat,dimSnowBand),missingValue,prec = "double", longname = "")
    var_Pfactor <- ncvar_def("Pfactor","fraction",list(dimLon,dimLat,dimSnowBand),missingValue,prec = "double", longname = "")
    
    ncid_old <- nc_open(paramsNcdfFile, write=T)
    #ncnew <- nc_create(paramsNcdfFile,ncid_old)
    
    ncid_old <- ncvar_add(ncid_old,var_cellnum)
    ncid_old <- ncvar_add(ncid_old,var_areaFract)
    ncid_old <- ncvar_add(ncid_old,var_elevation)
    ncid_old <- ncvar_add(ncid_old,var_Pfactor)
    
    ncvar_put(ncid_old,"cellnum",as.matrix(myFlip(gridCell)))
    
    areaFract[is.na(areaFract)] <- missingValue
    outArray <- as.array(myFlip(areaFract))
    ncvar_put(ncid_old,"AreaFract",outArray)
    
    elevation[is.na(elevation)] <- missingValue
    outArray <- as.array(myFlip(elevation))
    ncvar_put(ncid_old,"elevation",outArray)
    
    Pfactor[is.na(Pfactor)] <- missingValue
    outArray <- as.array(myFlip(Pfactor))
    ncvar_put(ncid_old,"Pfactor",outArray)
    nc_close(ncid_old)
  }

  
  ##########
  #
  # Vegetation Parameters
  #
  ##########
  
  processVeg <- F
  if(processVeg){
   infile <- paste(rdsDataDir,"globalLandCoverV2.westernHemisphere.RData",sep='')
    landCover <- get(load(file=infile))
    thisLandCover <- crop(landCover,vicGrid)
    createLandCoverImage(inRaster = thisLandCover,inBoundary=vicBdry)
    landCvrLegendData <- read.table(file=paste0(climaticsCodeDir,'Master/legendColorsLabels.csv'),header = T,sep=",",skip = 1)
    
    centerCoords <- coordinates(vicGrid)
    maskVals <- getValues(vicGrid)
    coarseDF <- data.frame(lon=centerCoords[,1],lat=centerCoords[,2],val=maskVals)
    coarseDF <- coarseDF[!is.na(coarseDF$val),]
    
    # make the grid reference: https://gis.stackexchange.com/questions/113963/create-a-rectangular-buffer-from-point-shapefile-using-qgis-grass-saga-gis-r
    # keywords: gridbox
    x <- coarseDF[,1]
    y <- coarseDF[,2]
    sp <- SpatialPoints(data.frame(x,y),proj4string = crs(vicGrid))
    topleftCorner <- bbox(sp)[,1]
    columns <- length(unique(x))
    rows <- length(unique(y))
    cellWidth <- res(vicGrid)[1]
    cellHeight <- res(vicGrid)[2]
    grd <- GridTopology(topleftCorner, c(cellWidth,cellHeight), c(columns,rows))
    grdPolys <- as.SpatialPolygons.GridTopology(grd) #rectangles
    crs(grdPolys) <- crs(vicGrid)
    # remove the polys which don't have a point
    grdPolys <- grdPolys[over(sp, grdPolys),]
    nmbrCells <- length(grdPolys)
    processVeg <- function(cellNmbr){
      print(cellNmbr)
      thisCell <- grdPolys[cellNmbr,]
      par(mfrow=c(1,1))
      plot(thisLandCover)
      plot(grdPolys, border="blue",add = T)
      plot(sp, col="black", add=T, pch="+",cex=2)
      plot(vicBdry, add=T, border="grey40",lwd=2)
      color <- c(2,2,3,4,5) 
      color_transparent <- adjustcolor(color, alpha.f = 0.5) 
      plot(thisCell,col=color_transparent,add=T)
      
      thisVeg <- mask(crop(thisLandCover,thisCell),vicBdry)
      # plot(thisElev)
      # plot(grdPolys, border="blue",add = T)
      # plot(initialRegion, add=T, border="grey40")
      allVals <- getValues(thisVeg)
      validVals <- allVals[!is.na(allVals)]
      fractionInCell <- length(validVals) / length(allVals)
      counts <- plyr::count(validVals)
      fraction <- signif(counts$freq / sum(counts$freq),4)
      residual <- 1 - sum(fraction)
      if(residual != 0) fraction[which.max(fraction)] <- fraction[which.max(fraction)] + residual # true up for fractions to equal one
      outTable <- matrix(c(seq(1,20),rep(NA,20)),nc=2,byrow=F)
      outTable[,2] <- NA
      thisTable <- cbind(counts$x,fraction)
      outTable[match(thisTable[,1], outTable[,1]),2]
      outTable[match(thisTable[,1],outTable[,1]),2] <- thisTable[,2]
      return(outTable[,2])
    }
    vegFracMatrix <- do.call(cbind, lapply(1:nmbrCells, processVeg))
    colSums(vegFracMatrix, na.rm = T) # should all be "1"
    colnames(vegFracMatrix) <- paste(coordinates(grdPolys)[,1],coordinates(grdPolys)[,2],sep = "_")
    rownames(vegFracMatrix) <- landCvrLegendData$Class.Name
  }
   
  # Write to the netcdf
  
  dimVegClasses <- ncdim_def("veg_class","",c(1:20),create_dimvar=FALSE )
  dimRootZones <- ncdim_def("root_zone","",1:2,create_dimvar=FALSE )
  dimMonth <- ncdim_def("month","",1:12,create_dimvar=FALSE )
  dimnchar <- ncdim_def("nchar", "", 1:40, create_dimvar=FALSE )
  
  var_month <- ncvar_def("month","",dimMonth,prec = "integer", longname = "month of year")
  var_veg_class <- ncvar_def("veg_class","N/A",dimVegClasses,prec = "integer", longname = "veg_class")
  var_veg_descr <- ncvar_def("veg_descr","N/A",list(dimnchar,dimVegClasses),prec = "char", longname = "Vegetation Class Description")
  var_root_zone <- ncvar_def("root_zone","",dimRootZones,prec = "integer", longname = "root zone")
  var_Nveg <- ncvar_def("Nveg","N/A",list(dimLon,dimLat),-999,
                        prec = "integer", longname = "Nveg")
  var_Cv <- ncvar_def("Cv","fraction",list(dimLon,dimLat,dimVegClasses),missingValue,
                      prec = "double", longname = "Fraction of grid cell covered by vegetation tile")
  var_root_depth <- ncvar_def("root_depth","m",list(dimLon,dimLat,dimRootZones,dimVegClasses),missingValue,
                              prec = "double", longname = "Root zone thickness (sum of depths is total depth of root penetration)")
  var_root_fract <- ncvar_def("root_fract","fraction",list(dimLon,dimLat,dimRootZones,dimVegClasses),missingValue,
                              prec = "double", longname = "Fraction of root in the current root zone")
  var_LAI <- ncvar_def("LAI","m2/m2",list(dimLon,dimLat,dimMonth,dimVegClasses),missingValue,
                       prec = "double", longname = "Leaf Area Index, one per month")
  var_overstory <- ncvar_def("overstory","N/A",list(dimLon,dimLat,dimVegClasses),missval = -999,
                             prec = "integer", longname = "Flag to indicate whether or not the current vegetation type has an overstory (1 for overstory present [e.g. trees], 0 for overstory not present [e.g. grass])")
  var_rarc <- ncvar_def("rarc","s/m",list(dimLon,dimLat,dimVegClasses),missingValue,
                        prec = "double", longname = "Architectural resistance of vegetation type (~2 s/m)")
  var_rmin <- ncvar_def("rmin","s/m",list(dimLon,dimLat,dimVegClasses),missingValue,
                        prec = "double", longname = "Minimum stomatal resistance of vegetation type (~100 s/m)")
  var_wind_h <- ncvar_def("wind_h","m",list(dimLon,dimLat,dimVegClasses),missingValue,
                          prec = "double", longname = "Height at which wind speed is measured")
  var_RGL <- ncvar_def("RGL","W/m2",list(dimLon,dimLat,dimVegClasses),missingValue,
                       prec = "double", longname = "Minimum incoming shortwave radiation at which there will be transpiration. For trees this is about 30 W/m^2, for crops about 100 W/m^2.")
  var_rad_atten <- ncvar_def("rad_atten","fraction",list(dimLon,dimLat,dimVegClasses),missingValue,
                             prec = "double", longname = "Radiation attenuation factor. Normally set to 0.5, though may need to be adjusted for high latitudes")
  var_wind_atten <- ncvar_def("wind_atten","fraction",list(dimLon,dimLat,dimVegClasses),missingValue,
                              prec = "double", longname = "Wind speed attenuation through the overstory. The default value has been 0.5")
  var_trunk_ratio <- ncvar_def("trunk_ratio","fraction",list(dimLon,dimLat,dimVegClasses),missingValue,
                               prec = "double", longname = "Ratio of total tree height that is trunk (no branches). The default value has been 0.2.")
  var_albedo <- ncvar_def("albedo","fraction",list(dimLon,dimLat,dimMonth,dimVegClasses),missingValue,
                          prec  = "double",longname = "Shortwave albedo for vegetation type")
  var_veg_rough <- ncvar_def("veg_rough","m",list(dimLon,dimLat,dimMonth,dimVegClasses),missingValue,
                             prec = "double", longname = "Vegetation roughness length (typically 0.123 * vegetation height)")
  var_displacement <- ncvar_def("displacement","m",list(dimLon,dimLat,dimMonth,dimVegClasses),missingValue,
                                prec = "double", longname = "Vegetation displacement height (typically 0.67 * vegetation height)")
  
  # the below are created using the createCombinedVegTable.R
  combinedVegTable <- read.csv(file=paste0(climaticsCodeDir,'Master/combinedVegTable.csv'))
  glcRemapTable <- read.csv(file=paste0(climaticsCodeDir,'Master/glcRemapTable.csv'))
  
  # Write Vegetation netcdf portion
  writeVeg <- T # I just put this if statement in to help with debugging
  if (writeVeg){
    ncid_old <- nc_open(paramsNcdfFile, write=T)
    ncid_old <- ncvar_add(ncid_old,var_veg_class)
    ncvar_put(ncid_old,"veg_class",c(1:20))
    
    # veg_descr
    ncid_old <- ncvar_add(ncid_old,var_veg_descr)
    ncvar_put(ncid_old,"veg_descr",glcRemapTable$glcVegDesc)
    
    # root_zone
    ncid_old <- ncvar_add(ncid_old,var_root_zone)
    ncvar_put(ncid_old,"root_zone",c(1:2))
    
    #month
    ncid_old <- ncvar_add(ncid_old,var_month)
    ncvar_put(ncid_old,"month",c(1:12))
    
    #Nveg
    ncid_old <- ncvar_add(ncid_old,var_Nveg)
    xyz <- data.frame(cbind(coordinates(vicGrid),z=rep(NA,dim(coordinates(vicGrid))[2])))
    xyMap <- paste(coordinates(vicGrid)[,1],coordinates(vicGrid)[,2],sep = "_") 
    indices <- match(xyMap,colnames(vegFracMatrix))
    indices <- indices[!is.na(indices)]
    xyz$z[indices] <- colSums(is.na(vegFracMatrix))
    nVegRaster <-rasterFromXYZ(xyz)
    outMatrix <- as.matrix(myFlip(nVegRaster))
    ncvar_put(ncid_old,"Nveg",outMatrix)
    
    # Cv: Fraction of grid cell covered by vegetation tile
    ncid_old <- ncvar_add(ncid_old,var_Cv)
    vicCoords <- coordinates(vicGrid)
    xyz <- data.frame(cbind(vicCoords,matrix(NA,nc=20,nr=dim(vicCoords)[1])))
    xyMap <- paste(coordinates(vicGrid)[,1],coordinates(vicGrid)[,2],sep = "_") 
    indices <- match(xyMap,colnames(vegFracMatrix))
    indices <- indices[!is.na(indices)]
    xyz[indices,3:dim(xyz)[2]] <- vegFracMatrix
    cvBrick <-rasterFromXYZ(xyz)
    outArray <- as.array(myFlip(cvBrick))
    outArray[outArray == missingValue] <- 0
    outArray <- outArray * vicGrid
    ncvar_put(ncid_old,"Cv",outArray)
    
    # root_depth
    ncid_old <- ncvar_add(ncid_old,var_root_depth)
    maskVector <- getValues(vicGrid)
    vicCoords <- coordinates(vicGrid)
    vecLength <- dim(vicCoords)[1]
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      rootDepth1 <- combinedVegTable$root_depth_rz1[glcRemapTable$remap][idx]
      rDVector1 <- rep(rootDepth1,vecLength) * maskVector
      rootDepth2 <- combinedVegTable$root_depth_rz2[glcRemapTable$remap][idx]
      rDVector2 <- rep(rootDepth2,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,rDVector1,rDVector2))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"root_depth",outArray)
   
    # root_fract
    ncid_old <- ncvar_add(ncid_old,var_root_fract)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      rootDepth1 <- combinedVegTable$root_fract_rz1[glcRemapTable$remap][idx]
      rDVector1 <- rep(rootDepth1,vecLength) * maskVector
      rootDepth2 <- combinedVegTable$root_fract_rz2[glcRemapTable$remap][idx]
      rDVector2 <- rep(rootDepth2,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,rDVector1,rDVector2))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"root_fract",outArray)
  
  
    # LAI (veg_class, month, lat, lon)
    ncid_old <- ncvar_add(ncid_old,var_LAI)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      laiMax <- combinedVegTable$LAImax[glcRemapTable$remap][idx]
      laiMaxVector <- rep(laiMax,vecLength) * maskVector
      xyzMax <- data.frame(cbind(vicCoords,laiMaxVector))
      xyzMaxRaster <- stack(replicate(12, rasterFromXYZ(xyzMax)))
      
      laiMin <- combinedVegTable$LAImin[glcRemapTable$remap][idx]
      laiMinVector <- rep(laiMin,vecLength) * maskVector
      xyzMin <- data.frame(cbind(vicCoords,laiMinVector))
      xyzMinRaster <- stack(replicate(12, rasterFromXYZ(xyzMin)))
      # vector for the monthly values to use the summer (max) weighting, '1 - this' equals the winter (min) weighting
      #                    J, F, M,  A,   M, J, J, A, S, O,  N, D
      summerWeighting <- c(0, 0, 0, .25, .5, 1, 1, 1, 1, .5, 0, 0)
      foo <- xyzMaxRaster * summerWeighting + xyzMinRaster * (1 - summerWeighting)
      return(foo)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"LAI",outArray)
  
    # overstory (veg_class, lat, lon) 
    ncid_old <- ncvar_add(ncid_old,var_overstory)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      r1 <- combinedVegTable$overstory[glcRemapTable$remap][idx]
      v1 <- rep(r1,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,v1))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    mode(outArray) <- "integer"
    ncvar_put(ncid_old,"overstory",outArray)
  
    # rarc (veg_class, lat, lon) 
    ncid_old <- ncvar_add(ncid_old,var_rarc)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      r1 <- combinedVegTable$rarc[glcRemapTable$remap][idx]
      v1 <- rep(r1,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,v1))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"rarc",outArray)
    
    # rmin (veg_class, lat, lon) 
    ncid_old <- ncvar_add(ncid_old,var_rmin)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      r1 <- combinedVegTable$rmin[glcRemapTable$remap][idx]
      v1 <- rep(r1,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,v1))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"rmin",outArray)
    
    # wind_h (veg_class, lat, lon) 
    ncid_old <- ncvar_add(ncid_old,var_wind_h)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      r1 <- combinedVegTable$wind_h[glcRemapTable$remap][idx]
      v1 <- rep(r1,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,v1))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"wind_h",outArray)
    
    # RGL (veg_class, lat, lon) 
    ncid_old <- ncvar_add(ncid_old,var_RGL)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      r1 <- combinedVegTable$RGL[glcRemapTable$remap][idx]
      v1 <- rep(r1,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,v1))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"RGL",outArray)
    
    # rad_atten (veg_class, lat, lon) 
    ncid_old <- ncvar_add(ncid_old,var_rad_atten)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      r1 <- combinedVegTable$rad_atten[glcRemapTable$remap][idx]
      v1 <- rep(r1,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,v1))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"rad_atten",outArray)
    
    # wind_atten (veg_class, lat, lon) 
    ncid_old <- ncvar_add(ncid_old,var_wind_atten)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      r1 <- combinedVegTable$wind_atten[glcRemapTable$remap][idx]
      v1 <- rep(r1,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,v1))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"wind_atten",outArray)
    
    # trunk_ratio (veg_class, lat, lon) 
    ncid_old <- ncvar_add(ncid_old,var_trunk_ratio)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      r1 <- combinedVegTable$trunk_ratio[glcRemapTable$remap][idx]
      v1 <- rep(r1,vecLength) * maskVector
      xyz <- data.frame(cbind(vicCoords,v1))
      outBrick <-rasterFromXYZ(xyz)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"trunk_ratio",outArray)
    
    # albedo (veg_class, month, lat, lon)
    ncid_old <- ncvar_add(ncid_old,var_albedo)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      roughRaster <- combinedVegTable$roughnessLengthM[glcRemapTable$remap][idx]
      roughVector <- rep(roughRaster,vecLength) * maskVector
      xyzDF <- data.frame(cbind(vicCoords,roughVector))
      xyzStack <- stack(replicate(12, rasterFromXYZ(xyzDF)))
      return(xyzStack)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"albedo",outArray)
    
    # veg_rough (veg_class, month, lat, lon)
    ncid_old <- ncvar_add(ncid_old,var_veg_rough)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      roughRaster <- combinedVegTable$roughnessLengthM[glcRemapTable$remap][idx]
      roughVector <- rep(roughRaster,vecLength) * maskVector
      xyzDF <- data.frame(cbind(vicCoords,roughVector))
      xyzStack <- stack(replicate(12, rasterFromXYZ(xyzDF)))
      return(xyzStack)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"veg_rough",outArray)
    
    # displacement (veg_class, month, lat, lon)
    ncid_old <- ncvar_add(ncid_old,var_displacement)
    outBrick <- do.call(addLayer, lapply(1:20, function(idx){
      displaceRaster <- combinedVegTable$displacement[glcRemapTable$remap][idx]
      displaceVector <- rep(displaceRaster,vecLength) * maskVector
      xyzDF <- data.frame(cbind(vicCoords,displaceVector))
      xyzStack <- stack(replicate(12, rasterFromXYZ(xyzDF)))
      return(xyzStack)
    }))
    outArray <- as.array(myFlip(outBrick))
    ncvar_put(ncid_old,"displacement",outArray)
    
    nc_close(ncid_old)
  }
  
}
