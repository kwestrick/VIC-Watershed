driverVICInputsClassic <- function(watershedPourPointName = "KeenleysideDam",
                                   resolutionKm = 16.){
  #
  # This is the main driver to set up a watershed and all required input files 
  # for VIC version 5. It is built around many of the routines that were R&D'd
  # in the Code/ directory.
  #
  # K. Westrick 16 March 2019
  #
  options(stringsAsFactors = F)
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,rgdal,raster,openxlsx,ggmap,rgeos)
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source("/Volumes/MiniPro/Climatics/Code/Master/utilities.R")
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICDomain.R"))
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICForcingClassic.R"))
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICSoilClassic.R"))
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICSnowbandsClassic.R"))
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICVegClassicNLDAS11.R"))
  # ,ggplot2,rgeos,tidyverse,rvest,devtools,git2r,ggmap,
  #               ggrepel,git2r, openxlsx)
  
  thisDTG <- substr(as.character(gsub("-","",Sys.time())),1,8)
  thisDTG <- "20190316"
  thisFilename <- paste0(tolower(watershedPourPointName),".",thisDTG)
  saveDir <- paste0(climaticsOpDir,"Data/",
                    watershedPourPointName,"_",resolutionKm,"km")
  if(!dir.exists(saveDir)) dir.create(saveDir)

  ################
  #
  # Determine and configure the watershed based on pour point
  #
  ################
  
  ##
  determinePourPoint <- function(searchTerm = watershedPourPointName){
    energyGpsLocs <- data.frame(read.xlsx(paste0(climaticsDataDir,"/Client/EnergyGPS/EGPS_HydroDams_BasinBreakout.xlsx")))
    hydroSites17 <- readRDS(paste0(climaticsDataDir,"RDS/columbiaRiverDamSitesSpatialPoints.rds"))
    damLocs <- readRDS(paste0(climaticsDataDir,"RDS/damLocs.rds"))
    damLocs <- damLocs[!damLocs$STATE %in% c("PR","HI","AK"),]
    damLocs <- spTransform(damLocs,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    damLocs <- damLocs[,c("DAMS00X020","DAM_NAME","LATITUDE","LONGITUDE","RIVER","PURPOSES","STATE","storM3")]
    names(damLocs) <- c("ID","NAME","LATITUDE","LONGITUDE","RIVER","PURPOSES","STATE","STORM3")
    damLocs$TYPE <- "Dam"
    purposeTable <- data.frame(id=c("C","D","F","H","I","N","O","P","R","S","T"),
                               text=c("flood control & storm water management",
                                      "debris control",
                                      "fish & wildlife pond",
                                      "hydroelectric",
                                      "irrigation",
                                      "navigation",
                                      "other",
                                      "fire protection, stock, or small farm pond",
                                      "recreation",
                                      "water supply",
                                      "tailings"))
    purposesString <- strsplit(damLocs$PURPOSES,"")
    damLocs$PURPOSES_STRING <- do.call(rbind, lapply(purposesString, function(x){
      paste(purposeTable$text[match(x,purposeTable$id)],collapse = "; ")
    }))
    gages2Locs <- readRDS(paste0(climaticsDataDir,"RDS/gages2.rds"))
    gages2Locs <- gages2Locs[!gages2Locs$STATE %in% c("PR","HI","AK"),]
    gages2Locs <- spTransform(gages2Locs,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    gages2Locs <- gages2Locs[,c("STAID","STANAME","LAT_GAGE","LNG_GAGE","STATE","ACTIVE09","HCDN_2009","DRAIN_SQKM")]
    names(gages2Locs) <- c("ID","NAME","LATITUDE","LONGITUDE","STATE","ACTIVE09","HCDN09","DRAIN_SQKM")
    gages2Locs$TYPE <- "Gage"
    plot(gages2Locs)
    thisSite <- gages2Locs[grepl(searchTerm,gages2Locs$NAME),]
    #thisSite <- damLocs[grepl(searchTerm,damLocs$NAME),]
    #thisSite <- data.frame(lat=thisSite$LATITUDE,lon=thisSite$LONGITUDE,name=thisSite$NAME)
    #thisSite <- hydroSites17[grepl(searchTerm,hydroSites17$DAM_NAME),]
    #thisSite <- energyGpsLocs[grepl(searchTerm,energyGpsLocs$DamName)]
    return(thisSite)
  }
  thisSiteSpPt <- determinePourPoint()
  thisSite <- data.frame(name=watershedPourPointName, lat=thisSiteSpPt$LATITUDE,lon=thisSiteSpPt$LONGITUDE)
  
  ## Canada dams in the upper Columbia
  canadaDams <- data.frame(name=c("Mica Dam","Revelstoke Dam","Keenleyside Dam"),
                           lat=c(52.077,51.049,49.338),lon=c(-118.5663,-118.1933,-117.7716))
  thisSite <- canadaDams[3,]
  
  columbiaBdry <- readRDS(paste0(climaticsDataDir,"RDS/columbiaBdryPolygonExpanded.rds"))
  
  ## this routine is abbreviated to use filled dems and other files created from the "full" routine (below)
  foo <- defineVICWatershedPourPointMethodAbbreviated(inBoundaryPoly=columbiaBdry,
    lat=thisSite$lat,lon=thisSite$lon, 
    title = thisSite$name)
  myWatershed <- spTransform(foo,projVals)
  raster::area(myWatershed) / 1e6
  
  ## the "full" routine
  # thisWatershedBdry <- defineVICWatershedPourPointMethodAbbreviated(inBoundaryPoly = columbiaBdry, 
  #                                   pourPointCoords = pourPoint[,1:2], 
  #                                   title = "STEHEKIN")
  thisWatershedBdry <- spTransform(readRDS(paste0(climaticsDataDir,"VICWatersheds/",
                                                  watershedPourPointName,"/shape_",
                                                  watershedPourPointName,".rds")),
                                   crs(projVals))
  watershedSizeErrorPercent <- round(abs(((raster::area(thisWatershedBdry) / 1e6) / thisSite$DRAIN_SQKM) - 1) * 100,2)
  cat(paste0("watershed size error is ",watershedSizeErrorPercent,"%"),fill=T)
  watershedPourPoint <- data.frame(lat=thisSite$lat,lon=thisSite$lon)
  #thisWatershedFcstPts <- readRDS(paste0(climaticsDataDir,"RDS/columbiaRiverFcstPoints.rds"))
  #foo <- over(thisWatershedBdry,thisWatershedFcstPts)
  plotPt <- fortify(as.data.frame(thisSite))
  c1 <- coordinates(gCentroid(thisWatershedBdry))
  googleApiKey <- "AIzaSyA_nSoEvdT_zfO4q_FUucG_XfOpFtqHduk"
  ggmap::register_google(key = googleApiKey)
  maptype <- "terrain" #"terrain-background"
  myMap <- ggmap(get_map(c1, zoom = 7, maptype = maptype)) 
  myMap +
    geom_polygon(aes(x = long, y = lat), data = thisWatershedBdry,
                 colour = 'red', fill = 'transparent', alpha = .15, size = .4, lwd=3) +
    geom_point(aes(x = thisSite$lon, y = thisSite$lat), shape=23, 
               fill = "orange", color= "black", size = 4, data = plotPt) +
    theme(legend.position="bottom") + ggtitle("Forecast and Validation Points") +
    annotate("text",x = c1[,1], y = 37, label = 
               "Copyright 2018 Climatics Co - Not for Reproduction", 
             colour = "blue",size=5.5) +
    annotate("text",x = c1[,1], y = 37.5, label = 
               "atop(italic('source code: createColumbiaFcstPoints.R'))", 
             colour = "black",size=4,parse = T)
  #
  ################
  #
  # Create the domain file for VIC, this includes making and resizing the raster mask
  # and returning this raster, based on the 
  #
  ################
  vicDomList <- createVICDomain(vicBdry = thisWatershedBdry, resolutionKm = resolutionKm)
  vicGrid <- vicDomList[["MASK"]]
  vicCellNmbr <- vicDomList[["GRID"]]
  LAT <- matrix(coordinates(vicGrid)[,2],nc=1)
  colnames(LAT) <- "LAT"
  LNG <- matrix(coordinates(vicGrid)[,1],nc=1)
  colnames(LNG) <- "LNG"
  RUN <- matrix(vicDomList[["RUN"]],nr=ncell(vicDomList[["RUN"]]))
  colnames(RUN) <- "RUN"
  runCell <- vicDomList[["RUN"]]
  runCell
  GRID <- matrix(vicDomList[["GRID"]],nr=ncell(vicDomList[["GRID"]]))
  colnames(GRID) <- "GRID"
  gridCell <- vicDomList[["GRID"]]
  AREA <- matrix(vicDomList[["AREA"]],nr=ncell(vicDomList[["AREA"]])) 
  colnames(AREA) <- "AREA"
  FRAC <- matrix(vicDomList[["FRAC"]],nr=ncell(vicDomList[["FRAC"]])) 
  colnames(FRAC) <- "FRAC"
  plot(vicGrid,colNA='brown')
  plot(thisWatershedBdry,border='black',lwd=2,add=T)
  stackData <- stack(unlist(vicDomList))
  names(stackData) <- names(vicDomList)
  saveStackFile <- paste0(saveDir,"/stack_allRasters_",watershedPourPointName,"_",resolutionKm,"km")
  writeRaster(stackData,file=saveStackFile,overwrite=T)
  saveWatershedFile <- paste0(saveDir,"/watershed_shape_",watershedPourPointName,"_",resolutionKm,"km.rds")
  saveRDS(thisWatershedBdry,file=saveWatershedFile)
  
  
  ################
  #
  # Create the soil file for VIC
  #
  ################
  
  soil <- createVICSoilClassic(vicBdry = thisWatershedBdry,vicGrid = runCell, gridCell=gridCell,
                      resolutionKm = resolutionKm)
  soilParamFile <- paste0(saveDir,"/soilParam_",thisFilename,
                          "_",resolutionKm,"km.txt")
  if(file.exists(soilParamFile)) unlink(soilParamFile)
  write.table(soil,soilParamFile,append = F, sep="\t",col.names = F,row.names = F)
  
  ################
  #
  # Create the snow file for VIC
  #
  ################
  
  bands <- createVICSnowbandsClassic(vicBdry = thisWatershedBdry, vicGrid = vicGrid, 
                                  vicCellNmbr = vicCellNmbr, runCell = runCell,
                                  resolutionKm = resolutionKm)
  snowBandsFile <- paste0(saveDir,"/snowBands_",thisFilename,
                         "_",resolutionKm,"km.txt")
  if(file.exists(snowBandsFile)) unlink(snowBandsFile)
  write.table(bands,snowBandsFile,append = F, sep="\t",col.names = F,row.names = F)
  
  ################
  #
  # Create the Veg file for VIC
  #
  ################
  
  vegList <- createVICVegClassicNLDAS11(vicBdry = thisWatershedBdry, vicGrid = vicGrid, 
                                vicCellNmbr = vicCellNmbr, runCell = runCell,
                                resolutionKm = resolutionKm)
  veg <- vegList[[1]]
  vegParamFile <- paste0(saveDir,"/vegParam_",thisFilename,
                            "_",resolutionKm,"km.txt")
  if(file.exists(vegParamFile)) unlink(vegParamFile)
  for(idx in 1:length(veg)){
    vals <- veg[[idx]]
    nmbrVegTiles <- dim(vals)[1]
    headerLine <- matrix(c(idx,nmbrVegTiles),nc=2,nr=1)
    write.table(headerLine,vegParamFile,append = T, sep="\t",col.names = F,row.names = F)
    valsDF <- data.frame(vals)
    veg1Line <- valsDF[,1:6]
    veg2Line <- valsDF[,7:18]
    for(vegTileIdx in 1:nmbrVegTiles){
      write.table(veg1Line[vegTileIdx,],vegParamFile,append = T, sep="\t",col.names = F,row.names = F)
      write.table(veg2Line[vegTileIdx,],vegParamFile,append = T, sep="\t",col.names = F,row.names = F)
      
    }
  }
  
  vegLib <- vegList[[2]]
  vegLibFile <- paste0(saveDir,"/vegLib_",thisFilename,
                         "_",resolutionKm,"km.txt")
  
  if(file.exists(vegLibFile)) unlink(vegLibFile)
  write.table(vegLib,vegLibFile,append = F, sep="\t",col.names = F,row.names = F)
 
  ################
  #
  # Create the met forcing file input for VIC
  #
  ################
  
  yrs <- as.character(seq(2010,2019,1))
  for(yr in yrs){
    createVICForcingClassic(vicGrid = vicGrid, 
                            yr = yr,
                            thisFilename = thisFilename, 
                            saveDir = saveDir)
  }
}
