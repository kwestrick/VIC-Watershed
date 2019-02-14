driverVICSetup <- function(watershedPourPointName = "STEHEKIN"){
  #
  # This is the main driver to set up a watershed and all required input files 
  # for VIC version 5. It is built around many of the routines that were R&D'd
  # in the Code/ directory.
  #
  # K. Westrick 29 Jan 2019
  #
  options(stringsAsFactors = F)
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,rgdal,raster,openxlsx,ggmap,rgeos)
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source("/Volumes/MiniPro/Climatics/Code/Master/utilities.R")
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICDomain.R"))
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICForcing.R"))
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICParameter.R"))
  # ,ggplot2,rgeos,tidyverse,rvest,devtools,git2r,ggmap,
  #               ggrepel,git2r, openxlsx)
  
  thisDTG <- substr(as.character(gsub("-","",Sys.time())),1,8)
  thisDTG <- "20190203"
  thisFilename <- paste0(tolower(watershedPourPointName),".",thisDTG)
  
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
  thisSite <- determinePourPoint()
  
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
                                                  watershedPourPointName,"/shape_STEHEKIN.rds")),
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
  myMap <- ggmap(get_map(c1, zoom = 10, maptype = maptype)) 
  myMap +
    geom_polygon(aes(x = long, y = lat), data = thisWatershedBdry,
                 colour = 'red', fill = 'transparent', alpha = .15, size = .3) +
    geom_point(aes(x = LONGITUDE, y = LATITUDE), shape=25, 
               fill = "green", color= "blue", size = 4, data = plotPt) +
    theme(legend.position="bottom") + ggtitle("Forecast and Validation Points") +
    annotate("text",x = c1[,1], y = 37, label = 
               "Copyright 2018 Climatics Co - Not for Reproduction", 
             colour = "blue",size=5.5) +
    annotate("text",x = c1[,1], y = 37.5, label = 
               "atop(italic('source code: createColumbiaFcstPoints.R'))", 
             colour = "black",size=4,parse = T)
  #### FOO ###
  # just use the Stehekin test data right now
  #
  readStehekin <- function(){
    stehekinDir <- "/Volumes/MiniPro/Climatics/Code/VIC"
    stehekinDomain <- paste0(stehekinDir,"/vic/drivers/image/parameters/domain.stehekin.20151028.nc")
    domainPtr <- nc_open(stehekinDomain)
    as.matrix(names(domainPtr$var))
    stehekinLand <- paste0(stehekinDir,"/vic/drivers/image/parameters/Stehekin_test_params_20160327.nc")
    domainPtr <- nc_open(stehekinLand)
    as.matrix(names(domainPtr$var))
    elev <- ncvar_get(domainPtr,varid = "elev")
    lats <- ncvar_get(domainPtr,varid = "lat")
    lons <- ncvar_get(domainPtr,varid = "lon")
    lonLats <- expand.grid(lons,lats)
    mask <- ncvar_get(domainPtr,varid = "mask")
    maskRaster <- rasterFromXYZ(cbind(lonLats,as.vector(mask)),crs=projVals)
    plot(maskRaster)
  }
  #readStehekin()
  #
  ################
  #
  # Create the domain file for VIC, this includes making and resizing the raster mask
  # and returning this raster, based on the 
  #
  ################
  resolutionKm = 14.
  vicGrid <- createVICDomain(vicBdry = thisWatershedBdry, resolutionKm = resolutionKm,
                  thisFilename = thisFilename)
  
  ################
  #
  # Create the met forcing file for VIC
  #
  ################
  
  createVICForcing(maskRaster = vicGrid, yr = "1980",
                   thisFilename = thisFilename)
  
  ################
  #
  # Create the parameters file for VIC
  #
  ################
  createVICParameters(vicBdry = thisWatershedBdry,vicGrid = vicGrid, thisFilename = thisFilename,
                      resolutionKm = resolutionKm)

}
