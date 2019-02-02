createVICParameter <- function(){
  # Requires various raw and prepared spatial fields 
  #
  #   Soil: processSoilForVic.R
  #   Veg:  processVegForVic.R
  #   DEM:  processDemForVic.R (also does snow fields)
  
  options(stringsAsFactors = F)
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,rgdal,raster,ggplot2,rgeos,tidyverse,rvest,devtools,git2r,ggmap,
                 ggrepel,git2r, openxlsx)
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source("/Volumes/MiniPro/Climatics/Code/Master/utilities.R")
  source('/Volumes/MiniPro/Climatics/Code/VIC-Setup/defineWatershed.R')
  source(paste0(climaticsCodeDir,"Reference/retrieveDamSites.R"))
  source(paste0(climaticsCodeDir,"Master/createLandCoverImage.R"))
  if(!requireNamespace("devtools")) install.packages("devtools")
  ### devtools::install_github("dkahle/ggmap", ref = "tidyup")  
  
  # in order to make the raster look like what's needed in netcdf we need to rotate the matrix 90 degrees clockwise. 
  rotate <- function(x) t(apply(x, 2, rev))
  # this routine applies the rotation to raster stacks or bricks
  stack2array <- function(inArray=0){
    return(sapply(1:dim(inArray)[3], function(x)
      rotate(inArray[,,x]), simplify = "array"))
  }
  
  # include the below for inclusion in the soil netcdf

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
  if(length(missingSoilParams) == 0) cat("All base files are available") else {
    for(i in 1:length(missingSoilParams)){
      cat(paste0("missing ",missingSoilParams[i],fill=T))
    }
    stop("run processSoilForVic.R")
  }
  
  # references: http://geog.uoregon.edu/bartlein/courses/geog490/week04-netCDF.html#reading-restructuring-and-writing-netcdf-files-in-r
  #
  #
  # define region (based on our current flawed geographic identification a
  #
  huc17 <- spTransform(readOGR(dsn=paste0(climaticsDataDir,"Watersheds_Shapefile/HUC17ShapeWithCanada/WBDHU6.shp")), projVals)
  #huc17Simp <- rgeos::gSimplify(huc17,tol=.002,topologyPreserve = T)
  plot(huc17, col = rainbow(n=19))
  extent(huc17)
  centroids.df <- as.data.frame(coordinates(huc17))
  names(centroids.df) <- c("long","lat")
  id.df <- huc17$Name
  idNames <- data.frame(centroids.df,id.df)
  idHucs <- data.frame(centroids.df,huc=huc17$HUC6)
  map <- ggplot() + geom_polygon(data = huc17, aes(x = long, y = lat, group = group), colour = "black", fill = NA) +
  geom_text(data = idNames,aes(label = id.df, x = long, y = lat)) +
  geom_text(data = idHucs,aes(label = huc, x = long, y = lat-0.25)) 
  map

  # get the locations for the damms on the Columbia from Wikipedia
  
  inURL <- "https://en.wikipedia.org/wiki/List_of_dams_in_the_Columbia_River_watershed"
  temp <- inURL %>% 
    read_html %>%
    html_nodes("table")
  inTable1 <- data.frame(html_table(temp[2])) ## Just the "legend" table
  inTable2 <- data.frame(html_table(temp[3])) ## Just the "legend" table
  inTable3 <- data.frame(html_table(temp[15])) ## Just the "legend" table
  inTable4 <- data.frame(html_table(temp[16])) ## Just the "legend" table
  inTable5 <- data.frame(html_table(temp[24])) ## Just the "legend" table
  inTable6 <- data.frame(html_table(temp[28])) ## Just the "legend" table
  inTable7 <- data.frame(html_table(temp[29])) ## Just the "legend" table
  inTable8 <- data.frame(html_table(temp[33])) ## Just the "legend" table
  inTable9 <- data.frame(html_table(temp[34])) ## Just the "legend" table
  inTable10 <- data.frame(html_table(temp[35])) ## Just the "legend" table
  inTable11 <- data.frame(html_table(temp[43])) ## Just the "legend" table
  inTable12 <- data.frame(html_table(temp[44])) ## Just the "legend" table
  
  
  inTable <- rbind(inTable1,inTable2)
                   #,inTable3,inTable4,
                   #inTable5,inTable6,inTable7,inTable8,
                   #inTable9,inTable10,inTable11,inTable12)  
  inTable <- inTable[!inTable$Capacity..MW. == "0",]

  inCoords <- trimws(matrix(unlist(strsplit(inTable$Coordinates, "/")),nc=3,byrow = T)[,2])
  lons <- -(as.numeric(substr(matrix(unlist(strsplit(inCoords," ")),nc=2,byrow=T)[,2],1,8)))
  lats <- as.numeric(substr(matrix(unlist(strsplit(inCoords," ")),nc=2,byrow=T)[,1],2,7)) # not sure why substr has to start at 2 but it works
  capacity <- as.numeric(gsub(",","",matrix(unlist(strsplit(inTable$Capacity..MW.,"\\[")),nc=2,byrow=T)[,1]))
  name <- inTable$Name
  dams <- data.frame(lons,lats,capacity,name)
  sum(dams$capacity)

  # dam locations from the USGS

  damLocs <- retrieveDamSites()
  damLocsDf <- data.frame(damLocs@data)
  ht(damLocsDf)
  
  # hydropower locations from Oak Ridge National Laboratory
  # https://hydrosource.ornl.gov/market-info-and-data/existing-hydropower-assets/existing-hydropower-assets-datasets/national
  
  hydroPlantsFile <- paste0(climaticsDataDir,"/HydroPlants/GIS/ORNL_EHAHydroPlantV2_FY18Q3_Operational.shp")
  hydroPlantsShp <- readOGR(hydroPlantsFile)
  hydroPlantsDF <- data.frame(hydroPlantsShp@data)
  ht(hydroPlantsDF)
  
  mergedSites <- merge(hydroPlantsDF,damLocsDf,by.x = "NID_ID1",by.y = "NIDID")
  
  # Energy GPS forecast locations
  
  energyGpsLocs <- read.xlsx(paste0(climaticsDataDir,"/Client/EnergyGPS/EGPS_HydroDams_BasinBreakout.xlsx"))
  mergedSites <- merge(mergedSites,energyGpsLocs,by.x = "DAM_NAME",by.y = "DamName",all.y = T)
  mergedSites <- mergedSites[!duplicated(mergedSites$DAM_NAME),]
  ### foo: we need to go into the energy GPS xlsx file and put in a field to better merge the dam locations. Maybe use the NIDID?
  ### there are 31 locations where the DAM_NAME did not match a DamName.
  pairedPts_sp <- SpatialPointsDataFrame(coords = cbind(mergedSites$lon,mergedSites$lat), data = mergedSites, proj4string = CRS(projVals))
  
   # NWRFC forecast points https://www.nwrfc.noaa.gov/rfc/
  
  rfcFile <- paste0(climaticsDataDir,"/RiverForecastCenter/river_summary.html")
  rfcFcstPts <- read.csv(rfcFile)
  rfc
  ht(rfcFcstPts)
  # graphical map using ggmap

  c1 <- coordinates(gCentroid(huc17))
  hucTgts <- paste0("17",c("0200","0101","0102","0103","0603","0602","0402","0401",
                        "0601","0502","0501","0300","0701","0702","0703","0501"))
  initialRegion <- subset(huc17,huc17@data$HUC6 %in% hucTgts)
  initialRegionPlot <- fortify(initialRegion)
  
  # trim our forecast points down to 
  huc17Mask <- substr(mergedSites$HUC_8,1,6) %in% hucTgts
  hydroSites17 <- mergedSites[huc17Mask,]
  hydroSites17Plot <- fortify(hydroSites17)
  hydroSites17sp <- SpatialPointsDataFrame(coords = cbind(hydroSites17$lon,hydroSites17$lat), data = hydroSites17, proj4string = CRS(projVals))
  
  

  myMap <- ggmap(get_map(c1, zoom = 5, maptype ="satellite")) 
  myMap +
    geom_polygon(aes(x = long, y = lat, group = group), data = huc17,
                 colour = 'red', fill = 'transparent', alpha = .15, size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = initialRegion,
                 colour = 'blue', fill = 'red', alpha = .5, size = .3) +
    geom_point(aes(x = lon, y = lat,  fill = "lightblue", size = 3, alpha=.5), data = hydroSites17Plot) +
    theme(legend.position="none") + ggtitle("Major Dams on Main Stem Columbia and Snake") +
    ggrepel::geom_label_repel(data = dams, aes(x = lons, y = lats, label = name), 
                     fill = "white", box.padding = unit(.4, "lines"),
                     label.padding = unit(.15, "lines"),
                     segment.color = "red", segment.size = 1) +
    annotate("text",x = -118, y = 37, label = "Copyright 2018 Climatics Co - Not for Reproduction", colour = "orange",size=5.5) +
    annotate("text",x = -118, y = 36, label = "atop(italic('source code: setupVicColumbia.R'))", colour = "white",size=4,parse = T)
  
  
  
  # The parent VIC grid
  
  maskInfile <- paste0(climaticsDataDir,"/RDS/VICNativeGridUSA.raster.rds")
  nativeMask <- readRDS(file = maskInfile)
  boundaryInfile <- paste0(climaticsDataDir,"/RDS/VICNativeExtentUSA.shape.rds")
  nativeBound <- readRDS(file = boundaryInfile)
  plot(nativeMask, col="deepskyblue", legend=F)
  plot(nativeBound, add=T, border="grey40")
  
  # create the Columbia grid from the VIC native grid
  
  initialExtent <- extent(initialRegion)
  t1 <- cellsFromExtent(nativeMask, initialExtent, expand = T)
  initialExtentSnap <- extentFromCells(nativeMask,t1)
  columbiaMask <- crop(mask(nativeMask, initialRegion),initialExtentSnap)
  plot(columbiaMask)
  plot(initialRegion,add=T)
  plot(damLocs, add=T)
  
  # Change the resolution to something coarser, must go in increments of km
  
  fillValue <- -214748364
  resolutionKm <- 450. # km
  if (resolutionKm > 1 && resolutionKm %% 1 == 0) vicGrid <- aggregate(columbiaMask,fact = resolutionKm)
  vicGrid[is.nan(vicGrid)] <- NA
  plot(vicGrid)
  #
  # create run_cell mask
  runCell <- vicGrid
  runCell[!is.na(runCell)] <- 1
  runCell[is.na(runCell)] <- 0
  plot(runCell)
  
  # gridcell; sequence of numbers for the grid cells
  gridcell <- runCell
  gridcell[gridcell == 1] <- seq(1,ncell(gridcell[gridcell == 1]))
  plot(gridcell)
  plot(initialRegion, add=T, border="grey40")
  
  # area
  gridArea <- area(gridcell) * 1e6
  totalVicArea <- area(columbiaMask) * 1e6
  vicArea <- aggregate(totalVicArea * columbiaMask, fun= sum, fact = resolutionKm)
  vicArea[is.nan(vicArea)] <- NA
  totalVicArea <- aggregate(totalVicArea, fun= sum, fact = resolutionKm)
  totalvicArea[is.nan(totalvicArea)] <- NA
  as.matrix(vicArea) / 1e6
  as.matrix(totalVicArea) / 1e6
  plot(vicArea,colNA="black")
  plot(initialRegion, add=T, border="grey40")
  vicArea <- round(vicArea,0)
  
  # fraction of grid cell active
  vicFrac <- signif(vicArea / totalVicArea,5)
  as.matrix(vicFrac)
  plot(vicFrac,colNA="black")
  plot(initialRegion, add=T, border="grey40")
  
  
  #################
  #
  # CREATE THE NETCDF FILES
  #
  #################
  
  #
  fillValue <- -214748364
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

  #
  # thisDTG <- substr(as.character(gsub("-","",Sys.time())),1,8)
  thisDTG <- "test"
  
  ##########
  #
  # create the domain file
  # 
  ##########
  
  var_lat <- ncvar_def("lat","degrees_north",dimLat,fillValue,"latitude of grid cell center",prec="single")
  var_lon <- ncvar_def("lon","degrees_east",dimLon,fillValue,"longitude of grid cell center",prec="single")
  var_mask <- ncvar_def("mask","",list(dimLon,dimLat),longname = "domain mask", prec = "integer")
  var_frac <- ncvar_def("frac","unitless",list(dimLon,dimLat),longname = "fraction of grid cell that is active", prec = "double")
  var_area <- ncvar_def("area","m2",list(dimLon,dimLat),longname = "area of grid cell", prec = "double")
  
  domainsNcdfFile <- paste0("columbia_domain_",thisDTG,".nc")
  domainsNcdf <- nc_create(domainsNcdfFile,list(var_mask,var_frac,var_area),verbose = F)
  ncatt_put(domainsNcdf,0,"_FillValue",fillValue)
  ncvar_put(domainsNcdf,var_mask,rotate(as.matrix(runCell)))
  vicFrac[is.na(vicFrac)] <- fillValue
  ncvar_put(domainsNcdf,var_frac,rotate(as.matrix(vicFrac)))
  gridArea[is.na(gridArea)] <- fillValue
  ncvar_put(domainsNcdf,var_area,rotate(as.matrix(gridArea)))
  nc_close(domainsNcdf)
  
  ############
  #
  # name and open the parameters file
  #
  ############
  
  # SOIL
  var_layer <- ncvar_def("layer","",list(dimLayers),longname = "soil_layer",prec = "integer")
  var_run_cell <- ncvar_def("run_cell","unitless",list(dimLon,dimLat), fillValue,"run_cell",prec = "integer")
  var_gridcell <- ncvar_def("gridcell","N/A",list(dimLon,dimLat),prec = "integer", longname = "grid_cell")
  var_lats <- ncvar_def("lats","degrees",list(dimLon,dimLat),prec = "double", longname = "Latitude of grid cell")
  var_lons <- ncvar_def("lons","degrees",list(dimLon,dimLat),prec = "double", longname = "Longitude of grid cell")
  var_infilt <- ncvar_def("infilt", "mm/day",list(dimLon,dimLat),fillValue,prec = "double", longname = "infilt")
  var_Ds <- ncvar_def("Ds","fraction",list(dimLon,dimLat),prec = "double", longname = "Ds")
  var_Dsmax <- ncvar_def("Dsmax","mm/day",list(dimLon,dimLat),prec = "double", longname = "Dsmax")
  var_Ws <- ncvar_def("Ws","fraction",list(dimLon,dimLat),prec = "double", longname = "Ws ")
  var_c <- ncvar_def("c","N/A",list(dimLon,dimLat),prec = "double", longname = "c")
  var_expt <- ncvar_def("expt","N/A",list(dimLayers, dimLon,dimLat),prec = "double", longname = "expt")
  var_Ksat <- ncvar_def("Ksat","mm/day",list(dimLayers,dimLon,dimLat),prec = "double", longname = "kSat")
  var_phi_s <- ncvar_def("phi_s","mm/mm",list(dimLayers,dimLon,dimLat),prec = "double", longname = "phi_s")
  var_init_moist <- ncvar_def("init_moist","mm",list(dimLayers,dimLon,dimLat),prec = "double", longname = "init_moist")
  var_elev <- ncvar_def("elev","",list(dimLon,dimLat),prec = "double", longname = "elev")
  var_depth <- ncvar_def("depth","m",list(dimLayers,dimLon,dimLat),prec = "double", longname = "depth")
  var_avg_T <- ncvar_def("avg_T","C",list(dimLon,dimLat),prec = "double", longname = "avg_T")
  var_dp <- ncvar_def("dp","m",list(dimLon,dimLat),prec = "double", longname = "dp")
  var_bubble <- ncvar_def("bubble","cm",list(dimLayers, dimLon,dimLat),prec = "double", longname = "bubble")
  var_quartz <- ncvar_def("quartz","fraction",list(dimLayers, dimLon,dimLat),prec = "double", longname = "quartz")
  var_bulk_density <- ncvar_def("bulk_density","kg/m3",list(dimLayers, dimLon,dimLat),prec = "double", longname = "bulk_density")
  var_soil_density <- ncvar_def("soil_density","kg/m3",list(dimLayers,dimLon,dimLat),prec = "double", longname = "soil_density")
  var_off_gmt <- ncvar_def("off_gmt","hours",list(dimLon,dimLat),prec = "double", longname = "off_gmt")
  var_Wcr_FRACT <- ncvar_def("Wcr_FRACT","fraction",list(dimLayers,dimLon,dimLat),prec = "double", longname = "Wcr_FRACT")
  var_Wpwp_FRACT <- ncvar_def("Wpwp_FRACT","fraction",list(dimLayers,dimLon,dimLat),prec = "double", longname = "Wpwp_FRACT")
  var_rough <- ncvar_def("rough","m",list(dimLon,dimLat),prec = "double", longname = "rough")
  var_snow_rough <- ncvar_def("snow_rough","m",list(dimLon,dimLat),prec = "double", longname = "snow_rough")
  var_annual_prec <- ncvar_def("annual_prec","mm",list(dimLon,dimLat),prec = "double", longname = "annual_prec")
  var_resid_moist <- ncvar_def("resid_moist","",list(dimLayers, dimLon,dimLat),prec = "double", longname = "resid_moist")
  var_fs_active <- ncvar_def("fs_active","binary ",list(dimLon,dimLat),prec = "integer", longname = "fs_active")
  var_July_Tavg <- ncvar_def("July_Tavg","C ",list(dimLon,dimLat),prec = "double", longname = "July_Tavg")
  
  
  paramsNcdfFile <- paste0("columbia_params_",thisDTG,".nc")
  varsForList <- c("var_run_cell","var_gridcell","var_lats","var_lons", paste0("var_",names(soilParams[which(soilParams == T)])))
  varList <- lapply(varsForList, get)

  paramsNcdf <- nc_create(paramsNcdfFile,varList,verbose = F)
  
  # process and place the basic fields in the netcdf
  # run_cell (0,1 mask)
  ncvar_put(paramsNcdf,var_run_cell,rotate(as.matrix(runCell)))
  ncatt_put(paramsNcdf,var_run_cell,"description","1 = Run Grid Cell, 0 = Do Not Run")
  
  # gridcell (sequental numbering)
  ncvar_put(paramsNcdf,var_gridcell,rotate(as.matrix(gridcell)))
  ncatt_put(paramsNcdf,var_gridcell,"description","Grid cell number")
  
  # lats and lons
  ncvar_put(paramsNcdf,var_lats,rotate(as.matrix(latRaster)))
  ncvar_put(paramsNcdf,var_lons,rotate(as.matrix(lonRaster)))
  
  # process and place all the soil fields in the netcdf
  
  myAggregate <- function(inParam="",inFun="mean",inRes=0){
    inFile <- paste0(climaticsDataDir,"VicInputGrids/Vic.soil.",inParam,".grd")
    tmp <- crop(stack(inFile),columbiaMask)
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
    tmpGrid <- mask(get(inParam), runCell, maskvalue = 0)
    if(dim(tmpGrid)[3] == 1) plot(tmpGrid,colNA="black",main=inParam) else plot(tmpGrid[[1]],colNA="black",main=inParam)
    plot(initialRegion, add=T, border="grey40")
    tmpGrid[is.na(tmpGrid)] <- fillValue
    if(dim(tmpGrid)[3] == 1) ncvar_put(paramsNcdf,get(varIn),rotate(as.matrix(tmpGrid))) else {
      tmpArray <- as.array(tmpGrid)
      t1 <- rotate(as.array(tmpArray[,,1]))
      t2 <- rotate(as.array(tmpArray[,,2]))
      t3 <- rotate(as.array(tmpArray[,,3]))
      outArray <- abind(t1,t2,t3,along = 3)
      ncvar_put(paramsNcdf,get(varIn),outArray) 
      ncatt_put(paramsNcdf,get(varIn),"description","")
    }
  }
  nc_close(paramsNcdf)
  
  ##########
  #
  # Snow Parameters
  #
  ##########
  
  # Logic: if the resolution is greater than resolutionThreshold create snow bands
  # snow bands are limited to maxSnowBands
  # if the elevation difference is smaller than the elevDiffThreshold go with something smaller
  
  # NOTE: right now the routine does NOT have the ability to subset
  
  resolutionThreshold <- 5
  if(resolutionKm > resolutionThreshold) {
    # SNOW
    maxSnowBands <- 6
    ElevDiffThreshold <- 280 # this is the elevation difference based on the binning on whether to use maxSnowBands as the number of valid bins, or use something smaller if the 
    # nmbrSnowBands <- min(resolutionKm / 5, maxSnowBands)
    #
    inParam <- "elev"
    inFile <- paste0(climaticsDataDir,"VicInputGrids/Vic.soil.",inParam,".grd")
    elevHiRes <- mask(crop(stack(inFile),columbiaMask),initialRegion)
    plot(elevHiRes)
 
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
      par(mfrow=c(1,2))
      plot(elevHiRes)
      plot(grdPolys, border="blue",add = T)
      plot(sp, col="black", add=T, pch="+",cex=2)
      plot(initialRegion, add=T, border="grey40")
      color <- c(2,2,3,4,5) 
      color_transparent <- adjustcolor(color, alpha.f = 0.5) 
      plot(thisCell,col=color_transparent,add=T)

      thisElev <- mask(crop(elevHiRes,thisCell),thisCell)
      # plot(thisElev)
      # plot(grdPolys, border="blue",add = T)
      # plot(initialRegion, add=T, border="grey40")
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
    AreaFract <- snowParamMatrix[,(2 + maxSnowBands):(maxSnowBands*2 + 1)]
    Pfactor <- snowParamMatrix[,(2 + maxSnowBands*2):(maxSnowBands*3 + 1)]
    grdPolysCoords <- gCentroid(grdPolys,byid = T)
    elevation <-rasterFromXYZ(cbind(coordinates(grdPolysCoords),elevBands))
    plot(elevation)
    AreaFract <- rasterFromXYZ(cbind(coordinates(grdPolysCoords),AreaFract))
    plot(AreaFract)
    Pfactor <- rasterFromXYZ(cbind(coordinates(grdPolysCoords),Pfactor))
      
    dimSnowBand <- ncdim_def("snow_band","",1:maxSnowBands,longname = "snow band")
    var_cellnum <- ncvar_def("cellnum","N/A",list(dimLon,dimLat),prec = "double", longname = "")
    var_AreaFract <- ncvar_def("AreaFract","fraction",list(dimSnowBand,dimLon,dimLat),prec = "double")
    var_elevation <- ncvar_def("elevation","m",list(dimSnowBand,dimLon,dimLat),prec = "double", longname = "")
    var_Pfactor <- ncvar_def("Pfactor","fraction",list(dimSnowBand,dimLon,dimLat),prec = "double", longname = "")
    
    ncid_old <- nc_open( paramsNcdfFile, write=T)
    #ncnew <- nc_create(paste0(outNcdfFile,".nc"),ncid_old)
    
    ncid_old <- ncvar_add(ncid_old,var_cellnum)
    ncid_old <- ncvar_add(ncid_old,var_AreaFract)
    ncid_old <- ncvar_add(ncid_old,var_elevation)
    ncid_old <- ncvar_add(ncid_old,var_Pfactor)
    
    ncvar_put(ncid_old,"cellnum",rotate(as.matrix(gridcell)))
 
    AreaFract[is.na(AreaFract)] <- fillValue
    outArray <- stack2array(as.array(AreaFract))
    ncvar_put(ncid_old,"AreaFract",outArray)

    elevation[is.na(elevation)] <- fillValue
    outArray <- stack2array(as.array(elevation))
    ncvar_put(ncid_old,"elevation",outArray)
    
    Pfactor[is.na(Pfactor)] <- fillValue
    outArray <- stack2array(as.array(Pfactor))
    ncvar_put(ncid_old,"Pfactor",outArray)
    nc_close(ncid_old)
  }
  
  ##########
  #
  # Vegetation Parameters
  #
  ##########
  
  infile <- paste(rdsDataDir,"globalLandCoverV2.westernHemisphere.RData",sep='')
  landCover <- get(load(file=infile))
  thisLandCover <- mask(crop(landCover,columbiaMask),initialRegion)
  createLandCoverImage(inRaster = thisLandCover)
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
    par(mfrow=c(1,2))
    plot(thisLandCover)
    plot(grdPolys, border="blue",add = T)
    plot(sp, col="black", add=T, pch="+",cex=2)
    plot(initialRegion, add=T, border="grey40")
    color <- c(2,2,3,4,5) 
    color_transparent <- adjustcolor(color, alpha.f = 0.5) 
    plot(thisCell,col=color_transparent,add=T)
    
    thisVeg <- crop(thisLandCover,thisCell)
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
  

  
  # VEG
  dimVegClasses <- ncdim_def("veg_classes","",1:20,longname = "veg_classes")
  var_veg_class <- ncvar_def("veg_class","N/A",dimVegClasses,prec = "integer", longname = "veg_class")
  var_veg_descr <- ncvar_def("descr","N/A",dimVegClasses,prec = "character", longname = "")
  var_root_zone <- ncvar_def("root_zone","",list(dimLon,dimLat),prec = "integer", longname = "")
  var_month <- ncvar_def("month","",list(dimLon,dimLat),prec = "integer", longname = "")
  var_Nveg <- ncvar_def("nveg","",list(dimLon,dimLat),prec = "", longname = "")
  var_Cv <- ncvar_def("Cv","",list(dimVegClasses,dimLon,dimLat),prec = "", longname = "")
  var_root_depth <- ncvar_def("root_depth","",list(dimLon,dimLat),prec = "", longname = "")
  var_root_fract <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_LAI <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_overstory <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_rarc <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_rmin <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_wind_h <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_RGL <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_rad_atten <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_wind_atten <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_trunk_ratio <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_albedo <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_veg_rough <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  var_displacement <- ncvar_def("","",list(dimLon,dimLat),prec = "", longname = "")
  
  
  
  ncatt_put(ncout,0,"title",title$value)
  ncatt_put(ncout,0,"institution",institution$value)
  ncatt_put(ncout,0,"source",datasource$value)
  ncatt_put(ncout,0,"references",references$value)
  history <- paste("P.J. Bartlein", date(), sep=", ")
  ncatt_put(ncout,0,"history",history)
  ncatt_put(ncout,0,"Conventions",Conventions$value)
  
  
  
  
  
  
  
  
  
