createVICVegClassicNLDAS11 <- function(vicBdry = 0, vicGrid = 0, vicCellNmbr = 0,
                                       runCell = runCell, resolutionKm = NULL){
  # k. westrick 2/5/2019
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source(paste0(climaticsCodeDir,"Master/createLandCoverImage.R"))
  if (!require("pacman")) install.packages("pacman")
  
  pacman::p_load(ncdf4,raster)
  
  ##########
  #
  # Vegetation Parameters
  #
  ##########
  
  combinedVegTable <- read.csv(file=paste0(climaticsCodeDir,'Master/combinedVegTable.csv'))
  glcRemapTable <- read.csv(file=paste0(climaticsCodeDir,'Master/glcRemapTable.csv'))
  infile <- paste(rdsDataDir,"globalLandCoverV2.westernHemisphere.RData",sep='')
  landCover <- get(load(file=infile))
  thisLandCover <- crop(landCover,vicGrid)
  createLandCoverImage(inRaster = thisLandCover,inBoundary=vicBdry)
  landCvrLegendData <- read.table(file=paste0(climaticsCodeDir,'Master/legendColorsLabels.csv'),header = T,sep=",",skip = 1)
  
  rcl <- cbind(glcRemapTable$glcVegCode,glcRemapTable$remap)
  thisLandCoverReclass <- reclassify(thisLandCover,rcl)
  plot(thisLandCoverReclass)
  thisLandCover <- thisLandCoverReclass
  
  centerCoords <- coordinates(vicGrid)
  maskVals <- getValues(vicCellNmbr)
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
  #grdPolys <- grdPolys[over(sp, grdPolys),]
  nmbrCells <- length(grdPolys)
  processVeg <- function(cellNmbr){
    print(cellNmbr)
    thisCell <- grdPolys[cellNmbr,]
    inGrid <- runCell[cellNmbr]
    if(inGrid == 1){
      par(mfrow=c(1,1))
      plot(thisLandCover)
      plot(grdPolys, border="blue",add = T)
      #plot(sp, col="black", add=T, pch="+",cex=2)
      plot(vicBdry, add=T, border="grey40",lwd=2)
      color <- c(2,2,3,4,5) 
      color_transparent <- adjustcolor(color, alpha.f = 0.5) 
      plot(thisCell,col=color_transparent,add=T)
      
      thisVeg <- mask(crop(thisLandCover,thisCell),vicBdry)
      allVals <- getValues(thisVeg)
      validVals <- allVals[!is.na(allVals)]
      fractionInCell <- length(validVals) / length(allVals)
      counts <- plyr::count(validVals)
      fraction <- signif(counts$freq / sum(counts$freq),4)
      residual <- 1 - sum(fraction)
      if(residual != 0) fraction[which.max(fraction)] <- fraction[which.max(fraction)] + residual # true up for fractions to equal one
      veg <- matrix(c(seq(1,11),rep(NA,11)),nc=2,byrow=F)
      veg <- data.frame(veg_type=counts$x,area_fract=fraction)
      remapIndices <- glcRemapTable$remap[match(veg$veg_type,glcRemapTable$glcVegCode)]
      # root zones
      zone_1_depth <- combinedVegTable$root_depth_rz1[remapIndices]
      zone_1_fract <- combinedVegTable$root_fract_rz1[remapIndices]
      zone_2_depth <- combinedVegTable$root_depth_rz2[remapIndices]
      zone_2_fract <- combinedVegTable$root_fract_rz2[remapIndices]
      # LAI
      laiMax <- matrix(combinedVegTable$LAImax[glcRemapTable$remap][remapIndices],nc=1)
      laiMaxVector <- matrix(rep(laiMax,12),nc=12) 
      laiMin <- matrix(combinedVegTable$LAImin[glcRemapTable$remap][remapIndices],nc=1)
      laiMinVector <- matrix(rep(laiMin,12),nc=12) 
      # vector for the monthly values to use the summer (max) weighting, '1 - this' equals the winter (min) weighting
      #                    J, F, M,  A,   M, J, J, A, S, O,  N, D
      summerWeighting <- c(0, 0, 0, .25, .5, 1, 1, 1, 1, .5, 0, 0)
      summerWeightingVector <- matrix(rep(summerWeighting,length(remapIndices)),nc=12,byrow = T)
      LAI <- laiMaxVector * summerWeightingVector + laiMinVector * (1 - summerWeightingVector)
      colnames(LAI) <- paste0("LAI_",month.abb)
      
      veg <- matrix(c(seq(1,20),rep(NA,20)),nc=2,byrow=F)
      veg <- cbind(veg_type=counts$x,area_fract=fraction,
                   zone_1_depth,zone_1_fract,
                   zone_2_depth,zone_2_fract,
                   LAI)
      
    } else{
      #return a "missing" pixel profile
      veg <- matrix(rep(-99,18),nr=1)
    }
    return(veg)
  }
  veg <- do.call(list, lapply(1:nmbrCells, processVeg))
  names(veg) <- c(1:nmbrCells)
  #
  #
  # veglib
  processVeglib <- function(){
    Class <- matrix(seq(1,12,1),nc=1)
    OvrStry <- matrix(combinedVegTable$overstory,nc=1)
    Rarc <- matrix(combinedVegTable$rarc,nc=1)
    Rmin <- matrix(combinedVegTable$rmin,nc=1)
    WIND_H <- matrix(combinedVegTable$wind_h,nc=1)
    RGL <- matrix(combinedVegTable$RGL,nc=1)
    SolAtn <- matrix(combinedVegTable$rad_atten,nc=1)
    WndAtn <- matrix(combinedVegTable$wind_atten,nc=1)
    Trunk <- matrix(combinedVegTable$trunk_ratio,nc=1)
    Comments <- matrix(combinedVegTable$vegDescription,nc=1)
    # LAI
    laiMax <- matrix(combinedVegTable$LAImax,nc=1)
    laiMaxVector <- matrix(rep(laiMax,12),nc=12) 
    laiMin <- matrix(combinedVegTable$LAImin,nc=1)
    laiMinVector <- matrix(rep(laiMin,12),nc=12) 
    # vector for the monthly values to use the summer (max) weighting, '1 - this' equals the winter (min) weighting
    #                    J, F, M,  A,   M, J, J, A, S, O,  N, D
    summerWeighting <- c(0, 0, 0, .25, .5, 1, 1, 1, 1, .5, 0, 0)
    summerWeightingVector <- matrix(rep(summerWeighting,length(Rarc)),nc=12,byrow = T)
    LAI <- laiMaxVector * summerWeightingVector + laiMinVector * (1 - summerWeightingVector)
    colnames(LAI) <- paste0(toupper(month.abb),".LAI")
    # Albedo
    ALBmat <- matrix(combinedVegTable$albedo,nc=1)
    ALB <- matrix(rep(ALBmat,12),nc=12) 
    colnames(ALB) <- paste0(toupper(month.abb),".ALB")
    # Roughness
    ROUmat <- matrix(combinedVegTable$roughnessLengthM,nc=1)
    ROU <- matrix(rep(ROUmat,12),nc=12) 
    colnames(ROU) <- paste0(toupper(month.abb),".ROU")
    # Displacement
    DISmat <- matrix(combinedVegTable$displacement,nc=1)
    DIS <- matrix(rep(DISmat,12),nc=12) 
    colnames(DIS) <- paste0(toupper(month.abb),".DIS")
    veglib <- data.frame(Class, OvrStry, Rarc, Rmin, LAI, ALB,
                         ROU, DIS, WIND_H, RGL,
                         SolAtn, WndAtn, Trunk, 
                         Comments)
    return(veglib)
  }
  veglib <- processVeglib()
  
  outList <- list(veg,veglib)
}
