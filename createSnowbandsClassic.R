createSnowbandsClassic <- function(vicBdry = 0, vicGrid = 0, vicCellNmbr = 0,
                                resolutionKm = NULL){
  # k. westrick 2/5/2019
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(raster)
  
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
    grdPolys <- grdPolys[over(sp, grdPolys),]
    nmbrCells <- length(grdPolys)
    retrieveSnowBand <- function(cellNmbr){
      
      # This needs to be fixed, right now when the elevation difference is small it still defaults to nmbrSnowBands
      # this should actually go to something smaller, as we don't want to process the snow bands if they are not needed
      # see the commented out lines in the below "else" statement.
      print(cellNmbr)
      thisCell <- grdPolys[cellNmbr,]
      inGrid <- vicGrid[cellNmbr]
      inGrid <- vicGrid[cellNmbr]
      if(inGrid == 1){
        thisElev <- mask(crop(elevHiRes,thisCell),thisCell)
        plot(elevHiRes)
        plot(grdPolys, border="blue",add = T)
        plot(vicBdry, add=T, border="grey40")
        color <- c(2,2,3,4,5) 
        color_transparent <- adjustcolor(color, alpha.f = 0.5) 
        plot(thisCell,col=color_transparent,add=T)
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
                          freq=FALSE,col="orange", plot=F,
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
                          freq=FALSE,col="orange",plot=F,
                          main="Histogram",xlab="x",ylab="f(x)",yaxs="i",xaxs="i")
          binBracket <- (histOut$breaks[length(histOut$breaks)] - histOut$breaks[1]) / (maxSnowBands - 1)
          elevBands <- round(histOut$mids,0)
          areaFrac <- round(histOut$counts / sum(histOut$counts),3)
          residual <- 1 - sum(areaFrac)
          if(residual != 0) areaFrac[which.max(areaFrac)] <- areaFrac[which.max(areaFrac)] + residual # true up for fractions to equal one
          pFactor <- areaFrac
        } 
      cellMeanElev <- mean(getValues(thisElev),na.rm=T)
      } else {
        elevBands <- areaFrac <- pFactor <- rep(0,maxSnowBands)
      }
      outList <- list(cellNmbr, elevBands,areaFrac,pFactor)
    }
    snowParamList <- do.call(list, lapply(1:nmbrCells, retrieveSnowBand))
    nmbrCols <- length(unlist(snowParamList[1]))
    snowParamMatrix <- matrix(unlist(snowParamList),nc=nmbrCols,byrow = T)
  return(snowParamMatrix)
  }
}