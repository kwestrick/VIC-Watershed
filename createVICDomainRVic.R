createVICDomainRVic <- function(vicBdry = 0, 
                            resolutionKm = 4.){
  
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
  #gridcell[gridcell == 1] <- seq(1,ncell(gridcell[gridcell == 1]))
  gridcell <- vicGrid
  gridcell<-setValues(vicGrid,seq(1,ncell(gridcell)))
  plot(gridcell,colNA="black")
  plot(vicBdry, add=T, border="grey40")
  
  # area
  gridArea <- raster::area(vicGrid) * 1e6
  totalVicArea <- raster::area(origMask) * 1e6
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
  outList <- list(vicGrid,vicArea,gridcell,runCell,vicFrac)
  names(outList) <- c("MASK","AREA","GRID","RUN","FRAC")
  return(outList)
}

