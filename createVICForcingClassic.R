createVICForcingClassic <- function(vicGrid = 0, yr = "",
                                    saveDir=""){
  
  # Read in a VIC domain netcdf file
  # Use Merra2 data to create the met forcing file
  #
  # k. westrick 1/18/2019
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICForcingMerra2.R"))
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,raster,VICmodel)
  forcingDir <- paste0(saveDir,"/Forcing/",yr,"/")
  if(!dir.exists(forcingDir)) dir.create(forcingDir,recursive = T)
  metList <- createVicForcingMerra2(vicGrid = vicGrid,yr = yr)
  for(x in 1:9){
    print(metList[x])
    print(range(getValues(metList[[x]])))
  }
  # resample from Merra2 native grid to the vic grid
  cat("resampling...",fill=T)
  metListResampled <- do.call(list,lapply(names(metList), function(inParam){
    return(raster::resample(metList[[inParam]],vicGrid))
  }))
  names(metListResampled) <- names(metList)
  coords <- trimws(format(coordinates(vicGrid),nsmall = 3))
  for(coordIdx in 1:dim(coords)[1]){
    cat("processing pixel ",coordIdx,", out of ",dim(coords)[1],fill=T)
    xy <- matrix(as.numeric(coords[coordIdx,]),nc=2)
    xyCh <- matrix(as.character(coords[coordIdx,]),nc=2)
    gridForcing <- do.call(cbind, lapply(names(metListResampled), function(inParam){
      # cat(inParam,fill = T)  
      return(as.vector(extract(metListResampled[[inParam]],xy)))
    }))
    forcingFileOrig <- paste0(forcingDir,"/forcings_orig_",
                          yr,"_",xy[,2],"_",xy[,1])
    retainCols <- !c(names(metListResampled) %in% c("PREC","TERHT"))
    write.table(gridForcing[,retainCols],forcingFileOrig,
                row.names = F,col.names = F,sep="\t")
    
    forcingFileBias <- paste0(forcingDir,"/forcings_biasCorr_",
                              yr,"_",xyCh[,2],"_",xyCh[,1])
    retainCols <- !c(names(metListResampled) %in% c("PRECBIASCORR","TERHT"))
    write.table(gridForcing[,retainCols],forcingFileBias,
                row.names = F,col.names = F, sep="\t")
  }
}
#createVICForcingClassic()
