createVICSoilClassic <- function(vicBdry = 0, vicGrid = 0,gridCell=0,
                                resolutionKm = NULL, nLayer = 3){
  # k. westrick 2/5/2019
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICForcingMerra2.R"))
  if(!requireNamespace("devtools")) install.packages("devtools")
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ncdf4,raster)
  vicGrid[is.nan(vicGrid)] <- NA
  # define and check for the availability of the required input files
  
  soilParams <- 
   data.frame(matrix(c(   
        #VIC name         levels      
        "infilt",         1, 
        "Ds",             1,
        "Dsmax",          1,
        "Ws" ,            1,
        "c" ,             1,
        "expt" ,          nLayer,
        "Ksat",           nLayer,
        "phi_s",          nLayer,
        "init_moist",     nLayer,
        "elev",           1,
        "depth",          nLayer,
        "avg_T",          1,
        "dp",             1,
        "bubble",         nLayer,
        "quartz",         nLayer,
        "bulk_density",   nLayer,
        "soil_density",   nLayer,
        "off_gmt",        1,
        "Wcr_FRACT",      nLayer,
        "Wpwp_FRACT",     nLayer,
        "rough",          1,
        "snow_rough",     1,
        "annual_prec",    1, 
        "resid_moist",    nLayer,
        "fs_active",      1,
        "July_Tavg",      1),
     nc=2,byrow = T))
  colnames(soilParams) <- c("vicName","dims")    
  defineCheckInput <- function(){
    l1 <- list.files(paste0(climaticsDataDir,"VicInputGrids/"),pattern = "^Vic.soil.*.grd")
    soilVicFileList <- matrix(unlist(strsplit(l1,"\\.")),ncol = 4, byrow = T)[,3]
    missingMask <- !soilParams$vicName[as.logical(soilParams$Include)] %in% soilVicFileList 
    missingSoilParams <- names(soilParams$vicRName[missingMask])
    if(length(missingSoilParams) == 0) cat("All base files are available\n") else {
      for(i in 1:length(missingSoilParams)){
        cat(paste0("missing ",missingSoilParams[i],"\n",fill=T))
      }
      stop("run processSoilForVic.R")
    }
    return(soilParams)
  }
  soilParams <- defineCheckInput()
  myAggregate <- function(inParam="",inFun="mean",inRes=0){
    inFile <- paste0(climaticsDataDir,"VicInputGrids/Vic.soil.",inParam,".grd")
    tmp <- crop(stack(inFile),vicBdry)
    if(inFun == "mean") outGrid <- aggregate(tmp,fact = inRes)
    if(inFun == "max") outGrid <- aggregate(tmp,fun = max, fact = inRes)
    outGrid[is.nan(outGrid)] <- NA
    return(outGrid)
  }
  soil <- do.call(cbind, lapply(soilParams$vicName, function(inParam){
      cat(inParam,fill = T)      
      varIn <- paste0("var_",inParam)
      assign(inParam, myAggregate(inParam,inRes=resolutionKm))
      tmpGrid <- crop(get(inParam), vicGrid, maskvalue = 0) 
      tmpGrid[is.na(tmpGrid)] <- -99
      if(dim(tmpGrid)[3] == 1) plot(tmpGrid,colNA="black",main=inParam) else plot(tmpGrid[[1]],colNA="black",main=inParam)
      plot(vicBdry, add=T, border="grey40")
      vicRName <- soilParams$vicRName[which(soilParams$vicName == inParam)]
      outFrame <- data.frame(getValues(tmpGrid))
      if(dim(tmpGrid)[3] > 1) outNames <- paste(inParam,1:3,sep = "_") else outNames <- inParam
      names(outFrame) <- outNames
      return(outFrame)
    }))
  soil <- cbind(data.frame(run_cell = getValues(vicGrid),
                     grid_cell = getValues(gridCell),
                     lat = coordinates(vicGrid)[,2],
                     lon = coordinates(vicGrid)[,1]),soil)
  return(soil)
}
