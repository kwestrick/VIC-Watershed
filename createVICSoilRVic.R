createVICSoilRVic <- function(vicBdry = 0, vicGrid = 0, thisFilename = "Test",
                                resolutionKm = NULL){
  # k. westrick 2/5/2019
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  source(paste0(climaticsOpDir,"Preprocessing/VIC-Watershed/createVICForcingMerra2.R"))
  source(paste0(climaticsCodeDir,"Master/createLandCoverImage.R"))
  if(!requireNamespace("devtools")) install.packages("devtools")
  if (!require("pacman")) install.packages("pacman")
  
  pacman::p_load(ncdf4,raster)
  vicGrid[is.nan(vicGrid)] <- NA
  # define and check for the availability of the required input files
  
  soilParams <- 
   data.frame(matrix(c(   
        #VIC name         VIC-R name       Include      levels      separator
        "infilt",         "INFILT",        T,           1,          "",
        "Ds",             "Ds",            T,           1,          "",
        "Dsmax",          "Ds_MAX",        T,           1,          "",
        "Ws" ,            "Ws",            T,           1,          "",
        "c" ,             "C",             T,           1,          "",
        "expt" ,          "EXPT",          T,           3,          "_",
        "Ksat",           "Ksat",          T,           3,          "_",
        "phi_s",          "PHI",           T,           3,          "_",
        "init_moist",     "MOIST" ,        T,           3,          "_",
        "elev",           "ELEV" ,         T,           1,          "",
        "depth",          "DEPTH" ,        T,           3,          "_",
        "avg_T",          "AVG_T" ,        T,           1,          "",
        "dp",             "DP" ,           T,           1,          "",
        "bubble",         "BUBLE" ,        T,           3,          "",
        "quartz",         "QUARZ" ,        T,           3,          "",
        "bulk_density",   "BULKDN" ,       T,           3,          "",
        "soil_density",   "PARTDN" ,       T,           3,          "",
        "off_gmt",        "OFF_GMT" ,      T,           1,          "",
        "Wcr_FRACT",      "WcrFT" ,        T,           3,          "",
        "Wpwp_FRACT",     "WpFT" ,         T,           3,          "",
        "rough",          "Z0_SOIL" ,      T,           1,          "",
        "snow_rough",     "Z0_SNOW" ,      T,           1,          "",
        "annual_prec",    "PRCP",          T,           1,          "",
        "resid_moist",    "RESM" ,         T,           3,          "",
        "fs_active",      "FS_ACTV" ,      T,           1,          "",
        "July_Tavg",      "JULY_TAVG" ,    T,           1,          ""),
     nc=5,byrow = T))
  colnames(soilParams) <- c("vicName","vicRName","Include","dims","dimSepChar")    
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
  
  soil <- do.call(rbind, lapply(soilParams$vicName, function(inParam){
      cat(inParam,fill = T)      
      varIn <- paste0("var_",inParam)
      assign(inParam, myAggregate(inParam,inRes=resolutionKm))
      #tmpGrid <- mask(get(inParam), runCell, maskvalue = 0) * vicGrid
      tmpGrid <- crop(get(inParam), vicGrid, maskvalue = 0) 
      if(dim(tmpGrid)[3] == 1) plot(tmpGrid,colNA="black",main=inParam) else plot(tmpGrid[[1]],colNA="black",main=inParam)
      plot(vicBdry, add=T, border="grey40")
      vicRName <- soilParams$vicRName[which(soilParams$vicName == inParam)]
      if(dim(tmpGrid)[3] == 1) {
        outList <- matrix(tmpGrid,nc=ncell(tmpGrid))
        rownames(outList) <- vicRName
      } else {
        levSep <- soilParams$dimSepChar[which(soilParams$vicName == inParam)]
        outList <- matrix(tmpGrid,nc=ncell(tmpGrid),byrow = T)
        rownames(outList) <- paste(vicRName,1:3,sep=levSep)
      }
      return(outList)
    }))
  soil <-data.frame(t(soil))
  
return(soil)
}
