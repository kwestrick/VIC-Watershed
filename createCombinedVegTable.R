createCombinedVegTable <- function(){
  
  # Creates a GLC to VIC veg parameters lookup table
  
  # K. Westrick 2/7/2019
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  
  #
  vegParams <- c("veg_class","veg_descr","root_fract1","root_fract2","root_depth1", "root_depth2", "LAI",
                "overstory","RGL",	"rad_atten",	"wind_atten",	"rarc",	"rmin",	"wind_h",	
                "trunk_ratio","albedo","veg_rough","displacement")
  glcTable <-read.table(file=paste0(climaticsCodeDir,'Master/legendColorsLabels.csv'),header = T,sep=",",skip = 1)[,c(1,5)]
  colnames(glcTable) <- c("glcVegCode","glcVegDesc")
  #template[,3:length(vegParams)] <- NA
  #colnames(template) <- vegParams
  # openxlsx::write.xlsx(template,file=paste0(climaticsCodeDir,'Master/GLCtoVICLookupTable.xlsx'))
  
  # fill in some of the various categories using the NLDAS lookup avaialble from:
  # https://ldas.gsfc.nasa.gov/nldas/NLDASmapveg.php
      
  #
  createTables <- function(indices=0,colName=""){
    inTable <- openxlsx::read.xlsx(paste0(climaticsCodeDir,'Master/web.veg.table.xlsx'),
                                   colNames = F,
                                   rows = indices)[,1:2]
    inTable$X1 <- substring(inTable$X1,4)
    inTable$X1 <- gsub("^\\s+|\\s+$", "", inTable$X1)
    inTable$X2[inTable$X2 == "N/A"] <- NA
    inTable$X2[inTable$X2 == "NONE"] <- NA
    colnames(inTable) <- c("vegParam",colName)
    return(inTable)
  }
  # Canopy Top Meters
  vegTable <- createTables(c(7:20),"canopyTopM")[,2]
  
  # Roughness Length 
  vegTable <- cbind(vegTable, roughnessLengthM = createTables(c(45:58),"roughnessLengthM")[,2])

  # overstory 
  vegTable <- cbind(vegTable, overstory = createTables(c(501:514),"overstory")[,2])
  
  # displacement 
  vegTable <- cbind(vegTable, displacement = createTables(c(387:400),"displacement")[,2])
  
  # albedo 
  vegTable <- cbind(vegTable, albedo = createTables(c(463:476),"albedo")[,2])
  
  # displacement 
  vegTable <- cbind(vegTable, displacement = createTables(c(387:400),"displacement")[,2])
  
  # LAI max (summer)
  vegTable <- cbind(vegTable, LAImax = createTables(c(311:324),"LAImax")[,2])
  
  # LAI min (winter)
  vegTable <- cbind(vegTable, LAImin = createTables(c(330:343),"LAImin")[,2])
  
  # For other parameters use the table provided by Diana Gergel at the UW Computational Hydrology (7 Feb 2019)
  
  uwTable <- openxlsx::read.xlsx(paste0(climaticsCodeDir,'Master/vegParametersFromUWGergel.xlsx'))
  
  combinedVegTable <- cbind(uwTable,vegTable[2:13,],deparse.level = 1)

  write.csv(combinedVegTable,file=paste0(climaticsCodeDir,'Master/combinedVegTable.csv'))
  
  # create a mapping from glc indices to the combined veg table values, did this by hand
  # the ordering of the vector are the glcTable values, the actual vector values are what 
  # these values should remap to from the combined lookup
  
  glcTable
  
  #  glcTable values: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
  glcTable$remap <- c(2,4,1,3,5,7,8,7,9,9, 11,12,11,12,12,12,12,12,12,12)
  
  glcTable$glcMapsTo <- combinedVegTable$vegDescription[glcTable$remap]

  write.csv(glcTable,file=paste0(climaticsCodeDir,'Master/glcRemapTable.csv'))
}
