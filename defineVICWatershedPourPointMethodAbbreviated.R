defineVICWatershedPourPointMethodAbbreviated <- function(inBoundaryPoly=0,
                                                         lat = 0, lon = 0,
                                                         title = ""){
  
  # This is the "abbreviated" version of code to create a watershed boundary
  # It uses files that were created for the entire Columbia used in the routine
  # defineVICWatershedPourPointMethod.
  options(stringsAsFactors = F)
  source("/Volumes/MiniPro/Climatics/Code/Master/masterClimatics.R")
  library(sp)
  library(raster)
  library(maptools)
  library(rgdal)
  # requires saga "sudo port install saga"
  devtools::install_github("r-spatial/RSAGA", dependencies = TRUE)
  library(RSAGA)
  library(GISTools)
  library(ggmap)
  library(PBSmapping)
  source('/Volumes/MiniPro/Climatics/Code/Master/colorRamps.R')
  # set up environment for using RSAGA
  myenv <- rsaga.env()
  
  par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))
  long2UTM <- function(long) {
    (floor((long + 180)/6) %% 60) + 1
  }
  pourPointId <- gsub(" ","",title)
  outdir <- paste0(climaticsDataDir,"VICWatersheds/",pourPointId)
  if(!dir.exists(outdir)) dir.create(outdir)
  
  
  projUtmString <- "+proj=lcc +lat_1=54 +lat_2=40 +lon_0=-122.167 +ellps=WGS84"
  inBoundaryShape.utm = spTransform(inBoundaryPoly, CRS(projUtmString))
  
  poill <- data.frame(lat=thisSite$lat,lon=thisSite$lon)
  coordinates(poill) = ~lon + lat
  projection(poill) = CRS(projVals) 
  pod = spTransform(poill, CRSobj = CRS(projUtmString))
  
  careaFile <- paste0("/Volumes/MiniPro/Climatics/Data/VICWatersheds/ColumbiaMaster/carea_Columbia.sdat")
  SFda = raster(gsub("sgrd","sdat",careaFile))
  projection(SFda) = CRS(projUtmString)
  plot(log10(SFda))
  plot(pod, add = T, pch = 21, bg = "cyan")
  contour(fill_dem, col="grey50",add = T,lwd=.5)
  plot(inBoundaryShape.utm,add=T)
  
  bufferLength <- 1000 # meters
  buff = raster::extract(x = SFda, y = pod, buffer = bufferLength, cellnumbers = TRUE)
  buff.df = as.data.frame(buff)
  snap_cell <- buff.df$cell[which.max(buff.df$value)]
  snap_xy <- xyFromCell(SFda, snap_cell)
  snap_xy = data.frame(y = snap_xy[1, 2], x = snap_xy[1, 1])
  coordinates(snap_xy) = ~ x + y
  
  # determine pixels that drain into the poi - save as grid
  # rsaga.get.usage("ta_hydrology",module = 4)
  areaFile <-  paste0(outdir,"/basin_",pourPointId,".sgrd")
  demfilledfile <- "/Volumes/MiniPro/Climatics/Data/VICWatersheds/ColumbiaMaster/dem.utmfilled_Columbia.sgrd"
  rsaga.geoprocessor(
    lib = 'ta_hydrology', module = 4, env = myenv,
    param = list(TARGET_PT_X = snap_xy@coords[1],
                 TARGET_PT_Y = snap_xy@coords[2],
                 ELEVATION = demfilledfile,
                 AREA = areaFile,
                 METHOD = 0
    )					 
  )
  
  # convert grid version of catchment to a polygon
  rsaga.get.usage("shapes_grid",module = 6)
  shapeFile <- paste0(outdir,"/shape_",pourPointId,".shp")
  rsaga.geoprocessor(
    lib = 'shapes_grid', module = 6, env = myenv,
    param = list(GRID = areaFile,
                 POLYGONS = shapeFile,
                 CLASS_ALL = 0,
                 CLASS_ID = 100,
                 SPLIT = 0
    )
  )	
  
  # read in the shapefile as a spatial polygons data frame
  basin.spdf = readOGR(shapeFile)
  projection(basin.spdf) <- CRS(projUtmString)
  # convert the spdf to a spatial lines object 
  basin.sl = as(basin.spdf, "SpatialLines")
  
  # create full map with all elements
  bb = bbox(fill_dem)
  xminmax = bb[1, ]
  yminmax = bb[2, ]
  
  xmin = min(xminmax)
  xmax = max(xminmax)
  ymin = min(yminmax)
  ymax = max(yminmax)
  xsb = xmax + 500
  ysb = ymin + 700
  
  plotNetwork <- log10(SFda)
  plotNetwork[plotNetwork < 6.5,] <- NA
  par(mar = c(1, 1, 4, 4), cex = 1.1)
  plot(dem.utm, 
       col = gray(seq(1, 0.2, by = -0.1), alpha = 1),
       legend.mar = 5,
       legend.width = 0.6, legend.shrink = 0.4,
       legend.args = list(text = '  Elevation (m)', side = 3, font = 2, line = 1, cex = 0.8),
       main = paste0(title, " Watershed")
  )
  contour(fill_dem, col="brown",add = T,lwd=.5)
  plot(inBoundaryShape.utm,add=T, border="orange",lwd=2)
  plot(basin.sl, col = "green", lwd= 2,add = T) 
  plot(pod, add = T, pch = 21, bg = "cyan")
  plot(plotNetwork,col="blue",add=T, legend =F)

  t1 <- spTransform(basin.sl, CRS(projVals))
  t2 <- SpatialLines2PolySet(t1)
  basin.ll <- PolySet2SpatialPolygons(t2)
  raster::area(basin.ll) / 1e6
  demTrimmed <- crop(dem,basin.ll)
  c1 <- coordinates(gCentroid(basin.ll))
  outPolyfile <- paste0(outdir,"/shape_",pourPointId,".rds")
  saveRDS(basin.ll,outPolyfile)
  return(basin.ll)
}
  
