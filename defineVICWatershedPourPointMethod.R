defineVICWatershedPourPointMethod <- function(inBoundaryPoly = 0, pourPointCoords = 0, title = ""){
  # Reference for below code:
  # http://ibis.geog.ubc.ca/~rdmoore/rcode/ShannonFallsMap.r
  
  # http://ibis.geog.ubc.ca/~rdmoore/Rcode.htm
  # load libraries of required functions
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
  shapeExtent <- extent(inBoundaryPoly)
  lat1 <- ceiling(shapeExtent[4,]) + 1
  lat2 <- floor(shapeExtent[3,]) - 1
  lon0 <- mean(shapeExtent[1,],shapeExtent[2,])
  projUtmString <- paste0("+proj=lcc +lat_1=",lat1," +lat_2=",lat2," +lon_0=",lon0," +ellps=WGS84")
  
  
  # read in dem as a geotiff and reproject to utm square grid with 90 m resolution
  
  inFile <- paste0(climaticsDataDir,"Hydrosheds/mosaic_3s_NA_dem.raster")
  cropExtent <- extent(inBoundaryPoly)
  dem <- crop(raster(inFile),cropExtent)
  dem.utm = projectRaster(dem, crs = CRS(projUtmString), res = 90, method = "bilinear") 
  # Pour Point
  xy <- pourPointCoords
  pourPointCoords = xy
  poill <- pourPointCoords 
  coordinates(poill) = ~lon + lat
  projection(poill) = CRS(projVals) 
  pod = spTransform(poill, CRSobj = CRS(projUtmString))
  
  inBoundaryShape.utm = spTransform(inBoundaryPoly, CRS(projUtmString))

  # plot dem, tidbit locations and poi to check that everything looks okay
  plot(dem.utm, col = ter.cols)
  plot(inBoundaryShape.utm,add=T, border="grey40")
  plot(pod, add = T, pch = 21, bg = "cyan")
  textCoords <- coordinates(pod)
  text(textCoords[,1],textCoords[,2]+10000, title,col="blue")

  ########################################################################
  # delineate catchment area polygon using SAGA
  ########################################################################
  
  # write dem to SAGA format
  pourPointId <- gsub(" ","",title)
  outdir <- paste0(climaticsDataDir,"VICWatersheds/",pourPointId)
  if(!dir.exists(outdir)) dir.create(outdir)
  demfile <- paste0(outdir,"/dem.utm_",pourPointId,".sgrd")
  #flowdirfile <- paste0(outdir,"/flowdir_",pourPointId,".sgrd")
  writeRaster(dem.utm, filename = demfile, format = "SAGA", overwrite = T)
  
  # fill sinks - note that the algorithm will create negative elevations in some areas
  demfilledfile <- paste0(outdir,"/dem.utmfilled_",pourPointId,".sgrd")
  rsaga.fill.sinks(
    in.dem = demfile, out.dem = demfilledfile,
    method = "xxl.wang.liu.2006", minslope = 0.01,
    env = myenv
  )

  # read in filled dem, specify projection and convert negative values to NA
  gc()           
  fill_dem = raster(gsub("sgrd","sdat",demfilledfile))
  projection(fill_dem) <- CRS(projUtmString)
  values(fill_dem) = ifelse(values(fill_dem) < 0, NA, values(fill_dem))
  
  #writeRaster(fill_dem,demfilledfile,overwrite = T)
  
  # determine upslope drainage areas and map
  #rsaga.get.usage("ta_hydrology","Flow Accumulation (Top-Down)")
  gc()
  careaFile <- paste0(outdir,"/carea_",pourPointId,".sgrd")
  rsaga.topdown.processing(in.dem = demfilledfile, 
                          out.carea = careaFile, 
                           env = myenv)
  SFda = raster(gsub("sgrd","sdat",careaFile))
  projection(SFda) = CRS(projUtmString)
  plot(log10(SFda))
  plot(pod, add = T, pch = 21, bg = "cyan")
  contour(fill_dem, col="grey50",add = T,lwd=.5)
  plot(inBoundaryShape.utm,add=T)
  
  # find maximum drainage area near the poi - assign as catchment outlet
  # i'm using a buffer of 500 m, but a different size may be appropriate
  # for different data sets
  bufferLength <- 500 # meters
  buff = raster::extract(x = SFda, y = pod, buffer = bufferLength, cellnumbers = TRUE)
  buff.df = as.data.frame(buff)
  snap_cell <- buff.df$cell[which.max(buff.df$value)]
  snap_xy <- xyFromCell(SFda, snap_cell)
  snap_xy = data.frame(y = snap_xy[1, 2], x = snap_xy[1, 1])
  coordinates(snap_xy) = ~ x + y
  
  # determine pixels that drain into the poi - save as grid
  # rsaga.get.usage("ta_hydrology",module = 4)
  areaFile <-  paste0(outdir,"/basin_",pourPointId,".sgrd")
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
  #plot(inBoundaryShape.utm,add=T, border="orange",lwd=2)
  plot(basin.sl, col = "green", lwd= 2,add = T) 
  plot(pod, add = T, pch = 21, bg = "cyan")
  #plot(pod, add = T, pch = 21, bg = "yellow")
  plot(plotNetwork,col="blue",add=T, legend =F)
  contour(fill_dem, col="brown",add = T,lwd=.5)
  plot(inBoundaryShape.utm,add=T)
  
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

  
  