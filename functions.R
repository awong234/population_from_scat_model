
getGPX = function(){
  
  gpxFiles = dir()[grep(pattern = '.gpx', x = dir(), perl = T)]
  
  gpxLayers = ogrListLayers(gpxFiles[1]) # Gets the layers that exist in the gpx object.
  
  if(!exists("out", where = .GlobalEnv)){
    out = lapply(X = gpxFiles, FUN = function(x){readOGR(dsn = x, layer = 'track_points')})
  }
  
  names(out) = gpxFiles
  
  # Convert to UTM for easy grid creation.
  
  out = lapply(X = out, FUN = function(x) spTransform(x = x, CRSobj = CRS("+proj=utm +zone=18 +datum=WGS84")))
  
  return(out)
}

convertPoints = function(){
  
  allPoints = foreach(i = seq_along(out), .combine = rbind) %do% {
    data.frame(out[[i]]@coords, Site = sites[i], Date = out[[i]]@data$time %>% as.Date %>% unique) %>% rename(Easting = coords.x1, Northing = coords.x2)
  }
  
  allPoints$Date = as.factor(allPoints$Date)
  
  # Set dates to "round" equivalents. I.e. round 1, round 2, etc.
  
  allPoints = allPoints %>% mutate(Round = factor(ifelse(test = Date %in% c("2017-07-07", "2017-07-17"), yes = "Round 1", 
                                                         ifelse(test = Date %in% c("2017-07-30", "2017-07-28"), yes = "Round 2", no = "Round 3"))),
                                   RoundBySite = interaction(Site, Round))
  
  coordinates(allPoints) = ~Easting + Northing
  
  return(allPoints)
  
}

getGrid = function(adj.bbox = 100){
  
  allPoints_bt = apply(X = allPoints@coords, MARGIN = 1, FUN = function(x){x*scaled_12B2_scale + scaled_12B2_center}) %>% t %>% data.frame
  
  coordinates(allPoints_bt) = ~Easting + Northing
  
  (bbox = allPoints_bt@bbox)
  
  bbox = bbox + c(rep(-adj.bbox, 2), rep(adj.bbox, 2))
  
  gridOverlay = expand.grid(seq(bbox[1,1], bbox[1,2], 50), seq(bbox[2,1], bbox[2,2], 50))
  
  scaledGrid = gridOverlay  %>% scale(center = T, scale = meanScale) %>% as.data.frame %>% rename(Easting = Var1, Northing = Var2)
  
  return(scaledGrid)
  
}

getBbox = function(grid){
  
  mins = c(min(grid[,1]), min(grid[,2]))
  maxs = c(max(grid[,1]), max(grid[,2]))
  
  data = cbind(mins, maxs)
  
  return(data)
  
}
