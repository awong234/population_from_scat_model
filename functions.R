
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
    data.frame(out[[i]]@coords, 
               Site = sites[i], 
               Date = out[[i]]@data$time %>% as.Date %>% unique, 
               Time = out[[i]]@data$time %>% strptime(format = '%Y/%m/%d %T')) %>% 
      rename(Easting = coords.x1, Northing = coords.x2)
  }
  
  allPoints$Date = as.factor(allPoints$Date)
  
  # Set dates to "round" equivalents. I.e. round 1, round 2, etc.
  
  allPoints = allPoints %>% mutate(Round = factor(ifelse(test = Date %in% c("2017-07-07", "2017-07-17"), yes = 1, 
                                                         ifelse(test = Date %in% c("2017-07-30", "2017-07-28"), yes = 2, no = 3))),
                                   RoundBySite = interaction(Site, Round))
  
  coordinates(allPoints) = ~Easting + Northing
  
  return(allPoints)
  
}

getScaledGrid = function(adj.bbox = 100){ # Obtain a scaled grid from the bounding box from a back-scaled allPoints dataset.
  
  allPoints_bt = apply(X = allPoints@coords, MARGIN = 1, FUN = function(x){x*scaled_12B2_scale + scaled_12B2_center}) %>% t %>% data.frame
  
  coordinates(allPoints_bt) = ~Easting + Northing
  
  (bbox = allPoints_bt@bbox)
  
  bbox = bbox + c(rep(-adj.bbox, 2), rep(adj.bbox, 2))
  
  gridOverlay = expand.grid(seq(bbox[1,1], bbox[1,2], 50), seq(bbox[2,1], bbox[2,2], 50))
  
  scaledGrid = gridOverlay %>% scale(center = T, scale = meanScale) %>% as.data.frame %>% rename(Easting = Var1, Northing = Var2)
  
  return(scaledGrid)
  
}

getBbox = function(grid){
  
  mins = c(min(grid[,1]), min(grid[,2]))
  maxs = c(max(grid[,1]), max(grid[,2]))
  
  data = cbind(mins, maxs)
  
  return(data)
  
}

countPointsInGrid = function(queryPoints, gridPoints){ 
  
  # May use for counting scats within grid cells, but also tracklog points.
  # Mostly good for finding 0's.
  
  gridX = gridPoints[,1] %>% unique %>% sort
  gridY = gridPoints[,2] %>% unique %>% sort
  
  binxy = data.frame(x = findInterval(queryPoints[,1], gridX - 0.5*(diff(gridX) %>% mean)),
                     y = findInterval(queryPoints[,2], gridY - 0.5*(diff(gridY) %>% mean)))
  
  results = table(binxy)
  
  d = as.data.frame.table(results) %>% mutate(x = as.character(x) %>% as.integer, y = as.character(y) %>% as.integer())
  
  return(d)
  
}

refPointsToGrid = function(queryPoints, gridPoints){
  
  # Offers references to grid points. Each row is a scat, with column 'x'
  # referring to the x'th slot in gridX, and column 'y' referring to the y'th.
  
  gridX = gridPoints[,1] %>% unique %>% sort
  gridY = gridPoints[,2] %>% unique %>% sort
  
  binxy = data.frame(x = findInterval(queryPoints[,1], gridX - 0.5*(diff(gridX) %>% mean)),
                     y = findInterval(queryPoints[,2], gridY - 0.5*(diff(gridY) %>% mean)))
  
  
  
  return(binxy)
}

addGridID_to_Points = function(queryPoints, refPointsToGrid_Output, gridLayer){
  
  gridX = gridLayer %>% pull(Easting) %>% unique %>% sort
  gridY = gridLayer %>% pull(Northing) %>% unique %>% sort
  
  queryPoints = queryPoints %>% mutate(gridID = data.frame(x = gridX[refPointsToGrid_Output$x], y = gridY[refPointsToGrid_Output$y]) %>% 
                                         left_join(y = gridLayer, by = c("x" = "Easting", "y" = "Northing")) %>% pull(ID))
  
  return(queryPoints)
  
}

simScats = function(scats_init = 500, gridLayer, siteToTest = "12B2", recruit_rate = 20, maxR = 3, debug = F, seed = 1){
  
  if(debug){browser()}
  
  # First, we will need to know which grids were searched at all. 
  
  # Then, we will need to know how much distance was covered within the grid cell,
  # and how much time was taken to cover it.
  
  # Then, we will need to tabulate the scats that are available for encounter, and
  # where they are.
  
  # Then, we calculate probability of encounter based on this metric, for each
  # scat within the grid cell. Those encountered are removed, and a new set of
  # 'recruited' scats are generated. Of course, they are independent of the
  # previous set, so it's likely just a matter of a new Poisson distributed population. 
  
  set.seed(seed)
  
  scatXY = cbind.data.frame(ID = 1:scats_init,
                            x = runif(n = scats_init, min = bbox_scaled[1,1], max = bbox_scaled[1,2]),
                            y = runif(n = scats_init, min = bbox_scaled[2,1], max = bbox_scaled[2,2]),
                            RoundDeposited = factor(x = 0, levels = c(0:maxR)), pEnc = 0, Removed = factor(x = 0, levels = c(0,1)),
                            RoundRemoved = factor(x = NA, levels = c(1:maxR)))
  
  gridX = gridLayer %>% pull(Easting) %>% unique %>% sort
  gridY = gridLayer %>% pull(Northing) %>% unique %>% sort
  
  
  scatsGridRef = refPointsToGrid(queryPoints = scatXY %>% select(x,y), gridPoints = gridLayer %>% select(Easting, Northing))
  
  scatXY = addGridID_to_Points(queryPoints = scatXY, refPointsToGrid_Output = scatsGridRef, gridLayer = gridLayer)
  
  scatXY.rec = list(scatXY)
  
  for(r in 1:maxR){
    
    siteTrackPoints = allPoints %>% as.data.frame %>% filter(Site == siteToTest, Round == r)
    
    gridsVisited = countPointsInGrid(queryPoints = siteTrackPoints, gridPoints = gridLayer %>% select(Easting, Northing)) %>% filter(Freq > 0)
    
    gridsVisitedID = cbind(x = gridX[gridsVisited$x], y = gridY[gridsVisited$y]) %>% as.data.frame %>% left_join(y = gridLayer, by = c("x" = "Easting", "y" = "Northing"))
    
    scatsAvail = scatXY %>% filter(Removed == 0, gridID %in% gridsVisitedID$ID)
    
    # To add later: a metric for encounter probability that models time and distance searched. For now, set as 1.
    scatsAvail$pEnc = 1 #later, some probability model ~ time, distance
    
    scatXY[scatXY$ID %in% scatsAvail$ID,"pEnc"] = 1
    
    # Simulate encounters. All scats encountered will be removed.
    scatXY[scatXY$ID %in% scatsAvail$ID,"Removed"] = rbinom(n = 1, size = 1, prob = scatsAvail$pEnc)
    
    scatXY = scatXY %>% mutate(RoundRemoved = ifelse(test = {ID %in% scatsAvail$ID & Removed == 1}, yes = r, no = RoundRemoved))
    
    if(debug){
      gplot = ggplot() +
        geom_tile(data = gridLayer, aes(x = Easting, y = Northing), fill = 'white', color = 'black') +
        geom_text(data = gridLayer, aes(x = Easting, y = Northing, label = ID)) +
        geom_path(data = siteTrackPoints, aes(x = Easting, y = Northing)) +
        geom_point(data = scatXY, aes(x = x, y = y, shape = Removed, color = RoundDeposited), size = 5) + 
        scale_shape_manual(values = c(16,1))
      gplot
    }
    
    # Simulate new scat deposition; no need on last round. 
    
    if(r < maxR){ 
      
      newScatN = rpois(n = 1, lambda = recruit_rate)
      
      if(r == 1) {
        depositionLog = vector(mode = 'integer', length = maxR)
        depositionLog[1] = scats_init
        depositionLog[2] = newScatN
      }
      
      depositionLog[r+1] = newScatN
      
      lastID = scatXY$ID[nrow(scatXY)]
      
      scatXY = scatXY %>% bind_rows(data.frame(ID = (lastID+1):(lastID+newScatN), 
                                               x = runif(n = newScatN, min = bbox_scaled[1,1], max = bbox_scaled[1,2]),
                                               y = runif(n = newScatN, min = bbox_scaled[2,1], max = bbox_scaled[2,2]),
                                               RoundDeposited = factor(r, levels = c(0:maxR)), pEnc = 0, Removed = factor(0, levels = c(0,1))))
      
      scatsGridRef = refPointsToGrid(queryPoints = scatXY %>% select(x,y), gridPoints = gridLayer %>% select(Easting, Northing))
      
      scatXY = addGridID_to_Points(queryPoints = scatXY, refPointsToGrid_Output = scatsGridRef, gridLayer = gridLayer)
      
    }
    
    
    
    # Validate grid id's
    
    if(debug){
      ggplot() +
        geom_tile(data = gridLayer, aes(x = Easting, y = Northing), fill = 'white', color = 'black') +
        geom_text(data = gridLayer, aes(x = Easting, y = Northing, label = ID)) +
        geom_path(data = siteTrackPoints, aes(x = Easting, y = Northing)) +
        geom_point(data = scatXY, aes(x = x, y = y, shape = Removed, color = RoundDeposited), size = 5) + 
        scale_shape_manual(values = c(16,1))
        
    }
    
    
    # Need a record of WHEN scats encountered. These are snapshots of the scat data from round to round.
    
    scatXY.rec[[r+1]] = scatXY
    
  }
  
  names(scatXY.rec) = paste("Round", 0:maxR)
  names(depositionLog) = paste("Round", 0:(maxR-1)) 
  
  toReturn = list("ScatRecords" = scatXY.rec, "DepositionRecords" = depositionLog)
  
  return(toReturn)
  
}