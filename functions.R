require(rgdal)
require(doParallel)
require(dplyr)
require(rgeos)

# Setup --------------------------


importGPX = function(startPath = NULL, outPath){
  
  # Function goes to the archive and copies over all of the .gpx files from
  # COMPLETE paths. THEre are a bunch of incomplete paths, and they seem well
  # marked in their respective folders. I will iterate through the directories
  # and subdirectories and grab only those gpx files in non-incomplete folders.
  
  # Get folder structure to iterate through. ONly to be used on work computer
  # due to hard location of archived data.
  
  if(is.null(startPath)){struct = 
    dir('D:/MooseArchive/Field Data 2017 (Archive)/CK9Moose2017_Cornell-20180218T030010Z-001/CK9Moose2017_Cornell/Tracklogs/', 
        pattern = '.gpx',
        all.files = T, 
        recursive = T, 
        include.dirs = T,
        full.names = T)
  } else {
    struct = dir(startPath, 
                 pattern = '.gpx',
                 all.files = T, 
                 recursive = T, 
                 include.dirs = T,
                 full.names = T)
  }
  
  # Move files into trackLogs, ideally without the folder structure
  
  filesWant = grepl(pattern = "^((?!hum).)*(gpx)*$", x = struct, ignore.case = T, perl = T)
  
  file.copy(from = struct[filesWant], to = outPath, recursive = F, overwrite = F)
  
  return(NULL)
  
}

siteInfoFromFileName = function(path = NULL){
  # browser()
  if(is.null(path)){
    fileNames = list.files(pattern = '.gpx')
  }else{fileNames = list.files(path = path, pattern = '.gpx')}
  
  siteNames = fileNames %>% {regmatches(x = ., m = regexec(pattern = "\\d+[A-Z]\\d", text = ., perl = T))} %>% as.character()
  siteDates = fileNames %>% {regmatches(x = ., m = regexec(pattern = "\\d+\\.\\d+\\.\\d+", text = ., perl = T))} %>% as.character() %>% as.Date(format = '%m.%d.%y')
  siteHandler = fileNames %>% {regmatches(x = ., m = regexec(pattern = '(\\d{1,2}\\.{1}\\d{1,2}\\.{1}\\d{1,4}_)(\\w{2})', text = ., perl = T))} %>% 
    sapply(X = ., FUN = `[`, 3)
  
  error = rep("CLEAN", length(fileNames))
  error[fileNames %>% grepl(pattern = '^.*(inc).*$', x = ., perl = T)] = 'INCOMPLETE'
  error[fileNames %>% grepl(pattern = '^.*(shar).*$', x = ., perl = T)] = 'SHARED' # Track is supposedly split among two or more .gpx files. Only a few in 2017. 
  
  siteInfo = data.frame(siteID = siteNames, Date = siteDates, Handler = siteHandler, Error = error %>% factor)
  
}

# Reformats data dates properly and includes definition of sampling rounds.

reDate = function(data){
  
  data$Date = tryCatch(as.Date(data$Date, format = '%m/%d/%Y'),
                       error = function(m){as.Date(data$Date, origin = '1970-01-01')})
  data$Round = NA
  data$RoundNo = NA
  
  clearIndex = data$Date >= "2017-06-11" & data$Date <= "2017-07-05"
  samp1 = data$Date > "2017-07-05" & data$Date <= "2017-07-17"
  samp2 = data$Date > "2017-07-17" & data$Date <= "2017-07-28"
  samp3 = data$Date > "2017-07-28" & data$Date <= "2017-08-07"
  samp4 = data$Date > "2017-08-07" & data$Date <= "2017-08-17"
  samp5 = data$Date > "2017-08-17"
  
  data$Round[clearIndex] = "Clearing"
  data$Round[samp1] = "Sample1"
  data$Round[samp2] = "Sample2"
  data$Round[samp3] = "Sample3"
  data$Round[samp4] = "Sample4"
  data$Round[samp5] = "Sample5"
  
  data$RoundNo[clearIndex] = 0
  data$RoundNo[samp1] = 1
  data$RoundNo[samp2] = 2
  data$RoundNo[samp3] = 3
  data$RoundNo[samp4] = 4
  data$RoundNo[samp5] = 5
  
  # for(i in 1:nrow(data)){
  #   if(data$Date[i] >= "2017-06-11" & data$Date[i] <= "2017-07-05"){data$Round[i] = "Clearing"; data$RoundNo[i] = 0}
  #   if(data$Date[i] > "2017-07-05" & data$Date[i] <= "2017-07-17") {data$Round[i] = "Sample1";   data$RoundNo[i] = 1}
  #   if(data$Date[i] > "2017-07-17" & data$Date[i] <= "2017-07-28") {data$Round[i] = "Sample2";   data$RoundNo[i] = 2}
  #   if(data$Date[i] > "2017-07-28" & data$Date[i] <= "2017-08-07") {data$Round[i] = "Sample3";   data$RoundNo[i] = 3}
  #   if(data$Date[i] > "2017-08-07" & data$Date[i] <= "2017-08-17") {data$Round[i] = "Sample4";   data$RoundNo[i] = 4}
  #   if(data$Date[i] > "2017-08-17")                                {data$Round[i] = "Sample5";    data$RoundNo[i] = 5}
  # }
  
  return(data)
  
}


# Gets gpx files from path, and assigns some information to it such as site ID, dates, based on function siteInfoFromFileName.

getGPX = function(path = NULL, debug = F, debugLim = 10, siteInfo){
  
  if(is.null(path)){gpxFiles = dir()[grep(pattern = '.gpx', x = dir(), perl = T)]}else{
    gpxFiles = dir(path = path, full.names = T)[grep(pattern = '.gpx', x = dir(path = path), perl = T)]
  }
  
  if(debug){
    gpxFiles = gpxFiles[1:debugLim]
    siteInfo = siteInfo[1:debugLim,]
  }
  
  # gpxLayers = ogrListLayers(gpxFiles[1]) # Gets the layers that exist in the gpx object.
  
  out = lapply(X = gpxFiles, FUN = function(x){readOGR(dsn = x, layer = 'track_points')})
  browser()
  names(out) = paste0(siteInfo$siteID, "_", siteInfo$siteDate)
  
  # Convert to UTM for easy grid creation.
  
  out = lapply(X = out, FUN = function(x) spTransform(x = x, CRSobj = CRS("+proj=utm +zone=18 +datum=WGS84")))
  
  
  
  return(out)
  
}

# Point utilities ----------------------------------------------------------------------

# Takes gpx files, makes a SpatialPointsDataFrame from it. 
# Formats dates properly, adds 'Round', and 'RoundBySite' for ID purposes.

convertPoints = function(gpx, siteInfo){
  
  sites = siteInfo$siteID %>% as.character()
  
  allPoints = foreach(i = seq_along(gpx), .combine = rbind) %do% {
    
    # ONLY project the ones with an existing projection. Set the others later.
    if(proj4string(gpx[[i]]) != '+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' & !is.na(proj4string(gpx[[i]]))){
      gpx[[i]] = spTransform(x = gpx[[i]], CRSobj = CRS('+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
    }
    
    data.frame(gpx[[i]]@coords, 
               Site = siteInfo$siteID[i], 
               Date = ifelse(test = all(is.na(gpx[[i]]@data$time %>% as.Date %>% unique)), 
                             yes = siteInfo$Date[i] %>% as.Date(origin = '1970-01-01'),
                             no = gpx[[i]]@data$time %>% as.Date %>% unique
                             ),
               Time = tryCatch(expr = {gpx[[i]]@data$time %>% strptime(format = '%Y/%m/%d %T')},
                               error = function(m){return(NA)}
                               ),
               Handler = siteInfo$Handler[i],
               Error = siteInfo$Error[i]) %>% 
      rename(Easting = coords.x1, Northing = coords.x2)
  }
  
  # allPoints$Date = as.factor(allPoints$Date)
  
  # Set dates to "round" equivalents. I.e. round 1, round 2, etc.
  
  allPoints = reDate(allPoints)
  
  allPoints = allPoints %>% mutate(RoundBySite = paste0(Site, ".", Round))
  
  coordinates(allPoints) = ~Easting + Northing
  
  return(allPoints)
  
}

getScaledData = function(transectPoints, adj.bbox = 100, gridSize = 50){ # Obtain a scaled grid from the bounding box from a back-scaled allPoints dataset.
  
  # Should fix this later to be more general once we move to testing on more sites.
  
  scaled_12B2 = transectPoints %>% as.data.frame() %>% filter(Site == "12B2") %>% select(Easting, Northing) %>% scale
  scaled_15A4 = transectPoints %>% as.data.frame() %>% filter(Site == "15A4") %>% select(Easting, Northing) %>% scale
  
  scaled_12B2_center = attr(x = scaled_12B2, which = 'scaled:center')
  scaled_12B2_scale = attr(x = scaled_12B2, which = 'scaled:scale')
  
  scaled_15A4_center = attr(x = scaled_15A4, which = 'scaled:center')
  scaled_15A4_scale = attr(x = scaled_15A4, which = 'scaled:scale')
  
  transectPoints@coords[transectPoints@data$Site == "12B2"] = scaled_12B2
  transectPoints@coords[transectPoints@data$Site == "15A4"] = scaled_15A4
  
  meanScale = scaled_12B2_scale %>% bind_rows(scaled_15A4_scale) %>% colMeans
  
  # Back-transforming the scaled & centered points to obtain a bounding box surrounding all the points.
  
  allPoints_bt = apply(X = transectPoints@coords, MARGIN = 1, FUN = function(x){x*scaled_12B2_scale + scaled_12B2_center}) %>% t %>% data.frame
  
  coordinates(allPoints_bt) = ~Easting + Northing
  
  bbox = allPoints_bt@bbox
  
  bbox = bbox + c(rep(-adj.bbox, 2), rep(adj.bbox, 2)) 
  
  gridOverlay = expand.grid(seq(bbox[1,1], bbox[1,2], gridSize), seq(bbox[2,1], bbox[2,2], gridSize))
  
  scaledGrid = gridOverlay %>% scale(center = T, scale = meanScale) %>% as.data.frame %>% rename(Easting = Var1, Northing = Var2) %>% mutate(ID = seq(1:nrow(.)))
  
  return(list(scaledGrid = scaledGrid, scaledTracks = transectPoints))
  
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
  
  
  
  return(binxy) # Returns indexes of the x and y in the grid layer. Use next function to add the gridID's to the point layer.
}

addGridID_to_Points = function(queryPoints, refPointsToGrid_Output, gridLayer){
  
  gridX = gridLayer %>% pull(Easting) %>% unique %>% sort
  gridY = gridLayer %>% pull(Northing) %>% unique %>% sort
  
  # Adding to point layer
  queryPoints = queryPoints %>% mutate(gridID = data.frame(x = gridX[refPointsToGrid_Output$x], y = gridY[refPointsToGrid_Output$y]) %>% 
                                         left_join(y = gridLayer, by = c("x" = "Easting", "y" = "Northing")) %>% pull(ID))
  
  return(queryPoints)
  
}

trackFixesCount = function(track, gridLayer){
  
  # Count number of GPS fixes within a grid cell, proxy to time*space interaction. Fixes are at 1 minute intervals *nearly always*. 
  
  refOut = refPointsToGrid(track, gridLayer)
  
  track = addGridID_to_Points(track, refPointsToGrid_Output = refOut, gridLayer = gridLayer)
  
  # Now, count up the first, second, etc. occurrences of each grid cell. 
  
  fixesPerCell = track$gridID %>% table %>% as.data.frame() %>% rename(gridID = '.') %>% mutate(gridID = as.integer(as.character(gridID)))
  
  return(fixesPerCell)
  
}

points2line = function(points, ident){
  
  if(class(points) != "SpatialPointsDataFrame"){stop("Points layer must be SpatialPointsDataFrame")}
  
  layer_split = split(points, points@data[,ident])
  
  layer_lines = lapply(layer_split, function(x){Lines(list(Line(coordinates(x)[,c(1,2)])), x@data[,ident][1])})
  
  layer_lines = SpatialLines(layer_lines)
  
  # browser()
  
  ref = rle(points@data$RoundBySite)
  
  sites = points@data %>% group_by_(eval(ident)) %>% summarise(Site = first(Site))
  dates = points@data %>% group_by_(eval(ident)) %>% summarise(Date = first(Date))
  rounds = points@data %>% group_by_(eval(ident)) %>% summarise(Round = first(Round))
  roundsNo = points@data %>% group_by_(eval(ident)) %>% summarise(RoundNo = first(RoundNo))
  
  data = data.frame(id = ref$values,
                    Site = sites$Site,
                    Date = dates$Date,
                    Round = rounds$Round,
                    RoundNo = roundsNo$RoundNo)
  
  rownames(data) = data$id
  
  out = SpatialLinesDataFrame(layer_lines, data)
  
  return(out)
  
}

# Simulation ----------------------------------------------------------------------------

simScats = function(transectPoints, scats_init = 500, gridLayer, siteToTest = "12B2", recruit_rate = 20, maxR = 3, debug = F, seed = 1, 
                    probForm = "indicator", p0 = 0.8, a_fixes = 3, b_fixes = 0.05){
  
  if(!probForm %in% c("indicator", "length", "constant", "fixes")) message("probForm not within parameters. Defaulting to indicator probability formulation.")
  
  if(p0 < 0 | p0 > 1){stop("p0 must be bounded by 0,1.")}
  
  # function for probForm
  
  # there are two forms for probForm - 'length' and 'indicator'. 
  
  # Form indicator
  
  # p(detect) = p0 * I(track-in-grid)
  
  # Form length
  
  # p(detect) = p0 * (len(track-in-grid))
  
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
  
  bbox_scaled = getBbox(gridLayer %>% select(Easting, Northing))
  
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
  
  gridsVisited.rec = list()
  
  for(r in 1:maxR){
    
    siteTrackPoints = transectPoints %>% as.data.frame %>% filter(Site == siteToTest, Round == r)
    
    gridsVisited = countPointsInGrid(queryPoints = siteTrackPoints, gridPoints = gridLayer %>% select(Easting, Northing)) %>% filter(Freq > 0)
    
    gridsVisitedID = cbind(x = gridX[gridsVisited$x], y = gridY[gridsVisited$y]) %>% as.data.frame %>% left_join(y = gridLayer, by = c("x" = "Easting", "y" = "Northing"))
    
    scatsAvail = scatXY %>% filter(Removed == 0, gridID %in% gridsVisitedID$ID)
    
    if(debug){ #show which scats are available to be detected
      
      refOut = refPointsToGrid(siteTrackPoints, gridLayer)
      siteTrackPoints = addGridID_to_Points(siteTrackPoints, refOut, gridLayer)
      
      print(
        ggplot() +
          geom_tile(data = gridLayer, aes(x = Easting, y = Northing), fill = 'white', color = 'black') +
          geom_text(data = gridLayer, aes(x = Easting, y = Northing, label = ID), size = 3) +
          geom_path(data = siteTrackPoints, aes(x = Easting, y = Northing), color = 'blue') +
          geom_point(data = siteTrackPoints %>% filter(gridID == 130), aes(x = Easting, y = Northing), color = 'red')+
          geom_point(data = scatsAvail, aes(x = x, y = y)) + coord_map()
      )
      
    }
    
    
    # Probability models here. 
    
    if(probForm == "constant"){
      scatXY[scatXY$ID %in% scatsAvail$ID, "Removed"] = 1
    }
    
    if(probForm == "indicator"){
      
      scatsAvail$pEnc = p0
      
      scatXY[scatXY$ID %in% scatsAvail$ID,"pEnc"] = p0
      
      # Simulate encounters. All scats encountered will be removed.
      scatXY[scatXY$ID %in% scatsAvail$ID,"Removed"] = rbinom(n = nrow(scatsAvail), size = 1, prob = p0)
      
    }
    
    if(probForm == 'fixes'){
      
      gridFixes = trackFixesCount(siteTrackPoints, gridLayer)
      
      gridsVisitedID = gridsVisitedID %>% left_join(gridFixes, by = c("ID" = "gridID"))
      
      if(debug){
        
        gridLayer2 = gridLayer %>% left_join(gridFixes, by = c('ID' = 'gridID'))
        
        print(
          ggplot() +
            geom_tile(data = gridLayer2, aes(x = Easting, y = Northing, fill = Freq), color = 'black') +
            geom_text(data = gridLayer2, aes(x = Easting, y = Northing, label = ID), size = 3) + 
            geom_path(data = siteTrackPoints, aes(x = Easting, y = Northing), color = 'blue')
        )
        
      }
      
      scatsAvail = scatsAvail %>% left_join(gridFixes, by = c("gridID" = "gridID"))
      
      scatsAvail$Freq[is.na(scatsAvail$Freq)] = 0
      
      scatsAvail = scatsAvail %>% mutate(pEnc = (1 / (1 + exp((a_fixes - b_fixes*Freq)))))
      
      scatXY[scatXY$ID %in% scatsAvail$ID, "pEnc"] = scatsAvail$pEnc
      
      scatXY[scatXY$ID %in% scatsAvail$ID,"Removed"] = rbinom(n = nrow(scatsAvail), size = 1, prob = scatsAvail$pEnc)
      
      # scatXY$Removed = factor(scatXY$Removed, levels = c(0,1))
      
    }
    
    if(probForm == "length"){
      
      # Need to measure length of track in grid. Add to gridVisitRecords
      
    }
    
    gridsVisited.rec[[r]] = gridsVisitedID
    
    scatXY = scatXY %>% mutate(RoundRemoved = ifelse(test = {ID %in% scatsAvail$ID & Removed == 1}, yes = r, no = RoundRemoved)) %>% mutate(RoundRemoved = factor(RoundRemoved, levels = c(NA,seq(maxR))))
    
    if(debug){ # This plot displays grid ID's, and points characterized by Removal status, and RoundDeposited.
      print(
        ggplot() +
          geom_tile(data = gridLayer, aes(x = Easting, y = Northing), fill = 'white', color = 'black') +
          geom_text(data = gridLayer, aes(x = Easting, y = Northing, label = ID), size = 3) +
          geom_path(data = siteTrackPoints, aes(x = Easting, y = Northing), color = 'blue') +
          geom_point(data = scatXY %>% filter(ID %in% scatsAvail$ID), aes(x = x, y = y, shape = Removed, color = RoundDeposited, size = Removed)) + 
          scale_shape_manual(values = c(16,1)) + scale_size_manual(values = c(1,5)) +
          coord_map()
      )
      
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
  
  toReturn = list("ScatRecords" = scatXY.rec, "DepositionRecords" = depositionLog, "GridVisitsRecords" = gridsVisited.rec)
  
  return(toReturn)
  
}

# New ideas ----------------------------------------

simScats_simple = function(gridsVisited, scats_avg = 5, propDup = 0.5, scats_recruit = 1, maxR = 3, posVis = c(1,2), p0 = 0.5){
  
  # OOOOOK need to redo.
  # N array of dim c(nrow(gridsVisited), maxR+1, maxV)
  # y array of dim c(nrow(gridsVisited), maxR, maxV)
  # theta matrix of dim ; c(nrow(gridsVisited), maxR)
  
  # browser()
  
  # Variables
  numSites = nrow(gridsVisited)
  maxV = max(posVis)
  
  # Containers
  N = array(NA, dim = c(numSites, maxR+1, maxV))
  vis = y = array(0, dim = c(numSites, maxR+1, maxV))
  theta = matrix(0, nrow = numSites, ncol = maxR+1)
  
  # Initial deposition
  
  N[,1,] = rpois(n = numSites, lambda = scats_avg)
  
  gridsVisitRec = list()
  gridsVisitRec[[1]] = gridsVisited
  
  for(r in 1:maxR){
    
    # Roll for visit counts
    gridsVisited$nVis = sample(x = posVis, size = numSites, replace = T, prob = c(1-propDup, propDup))
    
    # Sum up scat population from previous pop, extra deposition, and removals
    N[,r+1,1] = N[,r,maxV] + theta[,r] - y[,r,maxV] # looks good
    
    sample_df = list()
    
    # Each pass, take some portion of scats.
    
    for(v in 1:maxV){
      
      gridsSampled = gridsVisited %>% filter(nVis >= v)
      
      index = gridsVisited$ID %in% gridsSampled$ID
      
      vis[index,r+1,v] = 1
      
      sample = rbinom(n = gridsSampled %>% nrow, size = N[index,r+1,v], prob = p0)
      
      y[index,r+1,v] = sample
      
      if(v < maxV){
        N[index,r+1,v+1] = N[index,r+1,v] - sample
        N[!index,r+1,v+1] = N[!index,r+1,v]
      }
      
      
      
    }
    
    theta[,r+1] = rpois(n = numSites, lambda = scats_recruit)
    
    gridsVisitRec[[r+1]] = gridsVisited
    
    
  }
  
  return(list("y" = y, "Deposition" = theta, "N" = N, 'vis' = vis, "GridVisitData" = gridsVisitRec))
  
}

runFunc = function(comboSet, iteration, gridsVisited){
  
  # Given a set of combinations of settings, pass each row in here to simulate and analyze.
  
  nSites = nrow(gridsVisited)
  
  # Number rounds
  maxR = 3
  
  # Possible visit counts; default is either visit once or twice. 
  posVis = c(1,2)
  
  # Max number of visits
  
  maxV = max(posVis)
  
  
  # # # # Simulation # # # #
  data = simScats_simple(gridsVisited = gridsVisited, scats_avg = combos$lam[iteration], propDup = combos$probdup[iteration], maxR = maxR, posVis = posVis, p0 = 0.5)
  # # # # # # # # # # # # # 
  
  # Pull out data
  y = data$y
  theta = data$Deposition
  vis = data$vis
  gridVisitData = data$GridVisitData
  N = data$N
  
  popAvail = N[,,1] %>% colSums()
  
  maxT = maxR + 1
  
  inits = function(){list(R = cbind(rep(NA,nSites), matrix(data = 1, nrow = nSites, ncol = maxR)),
                          N1 = rowSums(y))}
  
  data = list(y = y, vis = vis, nSites = nSites, maxT = maxT, maxV = maxV)
  
  params = c("N_time", 'p00', "theta", "lambda")
  
  niter = 1e5
  nburn = niter/4
  
  jagsOut = jags(data = data, inits = inits, parameters.to.save = params, model.file = 'model.txt', n.chains = 4, n.iter = niter, n.burnin = nburn, parallel = T)
  
  relevantData = jagsOut$summary
  
  allOut = list("N" = N,
                "popAvail" = popAvail,
                "y" = y,
                "vis" = vis,
                "output" = relevantData)
  
  save('allOut', file = paste0('jagsOut/out_', iteration,'.Rdata'))
  
  
}

scatDists = function(scatPos, trackLocs){
  
  # Want to use rgeos::gDistance, but limit search to tracks in same day & site as the scat. 
  
  site = scatPos$Site
  date = scatPos$Date
  
  tpoints_local = trackLocs[trackLocs$Site == site & trackLocs$Date == date,]
  
  return(rgeos::gDistance(scatPos, tpoints_local))
  
}


# Output analysis -----------------------

extract = function(what){invisible(Map(f = function(x,y){assign(x = x, value = y, pos = 1)}, x = names(what), y = what))}

assessOutputs = function(f, file.list, analysisIDs, settings){
  
  load(file.list[f])
  
  analysisID = analysisIDs[f]
  
  extract(what = settings %>% filter(ID == analysisID) %>% select(probdup : p))
  
  extract(allOut)
  
  # Estimate bias
  
  N_diff = (output[1:4] - popAvail) / popAvail
  
  p_diff = (output[5] - p) / p
  
  theta_diff = (output[6] - 1)
  
  lambda_diff = (output[7] - lam) / lam
  
  # Estimate SD 
  
  sd = output[9:15]
  
  return(data.frame(Param = rownames(output)[1:7],
                    Est = output[1:7],
                    Diff = c(N_diff, p_diff, theta_diff, lambda_diff),
                    SD = sd,
                    lower95 = output[1:7,3],
                    upper95 = output[1:7,7],
                    ID = analysisID))
  
}