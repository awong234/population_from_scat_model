# Here begins the script intended to verify whether counting scats using a modified Jolly-Seber model is feasible. 

# The steps will be 

# 1. Simulate scat deposition with a fixed rate \theta ('recruitment')
# 2. Simulate scat removals (i.e. encounters + 'survival') with p(detect) ~ track in grid cell `s`
# 3. Model scat encounters as a modified spatial Jolly-Seber model.
# 3a. Need to implement a per-capita recruitment rate, as data augmentation implies a varying recruitment rate. 

library(ggplot2)
library(dplyr)
library(foreach)
library(broom)

# Format example dog tracks ------------------------------------------------------------------------------------------------------------

library(rgdal)

sites = c(rep("12B2",3), rep("15A4", 3))

getGPX = function(){
  
  gpxFiles = dir()[grep(pattern = '.gpx', x = dir(), perl = T)]
  
  gpxLayers = ogrListLayers(gpxFiles[1])
  
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
  
  # Set dates to "round" equivalents.
  
  allPoints = allPoints %>% mutate(Round = factor(ifelse(test = Date %in% c("2017-07-07", "2017-07-17"), yes = "Round 1", 
                                                         ifelse(test = Date %in% c("2017-07-30", "2017-07-28"), yes = "Round 2", no = "Round 3"))),
                                   RoundBySite = interaction(Site, Round))
  
  coordinates(allPoints) = ~Easting + Northing
  
  return(allPoints)
  
}

out = getGPX()

allPoints = convertPoints()


# Need to normalize points to 0,1, but they need to be PER SITE, not all
# together. In this way, the centroids of the sites are centered on 0, and state
# space settings, scat populations apply to both.

# Scale individual transect points. Obtain scaling parameters.

scaled_12B2 = allPoints %>% as.data.frame() %>% filter(Site == "12B2") %>% select(Easting, Northing) %>% scale
scaled_15A4 = allPoints %>% as.data.frame() %>% filter(Site == "15A4") %>% select(Easting, Northing) %>% scale

scaled_12B2_center = attr(x = scaled_12B2, which = 'scaled:center')
scaled_12B2_scale = attr(x = scaled_12B2, which = 'scaled:scale')

scaled_15A4_center = attr(x = scaled_15A4, which = 'scaled:center')
scaled_15A4_scale = attr(x = scaled_15A4, which = 'scaled:scale')

allPoints@coords[allPoints@data$Site == "12B2"] = scaled_12B2
allPoints@coords[allPoints@data$Site == "15A4"] = scaled_15A4

meanScale = scaled_12B2_scale %>% bind_rows(scaled_15A4_scale) %>% colMeans


# Assign new coordinates to data frame.

xlim = c(-3,3)
ylim = c(-2,2)

# Generate grid

getGrid = function(adj.bbox = 100){
  
  allPoints_bt = apply(X = allPoints@coords, MARGIN = 1, FUN = function(x){x*scaled_12B2_scale + scaled_12B2_center}) %>% t %>% data.frame
  
  coordinates(allPoints_bt) = ~Easting + Northing
  
  (bbox = allPoints_bt@bbox)
  
  bbox = bbox + c(rep(-adj.bbox, 2), rep(adj.bbox, 2))
  
  gridOverlay = expand.grid(seq(bbox[1,1], bbox[1,2], 50), seq(bbox[2,1], bbox[2,2], 50))
  
  scaledGrid = gridOverlay  %>% scale(center = T, scale = meanScale) %>% as.data.frame %>% rename(Easting = Var1, Northing = Var2)
  
  return(scaledGrid)
  
}

scaledGrid = getGrid()

getBbox = function(grid){
  
  mins = c(min(grid[,1]), min(grid[,2]))
  maxs = c(max(grid[,1]), max(grid[,2]))
  
  data = cbind(mins, maxs)
  
  return(data)
  
}

bbox_scaled = getBbox(scaledGrid)

plotTransects = ggplot() + 
  geom_tile(data = scaledGrid, aes(x = Easting, y = Northing), fill = 'white', color = 'black') + 
  geom_path(data = allPoints %>% as.data.frame(), aes(x = Easting, y = Northing, color = RoundBySite)) + 
  coord_cartesian(xlim = xlim, ylim = ylim) + coord_map()
plotTransects

# Simulation of scats & state space ----------------------------------------------------------------------------------------------------

NScat = 600

set.seed(1)

scatXY = cbind.data.frame(x = runif(n = NScat, min = bbox_scaled[1,1], max = bbox_scaled[1,2]),
                          y = runif(n = NScat, min = bbox_scaled[2,1], max = bbox_scaled[2,2]))

plotTransects + 
  geom_point(data = scatXY, aes(x = x, y = y), shape = 1)

# Simulate encounters of scats ------------------------------------------------------------------------------------------------------------


# Analyze encounters using JAGS ------------------------------------------------------------------------------------------------------------