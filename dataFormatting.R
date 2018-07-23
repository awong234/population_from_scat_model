# Setup ------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sp)

source('functions.R')

# Import data ------------------------------------------------------------

tracks2016 = rgdal::readOGR(dsn = 'gpxTracks2016_CLEANED', layer = 'tracks2016_clean', stringsAsFactors = F)

tracks2016_points = rgdal::readOGR(dsn = 'gpxTracks2016_points_CLEANED', layer = 'tracks2016_points_clean', stringsAsFactors = F)

scats2016 = read.csv(file = 'scatLocs2016_cleaned.csv', stringsAsFactors = F)

load('scatsData.Rdata')

attr(tracks2016_points@coords, 'dimnames') = list(NULL, c("Easting", "Northing"))

# Make grid ------------------------------------------------------------------------------------------------------------------------

# Make a grid over each of the tracks

grids = makeGrid(tracks2016, data = tracks2016@data)
siteOrder = tracks2016@data$Site %>% unique %>% as.character()

# SKIP

# Look at a grid and a corresponding track

skip = T

if(!skip){
  
  siteIndexes = which(tracks2016@data$Site == '12B2')
  
  track_12B2 = tracks2016_points[which(tracks2016_points@data$Site == '12B2'),]
  
  grid_12B2 = grids[[which(siteOrder == '12B2')]]
  
  # Looks good.
  ggplot() + 
    geom_tile(data = grid_12B2 %>% data.frame, aes(x = x, y = y), color = 'black', fill = NA) +
    geom_text(data = grid_12B2 %>% data.frame, aes(x = x, y = y, label = id), size = 1) +
    geom_point(data = track_12B2 %>% data.frame, aes(x = Easting, y = Northing)) + 
    coord_equal()
  
  # Can we identify cells crossed by tracks?
  
  intersect = rle(over(x = track_12B2, y = grid_12B2)$id)
  
  intersect$values
  
  ggplot() + 
    geom_tile(data = grid_12B2 %>% data.frame, aes(x = x, y = y), fill = NA) +
    geom_tile(data = grid_12B2 %>% data.frame %>% filter(id %in% intersect$values), aes(x = x, y = y), color = 'red', fill = NA) +
    geom_text(data = grid_12B2 %>% data.frame, aes(x = x, y = y, label = id), size = 1) +
    geom_point(data = track_12B2 %>% data.frame, aes(x = Easting, y = Northing)) + 
    geom_point(data = scatsReferenced %>% filter(Site == '12B2'), aes(x = Easting, y = Northing), shape = 1, color = 'blue') + 
    geom_point(data = scats2016 %>% filter(Site == '12B2'), aes(x = Easting, y = Northing), shape = 1, color = 'green') + 
    coord_equal()
  
  # Yes, looks good. We also get the order of visitation from intersect$values.
  
  # No duplicated ID#'s, that's good
  allIDs = foreach(i = 1:length(grids), .combine = c) %do% {
    grid = grids[[i]] %>% data.frame
    grid$id
  }
  
  identical(allIDs, seq_along(allIDs))
  
  rm(allIDs)
  
  # END SKIP
  
  
}


# Get records of visitations to each grid cell, by ID, and the quantity of scats picked up -----------------------------------------------------------

# Two things to do 

# FIRST - overall visitation to cells (using the track points data) -------------------------------------------------------

roundVisits = tracks2016_points@data %>% select(Site, RndBySt, Date) %>% unique %>% arrange(Date)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
rleTracks = rlePoints(visitInfo = roundVisits, points = tracks2016_points, grids = grids, plots = F)
names(rleTracks) = roundVisits$RndBySt
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# SECOND - find cells where scats were collected (using the referenced scat location data `scatsReferenced`) --------------------------------------------

scatsReferenced_spdf = scatsReferenced %>% arrange(Time)
coordinates(scatsReferenced_spdf) = ~Easting + Northing
proj4string(scatsReferenced_spdf) = proj4string(tracks2016)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
rleScats = rlePoints(visitInfo = roundVisits, points = scatsReferenced_spdf, grids = grids, plots = F, animation = F)
names(rleScats) = roundVisits$RndBySt
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# Get records of intervals between visits. For each site, what are the days between the visits?



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
diffDays = tracks2016 %>% data.frame %>% group_by(Site) %>% do(out = diff(.$Date %>% as.Date()))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# Some are only visited once, which are not useful for us. Which ones are these

whichSingleVisit = diffDays$out %>% lapply(FUN = length) %>% {which(. == 0)}

sites =  tracks2016 %>% data.frame %>% pull(Site) %>% unique

singleVisit = sites[whichSingleVisit]

# we don't want these. Save this for later.

# Time intervals - not consistent enough to use intersection$length as a covariate. Sometimes 1s, 5s, 30s, etc.

if(!skip){
  trackIntervals = tracks2016_points %>% data.frame %>% group_by(RndBySt) %>% do(out = diff(.$Time %>% as.POSIXct()))
  trackIntervals$out %>% lapply(FUN = function(x){summary(x %>% as.integer)})
}

# How many replicated observations? oof . . . 26 replicated visits out of 5

rleScats %>% lapply(FUN = function(x){x$values %>% table %>% as.integer}) %>% unlist %>% table

# Look at 266; 9A1.Sample1

grid_local = grids[[which(siteOrder == '09A3')]]
scats_local = scatsReferenced_spdf[scatsReferenced$RndBySt == '09A3.Sample1',]

rle(over(x = scats_local, y = grid_local)$id)

ggplot() + 
  geom_tile(data = grid_local %>% data.frame, aes(x = x, y = y), fill = NA, color = 'black') + 
  geom_text(data = grid_local %>% data.frame, aes(x = x, y = y, label = id), color = 'gray50', size = 3) + 
  geom_point(data = tracks2016_points %>% data.frame %>% filter(RndBySt == '09A3.Sample1'), aes(x = Easting, y = Northing), shape = 1) +
  geom_point(data = scats_local %>% data.frame, aes(x = Easting, y = Northing, color = 'red')) + 
  coord_equal(xlim = c(min(scats_local$Easting) + c(-100,100), max(scats_local$Easting) + c(-100,100)), 
              ylim = c(min(scats_local$Northing) + c(-100,100), max(scats_local$Northing) + c(-100,100)))
  


# Format as needed by model ------------------------------------------------------------------------------------------------------------------------------------------------

