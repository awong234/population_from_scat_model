# Setup ------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

source('functions.R')

# Import data ------------------------------------------------------------

tracks2016 = rgdal::readOGR(dsn = 'gpxTracks2016_CLEANED', layer = 'gpxTracks2016_CLEANED')
scats2016 = read.csv(file = 'scatLocs2016_cleaned.csv')

tracks2016_points = as(object = tracks2016, "SpatialPointsDataFrame")
attr(tracks2016_points@coords, 'dimnames') = list(NULL, c("Easting", "Northing"))

# Make a grid over each of the tracks

grids = makeGrid(tracks2016, data = tracks2016@data)
siteOrder = tracks2016@data$Site %>% unique %>% as.character()

# Look at a grid and a corresponding track

siteIndexes = which(tracks2016@data$Site == '12B2')

i = 1

track_12B2 = tracks2016_points[which(tracks2016_points@data$Site == '12B2'),]

grid_12B2 = grids[[which(siteOrder == '12B2')]]

# Looks good.
ggplot() + 
  geom_tile(data = grid_12B2 %>% data.frame, aes(x = x, y = y), color = 'black', fill = NA) +
  geom_text(data = grid_12B2 %>% data.frame, aes(x = x, y = y, label = id), size = 1) +
  geom_point(data = track_12B2 %>% data.frame, aes(x = Easting, y = Northing)) + 
  coord_equal()

# Can we identify cells crossed by tracks?

?over

intersect = rle(over(x = track_12B2, y = grid_12B2)$id)

intersect$values

ggplot() + 
  geom_tile(data = grid_12B2 %>% data.frame, aes(x = x, y = y), fill = NA) +
  geom_tile(data = grid_12B2 %>% data.frame %>% filter(id %in% intersect$values), aes(x = x, y = y), color = 'red', fill = NA) +
  geom_text(data = grid_12B2 %>% data.frame, aes(x = x, y = y, label = id), size = 1) +
  geom_point(data = track_12B2 %>% data.frame, aes(x = Easting, y = Northing)) + 
  geom_point(data = scats2016 %>% filter(Site == '12B2'), aes(x = Easting, y = Northing), shape = 1, color = 'blue') + 
  coord_equal()

# Yes, looks good. We also get the order of visitation from intersect$values.

# No duplicated ID#'s, that's good
allIDs = foreach(i = 1:length(grids), .combine = c) %do% {
  grid = grids[[i]] %>% data.frame
  grid$id
}

# Shoot, we need the scats to be temporally referenced as well. Need to get a) scat collection time and b) nearest transect point in space & time, and assign that transect ID to it.
# scat data still temporally referenced, but track time data was lost in cleaning. Have to hope that the times were aligned.

load('trackPoints_2016_unclean.Rdata')

# For EVERY track point, we need to find the closest one, and then inherit the time from it. 

registerDoParallel(cores = detectCores())

tracks2016_points@data$Date_Time = NA

foreach(i = 1:nrow(tracks2016_points)) %do% {
  
  site = tracks2016_points[i,]$Site %>% as.character
  
  tracksFromSite = tracks_points[which(tracks_points$Site == site),]
  
  dists = rgeos::gDistance(tracks2016_points[i,], tracksFromSite, byid = T)
  
  index = which.min(dists)
  
  tracks2016_points@data$Date_Time[i] = tracksFromSite[index,]$Time
  
}
