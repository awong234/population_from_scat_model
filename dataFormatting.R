# Setup ------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

source('functions.R')

# Import data ------------------------------------------------------------

tracks2016 = rgdal::readOGR(dsn = 'gpxTracks2016_CLEANED', layer = 'tracks2016_clean')

tracks2016_points = rgdal::readOGR(dsn = 'gpxTracks2016_points_CLEANED', layer = 'tracks2016_points_clean', stringsAsFactors = F)

scats2016 = read.csv(file = 'scatLocs2016_cleaned.csv')

load('scatsData.Rdata')

attr(tracks2016_points@coords, 'dimnames') = list(NULL, c("Easting", "Northing"))

# Make grid ------------------------------------------------------------------------------------------------------------------------

# Make a grid over each of the tracks

grids = makeGrid(tracks2016, data = tracks2016@data)
siteOrder = tracks2016@data$Site %>% unique %>% as.character()

# Look at a grid and a corresponding track

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

# Get records of visitations to each grid cell, by ID, and the quantity of scats picked up -----------------------------------------------------------



# Format as needed by model ------------------------------------------------------------------------------------------------------------------------------------------------