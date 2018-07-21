# Setup ------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

source('functions.R')

# Import data ------------------------------------------------------------

tracks2016 = rgdal::readOGR(dsn = 'gpxTracks2016_CLEANED', layer = 'gpxTracks2016_CLEANED')
scats2016 = read.csv(file = 'scatLocs2016_cleaned.csv')

tracks2016_points = as(object = tracks2016, "SpatialPointsDataFrame")

# Make a grid over each of the tracks

grids = makeGrid(tracks2016, data = tracks2016@data)

for(i in seq_along(tracks2016)){
  
  makeGrid(tracks2016@lines[[i]], data = tracks2016@data)
  
}