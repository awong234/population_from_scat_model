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

source('functions.R')

# Format example dog tracks ------------------------------------------------------------------------------------------------------------

library(rgdal)

# These are the sites that the gpx files are pulled from. Can probably automate in the future if needed. 

# Press F2 on any function to see its contents.

sites = c(rep("12B2",3), rep("15A4", 3))

out = getGPX() #loads gpx files

allPoints = convertPoints() #takes gpx files and converts to a complete dataset with points, dates, sites, and 'rounds'

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

scaledGrid = getScaledGrid() %>% mutate(ID = 1:nrow(.)) %>% select(ID, Easting, Northing)

bbox_scaled = getBbox(scaledGrid %>% select(Easting, Northing))

ggplot() + 
  geom_tile(data = scaledGrid, aes(x = Easting, y = Northing), fill = 'white', color = 'black') + 
  geom_path(data = allPoints %>% as.data.frame(), aes(x = Easting, y = Northing, color = RoundBySite)) + 
  coord_cartesian(xlim = xlim, ylim = ylim) + coord_map()


# Simulation of scats & state space ----------------------------------------------------------------------------------------------------

set.seed(1)

scats_init = rpois(n = 1, lambda = 500)


scatXY = cbind.data.frame(ID = 1:scats_init,
                          x = runif(n = scats_init, min = bbox_scaled[1,1], max = bbox_scaled[1,2]),
                          y = runif(n = scats_init, min = bbox_scaled[2,1], max = bbox_scaled[2,2]),
                          Round = 0, pEnc = 0, Removed = 0)

# Try out binning scats by grid. 

gridX = scaledGrid$Easting %>% unique %>% sort
gridY = scaledGrid$Northing %>% unique %>% sort

d = countPointsInGrid(queryPoints = scatXY %>% select(x,y), gridPoints = scaledGrid %>% select(Easting, Northing))

scatsGridRef = refPointsToGrid(queryPoints = scatXY %>% select(x,y), gridPoints = scaledGrid %>% select(Easting, Northing))

scatXY = scatXY %>% mutate(gridID = data.frame(x = gridX[scatsGridRef$x], y = gridY[scatsGridRef$y]) %>% 
                             left_join(y = scaledGrid, by = c("x" = "Easting", "y" = "Northing")) %>% pull(ID))



ggplot() + 
  geom_tile(data = scaledGrid, aes(x = Easting, y = Northing), fill = 'white', color = 'black') + 
  geom_text(data = d, aes(x = gridX[x], y = gridY[y], label = Freq)) + 
  geom_path(data = allPoints %>% as.data.frame(), aes(x = Easting, y = Northing, color = RoundBySite)) + 
  geom_point(data = scatXY, aes(x = x, y = y), shape = 1) + 
  coord_cartesian(xlim = xlim, ylim = ylim) + coord_map()

# Simulate encounters of scats -------------------------------------------------------------------------------------------------------------

# First, we will need to know which grids were searched at all. 

gridSearched = countPointsInGrid(queryPoints = allPoints@coords, gridPoints = scaledGrid)

ggplot() + 
  geom_tile(data = scaledGrid, aes(x = Easting, y = Northing), fill = 'white', color = 'black') + 
  geom_text(data = gridSearched, aes(x = gridX[x], y = gridY[y], label = Freq)) + 
  geom_path(data = allPoints %>% as.data.frame(), aes(x = Easting, y = Northing, color = RoundBySite)) + 
  geom_point(data = scatXY, aes(x = x, y = y), shape = 1) + 
  coord_cartesian(xlim = xlim, ylim = ylim) + coord_map()




# Then, we will need to know how much distance was covered within the grid cell,
# and how much time was taken to cover it.

# Then, we will need to tabulate the scats that are available for encounter, and
# where they are.

# Then, we calculate probability of encounter based on this metric, for each
# scat within the grid cell. Those encountered are removed, and a new set of
# 'recruited' scats are generated. Of course, they are independent of the
# previous set, so it's likely just a matter of a new Poisson distributed population.

scatsAvail = scatXY %>% filter(Removed == 0)

scatSim = simScats(gridLayer = scaledGrid, scats_init = 500, recruit_rate = 20, maxR = 3, debug = F)


# Analyze encounters using JAGS ------------------------------------------------------------------------------------------------------------