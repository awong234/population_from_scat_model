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

scaledGrid = getGrid()

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