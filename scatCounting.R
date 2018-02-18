# Here begins the script intended to verify whether counting scats using a modified Jolly-Seber model is feasible. 

# The steps will be 

# 1. Simulate scat deposition with a fixed rate \theta ('recruitment')
# 2. Simulate scat removals (i.e. encounters + 'survival') with p(detect) ~ track in grid cell `s`
# 3. Model scat encounters as a modified spatial Jolly-Seber model.
# 3a. Need to implement a per-capita recruitment rate, as data augmentation implies a varying recruitment rate. 

library(ggplot2)
library(dplyr)
library(foreach)

# Format example dog tracks ------------------------------------------------------------------------------------------------------------

library(rgdal)

gpxFiles = dir()[grep(pattern = '.gpx', x = dir(), perl = T)]

sites = c(rep("12B2",3), rep("15A4", 3))

gpxLayers = ogrListLayers(gpxFiles[1])

out = lapply(X = gpxFiles, FUN = function(x){readOGR(dsn = x, layer = 'track_points')})

names(out) = gpxFiles

allPoints = foreach(i = seq_along(out), .combine = rbind) %do% {
  data.frame(out[[i]]@coords, Site = sites[i], Date = out[[i]]@data$time %>% as.Date %>% unique) %>% rename(lon = coords.x1, lat = coords.x2)
}

allPoints$Date = as.factor(allPoints$Date)

# Need to normalize points to 0,1, but they need to be PER SITE, not all
# together. In this way, the centroids of the sites are centered on 0, and state
# space settings, scat populations apply to both.

allPoints[allPoints$Site == "12B2",c(1:2)] = allPoints %>% filter(Site == "12B2") %>% select(lon, lat) %>% scale
allPoints[allPoints$Site == "15A4",c(1:2)] = allPoints %>% filter(Site == "15A4") %>% select(lon, lat) %>% scale

# Set dates to "round" equivalents.



xlim = c(-3,3)
ylim = c(-2,2)

plot = ggplot(allPoints) + 
  geom_path(aes(x = lon, y = lat, color = Date)) + 
  coord_cartesian(xlim = xlim, ylim = ylim)
plot
  

# Simulation of scats & state space ----------------------------------------------------------------------------------------------------

NScat = 600

set.seed(1)

scatXY = cbind.data.frame(x = runif(n = NScat, min = xlim[1], max = xlim[2]),
                          y = runif(n = NScat, min = ylim[1], max = ylim[2]))

plot + 
  geom_point(data = scatXY, aes(x = x, y = y), shape = 1)

# Simulate encounters of scats ------------------------------------------------------------------------------------------------------------

# Analyze encounters using JAGS ------------------------------------------------------------------------------------------------------------