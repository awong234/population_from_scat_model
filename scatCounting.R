# Here begins the script intended to verify whether counting scats using a modified Jolly-Seber model is feasible. 

# The steps will be 

# 1. Simulate scat deposition with a fixed rate \theta ('recruitment')
# 2. Simulate scat removals (i.e. encounters + 'survival') with p(detect) ~ track in grid cell `s`
# 3. Model scat encounters as a modified spatial Jolly-Seber model.
# 3a. Need to implement a per-capita recruitment rate, as data augmentation implies a varying recruitment rate. 

library(ggplot2)
library(dplyr)

# Format example dog tracks ------------------------------------------------------------------------------------------------------------

library(rgdal)

gpxFiles = dir()[grep(pattern = '.gpx', x = dir(), perl = T)]

gpxLayers = ogrListLayers(gpxFiles[1])

out = lapply(X = gpxFiles, FUN = function(x){readOGR(dsn = x, layer = 'track_points')})

names(out) = c("12B2", "15A4")

points_15A4 = out$`15A4`@coords %>% data.frame %>% rename(lon = coords.x2, lat = coords.x1)
points_12B2 = out$`12B2`@coords %>% data.frame

ggplot() + 
  geom_point(data = points_15A4, aes(x = ))

# Need to normalize points to 0,1



# Simulation of scats & state space ----------------------------------------------------------------------------------------------------

# Simulate encounters of scats ------------------------------------------------------------------------------------------------------------

# Analyze encounters using JAGS ------------------------------------------------------------------------------------------------------------