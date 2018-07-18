# Intention is to look at the collections and the tracks to see 

# A) how many grid cells are duplicated per track
# B) Of those, how many have collections on the first, and/or second duplication

# This is necessary to observe the potential of the pseudo-replication method.

# CLEANING TO BE DONE

# Some transects are mislabeled.
# At the moment, the RoundBySiteID is wrong. Fixing now.

source('functions.R')

library(dplyr)
library(rgdal)
library(sp)
library(ggplot2)

# Setup ------------------------------------------------------------------------------------------------------------------------

# Import gpx files that AREN't incomplete.
# importGPX() 

names = list.files(path = 'trackLogs/')

# Only reload all tracks if the .Rdata file somehow gets lost.

if(!"trackPoints.Rdata" %in% dir()){
  
  tracks = getGPX(path = 'trackLogs/', debug = F)
  tracks_points = convertPoints(gpx = tracks, siteInfo = siteInfo)
  sp::proj4string(tracks_points) = '+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
  save('tracks_points', file = 'trackPoints.Rdata')
  
} else {
  
  if(!exists("tracks_points")){
    
    load('trackPoints.Rdata')
    
  }
  
}

# Import observations. 

scatLocs = read.csv(file = 'cleanedAllScatData.csv', stringsAsFactors = F)
scatLocs = scatLocs %>% rename(Date = Day) %>% filter(Species == 'Moose')

# Convert to UTM

coordinates(scatLocs) = ~Longitude + Latitude
proj4string(scatLocs) = CRS("+init=epsg:4326")
scatLocs = spTransform(scatLocs, CRSobj = CRS(projargs = '+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))


# Trim orphaned scats ------------------------------------------------------------------------------------------------------------------------

# Should be able to match track date and name with scat observations. Check to see that this is so. Get rid of scat observations without a track for now.

# Look at 12B2 - looks good, but there are some weird collections

tracks_points %>% data.frame %>% filter(Site == '12B2') %>% pull("Date") %>% factor %>%  table
scatLocs %>% data.frame %>% filter(Site == '12B2') %>% pull(Date) %>% table

tracks_points %>% data.frame %>% filter(Site == '12B2') %>% pull("Date") %>% factor %>% table
scatLocs %>% data.frame %>% filter(Site == '12B2') %>% pull(Date) %>% table

# Okay now look at all site/date combos

dates_tracks = tracks_points %>% data.frame %>% select(Site, Date) %>% unique
dates_tracks = dates_tracks %<>% mutate(Date = as.character(as.Date(Date, format = "%Y-%m-%d")),
                                       Site = as.character(Site))

str(dates_tracks)



dates_scats = scatLocs %>% data.frame %>% select(Site, Date)
dates_scats = dates_scats %>% mutate(Date = as.character(as.Date(Date, format = "%Y-%m-%d")),
                                     Site = as.character(Site))

str(dates_scats)

dates_tracks$ID = apply(X = dates_tracks, MARGIN = 1, FUN = digest::sha1)
dates_scats$ID =  apply(X = dates_scats, MARGIN = 1, FUN = digest::sha1)

(datesLogic = !dates_scats$ID %in% dates_tracks$ID) %>% table

# Orphan scats
scatLocs[datesLogic,] %>% data.frame %>% nrow
nrow(scatLocs %>% data.frame)

# Good enough. The missing tracks are probably the incomplete ones that will need to  be stitched together. No need to do this at the moment.

scatLocs = scatLocs[!datesLogic,] # Get rid of the ones without a corresponding track.

# Track & scat overlay ------------------------------------------------------------------------------------------------------------------

# There should *always* be a track within the vicinity of a scat. See if there are any scats > 100m from a track. Use rgeos::gDistance to get the smallest distance from scat to track. 

allScatDists = vector(mode = 'numeric', length = nrow(scatLocs))


for(i in 1:nrow(scatLocs)){
  
  allScatDists[i] = scatDists(scatLocs[i,], tracks_points)
  
}

scatLocs[which(allScatDists > 100),]

# What's wrong with 6C1? 

scats_6C1 = scatLocs %>% data.frame %>% filter(Site == '06C1')
tracks_6C1 = tracks_points %>% data.frame %>% filter(Site == '06C1')


ggplot() + 
  geom_point(data = scats_6C1, aes(x = Longitude, y = Latitude, color = Date)) + 
  geom_path(data = tracks_6C1, aes(x = Easting, y = Northing, color = Date)) + 
  coord_equal()

# ONE OF THE SITES MISLABELED. EXPORT TO ARCGIS AND PERFORM SPATIAL JOIN WITH
# ORIGINAL TRANSECT LINES TO CORRECTLY IDENTIFY AND THEN RE-EXPORT ---------------------------

tracks_lines = points2line(tracks_points, ident = 'RoundBySite')

rgdal::writeOGR(obj = tracks_points, dsn = 'gpxTracksExport2017/', layer = 'gpxTracksExport2017_unclean', driver = 'ESRI Shapefile')