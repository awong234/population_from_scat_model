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
# importGPX(outPath = 'trackLogs_2017/') # default import from 2017 folder
# importGPX(startPath = 'D:/MooseArchive/Field Data 2016 (Archive)/Moose 2016/Dog Tracklogs/', outPath = 'trackLogs_2016/') # import from 2016 folder



# Initial Cleaning 2017 data ------------------------------------------------------------------------------------------------------------------------

# There were some errors in naming, found in commented section below. Fix these HERE, explicitly.

# MISLABELS

# 06C1-2017-08-27 should be 06C2. 
# 06C2-2017-08-27 should be 06C1.
# 08B3-2017-08-15 should be 08B4.
# 08B4-2017-08-15 should be 08B3.
# 10B5-2017-08-11 should be 10B3.
# 12A6-2017-07-09 should be 12A4.

file.rename(from = 'trackLogs_2017/06C1_08.27.17_JL.gpx', to = 'trackLogs_2017/06C2_08.27.17_JL_correct.gpx')
file.rename(from = 'trackLogs_2017/06C2_08.27.17_JL.gpx', to = 'trackLogs_2017/06C1_08.27.17_JL_correct.gpx')
file.rename(from = 'trackLogs_2017/08B3_8.15.17_SM.gpx', to = 'trackLogs_2017/08B4_8.15.17_SM_correct.gpx')
file.rename(from = 'trackLogs_2017/08B4_8.15.17_SM.gpx', to = 'trackLogs_2017/08B3_8.15.17_SM_correct.gpx')
file.rename(from = 'trackLogs_2017/10B5_08.11.17_JL.gpx', to = 'trackLogs_2017/10B3_08.11.17_JL_correct.gpx')
file.rename(from = 'trackLogs_2017/12A6_07.09.17_SM.gpx', to = 'trackLogs_2017/12A4_07.09.17_SM_correct.gpx')

# A bunch of Jake's gpx files were not separated out. Those are the weird ones; fixed in ArcMap, now just need to load the shapefiles. Delete from set

suppressWarnings(file.remove('trackLogs_2017/07.16-20.17_JL.gpx'))
suppressWarnings(file.remove('trackLogs_2017/07.22-25.17_JL.gpx'))

# Metadata

siteInfo = siteInfoFromFileName(path = 'trackLogs_2017/')

# Only reload all tracks if the .Rdata file somehow gets lost.

if(!"trackPoints_2017.Rdata" %in% dir()){
  
  tracks = getGPX(path = 'trackLogs_2017/', siteInfo = siteInfo, debug = F)
  tracks_points = convertPoints(gpx = tracks, siteInfo = siteInfo)
  sp::proj4string(tracks_points) = '+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' # May not be necessary if done in-function.
  save('tracks_points', file = 'trackPoints_2017.Rdata')
  
} else {
  
  if(!exists("tracks_points")){
    
    load('trackPoints_2017.Rdata')
    
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

skip = T

if(!skip){
  
  # Look at 12B2 - looks good, but there are some weird collections
  
  tracks_points %>% data.frame %>% filter(Site == '12B2') %>% pull("Date") %>% factor %>%  table
  scatLocs %>% data.frame %>% filter(Site == '12B2') %>% pull(Date) %>% table
  
  tracks_points %>% data.frame %>% filter(Site == '12B2') %>% pull("Date") %>% factor %>% table
  scatLocs %>% data.frame %>% filter(Site == '12B2') %>% pull(Date) %>% table
  
}


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

# Are there any scats with site/date ID

# Run once

any(!dates_scats$ID %in% dates_tracks$ID)
  
datesLogic = !dates_scats$ID %in% dates_tracks$ID

# Orphan scats
scatLocs[datesLogic,] %>% data.frame %>% nrow
nrow(scatLocs %>% data.frame)
scatLocs[datesLogic,]

scatLocs = scatLocs[!datesLogic,] # Get rid of the ones without a corresponding track.

# Good enough. The missing tracks are probably the incomplete ones that will need to  be stitched together. No need to do this at the moment.



# Track & scat overlay ------------------------------------------------------------------------------------------------------------------

# There should *always* be a track within the vicinity of a scat. See if there are any scats > 100m from a track. Use rgeos::gDistance to get the smallest distance from scat to track.

# NOT RUN

if(!skip){
  
  allScatDists = vector(mode = 'numeric', length = nrow(scatLocs))
  
  
  for(i in 1:nrow(scatLocs)){
    
    allScatDists[i] = scatDists(scatLocs[i,], tracks_points)
    
  }
  
  scatLocs[which(allScatDists > 100),]
  
  # What's wrong with 6C1? 
  
  scats_6C1 = scatLocs %>% data.frame %>% filter(Site == '06C1')
  tracks_6C1 = tracks_points %>% data.frame %>% filter(Site == '06C1')
  
  
  ggplot() + 
    geom_point(data = scats_6C1, aes(x = Longitude, y = Latitude, color = Date %>% factor)) + 
    geom_path(data = tracks_6C1, aes(x = Easting, y = Northing, color = Date %>% factor)) + 
    coord_equal()
  
}

# # Deprecated - problems fixed above
# 
# # ONE OF THE SITES MISLABELED. EXPORT TO ARCGIS AND PERFORM SPATIAL JOIN WITH
# # ORIGINAL TRANSECT LINES TO CORRECTLY IDENTIFY AND THEN RE-EXPORT ---------------------------
# 
tracks_lines = points2line(tracks_points, ident = 'RoundBySite')

proj4string(tracks_lines) = proj4string(tracks_points)

sites = tracks_points@data$Site %>% unique

siteIndex = 1

tracks_lines %>% fortify %>% filter(grepl(id, pattern = sites[siteIndex] %>% as.character)) %>%
  ggplot() +
    geom_path(aes(x = long, y = lat, group = group)) +
    geom_text(aes(x = min(long), y = min(lat), label = sites[siteIndex] %>% as.character))

siteIndex = siteIndex + 1

rgdal::writeOGR(obj = tracks_lines, dsn = 'gpxTracksExport2017_unclean', layer = 'gpxTracksExport2017_unclean', driver = 'ESRI Shapefile')



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# Initial cleaning 2016 Data ------------------------------------------------------------------------------------------------------------------------

