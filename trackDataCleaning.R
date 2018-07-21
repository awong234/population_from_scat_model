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
library(doParallel)

# Setup ------------------------------------------------------------------------------------------------------------------------

# Import gpx files that AREN't incomplete.
# importGPX(outPath = 'trackLogs_2017/') # default import from 2017 folder


# Initial Cleaning 2017 data ------------------------------------------------------------------------------------------------------------------------

# Summary: 

# All edits have been made to ensure clean-ness of 2017 track data. A couple of naming errors were rectified below, but the errors were identified in ArcMap, not R. 
# Further edits have been made in ArcMap, and the final clean product exported to this directory.

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
# file.rename(from = 'trackLogs_2017/08B3_8.15.17_SM.gpx', to = 'trackLogs_2017/08B4_8.15.17_SM_correct.gpx')
# file.rename(from = 'trackLogs_2017/08B4_8.15.17_SM.gpx', to = 'trackLogs_2017/08B3_8.15.17_SM_correct.gpx')
file.rename(from = 'trackLogs_2017/10B5_08.11.17_JL.gpx', to = 'trackLogs_2017/10B3_08.11.17_JL_correct.gpx')
file.rename(from = 'trackLogs_2017/12A6_07.09.17_SM.gpx', to = 'trackLogs_2017/12A4_07.09.17_SM_correct.gpx')
file.rename(from = 'trackLogs_2017/12A6_07.10.17_JL.gpx', to = 'trackLogs_2017/12A4_07.10.17_JL_correct.gpx')

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
    
    load('trackPoints_2017_unclean.Rdata')
    
  }
  
}

# Fix 12A6 site that leaks into 12A4.

# Import observations. 

scatLocs = read.csv(file = 'scats2017.csv', stringsAsFactors = F)
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
once = 1
# Run once
if(once == 1){

  any(!dates_scats$ID %in% dates_tracks$ID)
    
  datesLogic = !dates_scats$ID %in% dates_tracks$ID

  # Orphan scats
  scatLocs[datesLogic,] %>% data.frame %>% nrow
  nrow(scatLocs %>% data.frame)
  scatLocs[datesLogic,]

  scatLocs = scatLocs[!datesLogic,] # Get rid of the ones without a corresponding track.

  # Good enough. The missing tracks are probably the incomplete ones that will need to  be stitched together. No need to do this at the moment.

  once = once + 1

}

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

siteIndex = which(sites == '12A4')

tracks_lines %>% fortify %>% filter(grepl(id, pattern = sites[siteIndex] %>% as.character)) %>%
  ggplot() +
    geom_path(aes(x = long, y = lat, group = group)) +
    geom_text(aes(x = min(long), y = min(lat), label = sites[siteIndex] %>% as.character)) + 
    geom_point(data = scatLocs %>% data.frame %>% filter(Site == sites[siteIndex]), aes(x = Longitude, y = Latitude), color = 'red') +
    coord_equal()

siteIndex = siteIndex + 1

# In depth look

sitesByRound = tracks_lines %>% fortify %>% filter(grepl(id, pattern = '12A4')) %>% pull(id) %>% unique

siteIndex = 0

siteIndex = siteIndex + 1

tracks_lines %>% fortify %>% filter(grepl(id, pattern = sitesByRound[siteIndex] %>% as.character)) %>% 
  ggplot() +
  geom_path(aes(x = long, y = lat, group = group)) +
  geom_text(aes(x = min(long), y = min(lat), label = sitesByRound[siteIndex] %>% as.character)) + 
  coord_equal()

tracks_lines@data %>% filter(id == sitesByRound[siteIndex])

list.files(path = 'trackLogs_2017/', pattern = '12A4|12A6')

# Export to shapefile 

rgdal::writeOGR(obj = tracks_lines, dsn = 'gpxTracksExport2017_unclean', layer = 'gpxTracksExport2017_unclean', driver = 'ESRI Shapefile')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# Initial cleaning 2016 Data ------------------------------------------------------------------------------------------------------------------------

# importGPX(startPath = 'D:/MooseArchive/Field Data 2016 (Archive)/Moose 2016/Dog Tracklogs/', outPath = 'trackLogs_2016/') # import from 2016 folder

# Some naming errors. Correct here.

file.rename(from = 'trackLogs_2016/11B2_06_01_2016_JB_incomplete.gpx', to = 'trackLogs_2016/11B2_06.01.16_JB_incomplete.gpx')
file.rename(from = 'trackLogs_2016/3A1_08,24.16_SM.gpx', to = 'trackLogs_2016/3A1_08.24.16_SM.gpx')
file.rename(from = 'trackLogs_2016/10A3_8.14.16_SM.gpx', to = 'trackLogs_2016/10A3_08.14.16_SM.gpx')
file.rename(from = 'trackLogs_2016/10A2_8.14.16_SM.gpx', to = 'trackLogs_2016/10A2_08.14.16_SM.gpx')

# Find more naming errors

fileNames = list.files(path = 'trackLogs_2016/', pattern = '.gpx')

siteDates = fileNames %>% {regmatches(x = ., m = regexec(pattern = "\\d+\\.\\d+\\.\\d+", text = ., perl = T))} %>% as.character()

# Sites should be formatted like mm.dd.yy

correctForm = grepl(pattern = '^\\d{2}\\.\\d{2}\\.\\d{2}$', x = siteDates)

siteDates[correctForm] = as.Date(siteDates[correctForm], format = '%m.%d.%y')
siteDates[!correctForm] = as.Date(siteDates[!correctForm], format = '%m.%d.%Y')

siteDates = siteDates %>% as.integer %>% as.Date(origin = '1970-01-01')

siteInfo = siteInfoFromFileName(path = 'trackLogs_2016/', optDates = siteDates) %>% mutate(siteID = as.character(siteID))

cbind.data.frame(siteInfo, fileNames) %>% sample_n(size = 5)

# Rename site names to conform to 2017 style 'ddLd'.

index = grepl(pattern = '^\\d{1}\\w\\d', x = siteInfo$siteID)

siteInfo$siteID[index] = paste0('0', siteInfo$siteID[index])

if(!"trackPoints_2016_unclean.Rdata" %in% dir()){
  
  tracks = getGPX(path = 'trackLogs_2016/', siteInfo = siteInfo, debug = F)
  tracks_points = convertPoints(gpx = tracks, siteInfo = siteInfo, survey_year = 2016)
  sp::proj4string(tracks_points) = '+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' # May not be necessary if done in-function.
  save('tracks_points', file = 'trackPoints_2016_unclean.Rdata')
  
} else {
  
  if(!exists("tracks_points")){
    
    load('trackPoints_2016_unclean.Rdata')
    
  }
  
}

# Export over to shapefile

tracks_lines = points2line(tracks_points, ident = "RoundBySite")

proj4string(tracks_lines) = proj4string(tracks_points)

sites = tracks_points@data$Site %>% unique

rgdal::writeOGR(obj = tracks_lines, dsn = 'gpxTracksExport2016_unclean', layer = 'gpxTracksExport2016_unclean', driver = 'ESRI Shapefile')

# Load scat data, format --------------------------------------------------------------------------------------------------------------

# If clean tracks exist, bring those in.

tracks2016_clean = readOGR(dsn = 'gpxTracks2016_CLEANED', layer = 'gpxTracks2016_CLEANED')

tracks2016_clean@lines[[1]]@Lines[[1]]@coords %>% cbind.data.frame(., 1:nrow(.)) %>% rename(Easting = `1`, Northing = `2`, ID = `1:nrow(.)`) %>% 
  ggplot() + geom_point(aes(x = Easting, y = Northing, color = ID))

if(file.exists('scatLocs2016_cleaned.csv')){
  
  scatLocs = read.csv('scatLocs2016_cleaned.csv', stringsAsFactors = F)
  
} else {
  
  scatLocs = read.csv(file = 'scats2016.csv', stringsAsFactors = F)
  scatLocs$Date = as.Date(scatLocs$Date_Time)
  
  # Rename sites to conform to 2017 structure.
  
  index = grepl(pattern = '^\\d{1}\\w\\d', x = scatLocs$Site)
  scatLocs$Site[index] = paste0('0', scatLocs$Site[index])
  
  # Inspection 
  
  # Look at both overlaid
  
  skip = T
  
  if(!skip){
    
    siteIndex = 1
    
    siteIndex = siteIndex + 1
    
    tracks_lines %>% fortify %>% filter(grepl(id, pattern = sites[siteIndex] %>% as.character)) %>%
      ggplot() +
      geom_path(aes(x = long, y = lat, group = group)) +
      geom_text(aes(x = min(long), y = min(lat), label = sites[siteIndex] %>% as.character)) + 
      geom_point(data = scatLocs %>% data.frame %>% filter(Site == sites[siteIndex]), aes(x = Easting, y = Northing), color = 'red') +
      coord_equal()
    
    
  }
  
  
  # Find tsect nearest to scats to validate site ID
  
  scatLocs_spdf = scatLocs
  coordinates(scatLocs_spdf) = ~Easting + Northing
  
  scatTrackCompare = foreach(i = 1:nrow(scatLocs), .combine = rbind.data.frame) %do% {
    
    siteScat = scatLocs_spdf[i,"Site"]
    
    allDist = rgeos::gDistance(scatLocs_spdf[i,], tracks_lines, byid = T)
    minDistTrack = tracks_lines@data[which.min(allDist),"Site"]
    
    out = data.frame(scatSite = siteScat, nearestSite = minDistTrack)
    
    return(out)
    
  }
  
  compare = scatTrackCompare %>% rowwise %>% summarize(compare = scatSite.Site == nearestSite) 
  compare = compare$compare
  
  badScats = scatTrackCompare[!compare,]
  
  scatLocs[!compare,]$Site = scatTrackCompare[!compare, "nearestSite"] %>% as.character()
  
  # any orphan scats?
  
  tracks_points = tracks2016_clean %>% as("SpatialPointsDataFrame")
  
  dates_tracks = tracks_points %>% data.frame %>% select(Site, Date) %>% unique
  dates_tracks = dates_tracks %>% mutate(Date = as.character(as.Date(Date, format = "%Y-%m-%d")),
                                         Site = as.character(Site))
  
  str(dates_tracks)
  
  dates_scats = scatLocs %>% data.frame %>% select(Site, Date_Time)
  dates_scats = dates_scats %>% mutate(Date = as.character(as.Date(Date_Time, format = "%Y-%m-%d")),
                                       Site = as.character(Site)) %>% select(Site, Date)
  
  str(dates_scats)
  
  dates_tracks$ID = apply(X = dates_tracks, MARGIN = 1, FUN = digest::sha1)
  dates_scats$ID =  apply(X = dates_scats, MARGIN = 1, FUN = digest::sha1)
  
  any(!dates_scats$ID %in% dates_tracks$ID)
  
  datesLogic = !dates_scats$ID %in% dates_tracks$ID
  
  # Orphan scats
  scatLocs[datesLogic,] %>% data.frame %>% nrow
  nrow(scatLocs %>% data.frame)
  scatLocs[datesLogic,]
  
  # Remove orphan scats ; no use without a track
  
  scatLocs = scatLocs[!datesLogic,]
  
  write.csv(x = scatLocs,file = 'scatLocs2016_cleaned.csv')
  
  
}

