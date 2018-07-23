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
rgdal::writeOGR(obj = tracks_points, dsn = 'gpxTracksExport2016_unclean_points', layer = 'gpxTracksExport2016_unclean_points', driver = "ESRI Shapefile")

# Unclean shapefile was cleaned in ArcMap, but lost time information. The points data were exported to shapefile (above), and cleaned again in ArcMap.


# Load scat data, format --------------------------------------------------------------------------------------------------------------

# If clean tracks exist, bring those in.

tracks_lines = readOGR(dsn = 'gpxTracks2016_CLEANED', layer = 'tracks2016_clean')

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
  
  scatLocs = scatLocs %>% reDate(year = 2016)
  
  scatLocs = scatLocs %>% mutate(RndBySt = paste0(Site, ".", Round), Site = as.character(Site)) %>% rename(Time = Date_Time)
  
  write.csv(x = scatLocs,file = 'scatLocs2016_cleaned.csv')
  
  
}

# 2016 data, scat and transect time ------------------------------------------------------------------------------------------------------------------------------------


tracks2016_points = rgdal::readOGR(dsn = 'gpxTracks2016_points_CLEANED', layer = 'tracks2016_points_clean', stringsAsFactors = F)
attr(tracks2016_points@coords, 'dimnames') = list(NULL, c("Easting", "Northing"))

scats2016 = read.csv(file = 'scatLocs2016_cleaned.csv', stringsAsFactors = F)

tracks2016_points@data$Time = tracks2016_points@data$Time %>% as.POSIXct(format = "%Y-%m-%d %H:%M:%S")
scats2016$Time = scats2016$Time %>% as.POSIXct(format = "%Y-%m-%d %H:%M:%S")


# Okay but what about the times? Can we compare the distribution of scat and transect times? 

# Skip

timeData = data.frame(time = c(tracks2016_points@data$Time, scats2016$Time) %>% as.POSIXct(format = "%H:%M:%S"), 
                      RndBySt = c(tracks2016_points@data$RndBySt, scats2016$RndBySt), 
                      Site = c(tracks2016_points@data$Site, scats2016$Site),
                      type = c(rep('tracks', nrow(tracks2016_points)), rep('scats', nrow(scats2016)))
                      )

# Just looking at this, the scats have reasonable times - the tracks do not.
timeData %>% filter(RndBySt == '07A1.Clearing') %>% group_by(type) %>% summarize(maxT = max(time),
                                                                                 minT = min(time))
timeData %>% filter(RndBySt == '12B2.Clearing') %>% group_by(type) %>% summarize(maxT = max(time),
                                                                                  minT = min(time))

timeData %>% 
  ggplot() + 
  geom_errorbar(aes(x = RndBySt, ymin = time %>% min, ymax = time %>% max, color = RndBySt, linetype = type)) + 
  theme_bw()

num = 12
(start = -(num)+1)
(end = start + num - 1)

id = unique(timeData$RndBySt)

(start = end + 1)
(end = start + num - 1)

timeData %>% group_by(RndBySt) %>% filter(RndBySt %in% id[start:end]) %>% 
ggplot() + 
  geom_errorbar(aes(x = "Time", ymax = TimeMax, ymin = TimeMin, color = type))

# End skip

# Problem - none of the scat collection times are aligned to the track times. -----------------------------------

# Solution 1 - first scat of a transect inherits time of the closest track
# point. Offset is measured, and all future scats time values within that
# transect is offset by that amount. Scats then inherit the location of the
# closest track point *in time*.

# Solution 2 - all scats inherit the position of the closest track point.
# Probably the most practical and straightforward method, but may result in a
# few scats here or there being placed in the wrong replicate grid cell visit.
# Perhaps not a big deal. 

# DEPRECATED Solution 3 - The devices probably have a set deviation. What is this deviation? ------------------------------------------------------------------------

grouped_timeData = timeData %>% group_by(RndBySt)

groupsWithScat = attr(grouped_timeData, 'group_sizes') > 1

grouped_timeDataWithScat = grouped_timeData[(attr(grouped_timeData, "indices")[groupsWithScat] %>% unlist) + 1,]

timeDiffs = grouped_timeDataWithScat %>% summarize(timeDiffMin = diff(TimeMin),
                                                   timeDiffMax = diff(TimeMax))

timeDiffs_melt = reshape2::melt(timeDiffs)

ggplot(data = timeDiffs_melt) + 
  geom_point(aes(x = RndBySt, y = value, color = variable)) + 
  scale_color_manual(values = c('blue', 'red'), name = 'Diff Type', labels = c("Difference in earliest time at the site",
                                                                               "Difference in latest time at the site")) +
  theme_bw() + ylab('Time Difference (h)') + xlab("Site Visit") + ggtitle("Time Disparity Between Collection & Track Times") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# An idea: optimization problem, where the goal is to select offset the scat
# collection times and take the track points nearest to those offset-times in
# time, and minimize the total difference in space. Gonna be really time
# consuming, basically solving an analysis on millions of points.

# Time differences *usually* stable, but not good enough to employ over all scat/track combos. Go with solution 2; if this does not work, see if solution 1 works.

# Solution 1 - Scats inherit the spatio-temporal information of the nearest track point in space. ------------------------------------------------------------------------

# In doing so, we can tabulate from the get-go the order of visitation and
# collection of the scats - if the scats contain the same spatio-temporal
# information as the tracks, then we are essentially subsetting the tracks data
# by **only** those with scats. 

# Then, we can use rle() to get the order of visitation of the grid cells
# (rle()$values), and since there are no tracks, just scats, the rle()$lengths
# are the quantity of scats in those grid cells.

# Fields faster
microbenchmark::microbenchmark(rgeos = {rgeos::gDistance(scats2016_spdf[1,], tracks2016_points[tracks2016_points$RndBySt == '09A1.Clearing',], byid = T)},
                               fields = {fields::rdist(cbind(scats2016[1,"Easting"], scats2016[1,"Northing"]), coordinates(tracks2016_points[tracks2016_points$RndBySt == '09A1.Clearing',]))},
                               times = 500)

# Nearest track info

nearestTracks = foreach(i = 1:nrow(scats2016), .combine = rbind.data.frame) %dopar% {
  
  scatRndBySt = scats2016[i,]$RndBySt
  
  trackIndex = tracks2016_points@data$RndBySt == scatRndBySt
  
  track_local = tracks2016_points[trackIndex,]
  
  dists = fields::rdist(cbind(scats2016[i,"Easting"], scats2016[i,"Northing"]), 
                        coordinates(track_local))
  
  nearestTrackPt = track_local[which.min(dists),] %>% data.frame
  
  if(min(dists) > 100){
    
    dists_alt = fields::rdist(cbind(scats2016[i,"Easting"], scats2016[i,"Northing"]), coordinates(tracks2016_points))
    
    nearestTrackPt = tracks2016_points[which.min(dists_alt),] %>% data.frame
    
    if(min(dists_alt) < min(dists)){
      
      error = T
      
      nearestTrackPt$Dist_alt = dists_alt[which.min(dists_alt)]
      
    } 
    
  } else { 
    
    error = F
    
    nearestTrackPt$Dist_alt = NA
    
  }
  
  
  nearestTrackPt$Dist = dists[which.min(dists)]
  nearestTrackPt$ScatEasting = scats2016[i,"Easting"]
  nearestTrackPt$ScatNorthing = scats2016[i,"Northing"]
  nearestTrackPt$Error = error
  
  
  return(nearestTrackPt)
  
}

nearestTracks %>% head

nearestTracks$Dist %>% quantile(prob = seq(0,1,by = 0.01))

nearestTracks[nearestTracks$Error == T,]

scats2016[nearestTracks$Error == T,]

# Some errors found 

ggplot() + 
  geom_point(data = tracks2016_points %>% data.frame %>% filter(Site == '07A1'), aes(x = Easting, y = Northing)) + 
  # geom_text(data = tracks2016_points %>% data.frame %>% group_by(Site) %>% sample_n(1), aes(x = Easting, y = Northing, label = Site)) + 
  geom_point(data = scats2016[nearestTracks$Error,] %>% filter(Site == '03A3'), aes(x = Easting, y = Northing), color = 'blue')
