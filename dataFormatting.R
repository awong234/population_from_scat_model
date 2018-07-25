# Setup ------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sp)

source('functions.R')

# Import data ------------------------------------------------------------

tracks2016 = rgdal::readOGR(dsn = 'gpxTracks2016_CLEANED', layer = 'tracks2016_clean', stringsAsFactors = F)

tracks2016_points = rgdal::readOGR(dsn = 'gpxTracks2016_points_CLEANED', layer = 'tracks2016_points_clean', stringsAsFactors = F)

scats2016 = read.csv(file = 'scatLocs2016_cleaned.csv', stringsAsFactors = F)

load('scatsData_referenced.Rdata')

attr(tracks2016_points@coords, 'dimnames') = list(NULL, c("Easting", "Northing"))

# Make grid ------------------------------------------------------------------------------------------------------------------------

# Make a grid over each of the tracks

grids = makeGrid(tracks2016, data = tracks2016@data)
siteOrder = tracks2016@data$Site %>% unique %>% as.character()

# SKIP

# Look at a grid and a corresponding track

skip = T

if(!skip){
  
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
  
  rm(allIDs)
  
  # END SKIP
  
  
}


# Get records of visitations to each grid cell, by ID, and the quantity of scats picked up -----------------------------------------------------------


# Two things to do 

# FIRST - overall visitation to cells (using the track points data) -------------------------------------------------------

roundVisits = tracks2016_points@data %>% select(Site, Date, RndBySt) %>% unique %>% arrange(Site, Date)

roundVisits$VisitRank = siteVisitRank(sortedSites = roundVisits$Site)

roundVisits = roundVisits %>% arrange(Date)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
rleTracks = rlePoints(visitInfo = roundVisits, points = tracks2016_points, grids = grids, plots = F)
names(rleTracks) = roundVisits$RndBySt
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# SECOND - find cells where scats were collected (using the referenced scat location data `scatsReferenced`) --------------------------------------------

scatsReferenced_spdf = scatsReferenced %>% arrange(Time)
coordinates(scatsReferenced_spdf) = ~Easting + Northing
proj4string(scatsReferenced_spdf) = proj4string(tracks2016)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
rleScats = rlePoints(visitInfo = roundVisits, points = scatsReferenced_spdf, grids = grids, plots = F, animation = F)
names(rleScats) = roundVisits$RndBySt
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if(!skip){
  
  # Time intervals - not consistent enough to use intersection$length as a covariate. Sometimes 1s, 5s, 30s, etc.
  
  if(!skip){
    trackIntervals = tracks2016_points %>% data.frame %>% group_by(RndBySt) %>% do(out = diff(.$Time %>% as.POSIXct()))
    trackIntervals$out %>% lapply(FUN = function(x){summary(x %>% as.integer)})
  }
  
  # How many replicated observations? oof . . . 27 replicated visits out of 5
  
  rleScats %>% lapply(FUN = function(x){x$values %>% table %>% as.integer}) %>% unlist %>% table
  tracksDup = {rleTracks %>% lapply(FUN = function(x){x$values %>% table %>% as.integer}) %>% unlist %>% table %>% data.frame}
  tracksDup$Dups = tracksDup$.
  tracksDup = tracksDup %>% select(Freq, Dups) %>% mutate(Dups = Dups %>% as.integer)
  
  sum(tracksDup$Freq[tracksDup$Dups > 1])/sum(tracksDup$Freq)
  
  # Look at 266; 9A1.Sample1
  
  grid_local = grids[[which(siteOrder == '09A3')]]
  scats_local = scatsReferenced_spdf[scatsReferenced$RndBySt == '09A3.Sample1',]
  
  rle(over(x = scats_local, y = grid_local)$id)
  
  ggplot() + 
    geom_tile(data = grid_local %>% data.frame, aes(x = x, y = y), fill = NA, color = 'black') + 
    geom_text(data = grid_local %>% data.frame, aes(x = x, y = y, label = id), color = 'gray50', size = 3) + 
    geom_point(data = tracks2016_points %>% data.frame %>% filter(RndBySt == '09A3.Sample1'), aes(x = Easting, y = Northing), shape = 1) +
    geom_point(data = scats_local %>% data.frame, aes(x = Easting, y = Northing, color = 'red')) + 
    coord_equal(xlim = c(min(scats_local$Easting) + c(-100,100), max(scats_local$Easting) + c(-100,100)), 
                ylim = c(min(scats_local$Northing) + c(-100,100), max(scats_local$Northing) + c(-100,100)))

  
}

# Format as needed by model ------------------------------------------------------------------------------------------------------------------------------------------------

# What we need: 
# * y[i,t,v]
# * vis[i,t,v]
# * days[i,t-1]

# So what we're looking at, mainly, is a model of per-grid-cell initial deposition
# (lambda), per-grid-cell daily deposition (theta), probability of detection
# (p).

# There is a matrix of scats available to be sampled N. It is of dimension N[i,t,v], where 
# i = grid cell index
# t = index for sampling round (0th, 1st, 2nd, 3rd). There is no sampling performed in the 0th time period. The clearing round is, then, round 1, indexed by t = 2
# v = index for replicated observation of grid cell i
# This is unobserved, but must be initialized (?)

# There is a matrix of observations y[i,t,v] such that y[i,t,v] ~ Bin(N[i,t,v], p) , at least for the homogeneous p variety.
# y[1,t,v] is all 0's. If there are four sampling occasions max, maxT = 5. 
# Say a scat was found in a grid cell i, then another upon return to grid cell i, during time t. The structure would be y[i,t,1:2] = 1, and y[i,t,3:maxV] = 0
# I plan to start with a huge matrix of all grid cell ID's i, then cut out the indexes that don't matter (singleVisit sites, or no observations), then re-index, saving this index for reference as a data frame.

# There is a matrix of depositions R[i,t], which is applied after sampling is performed such that 
# N[i,t,1] = N[i,t-1,2] - y[i,t-1,2] + R[i,t-1]   :  Meaning, at the start of a new sampling round, what is there depends on what was there after 
# N[i,t,2] = N[i,t,1] - y[i,t,1]

# There is a matrix of visitation, indicating whether grid cells i were visited in round t, with replicate v; all 1's or 0's. If the grid cell was not visited, then p(observation) = 0 by way of multiplication: p[i,t,v] = p0[i,t,v] * vis[i,t,v]. P(observation) in t = 1 (initial deposition) is 0. 

# Site metadata  --------------------------------------------------------------------------------------------------------------------------------

# Get maxV --------------------------------------------------------
# Why am I using tracks in the next line? Well, I need to format y such that, if a grid cell was visited 20 times

maxV = rleTracks %>% lapply(FUN = function(x){x$values %>% table %>% as.integer}) %>% unlist %>% max 

# Skip

if(!skip){

    # 9A1.Clearing has a grid cell that was visited 20 times...what does this look like?
  rleTracks %>% lapply(FUN = function(x){x$values %>% table %>% as.integer}) %>% sapply(X = ., FUN = function(x){any(x == 20)}) %>% {which(.)}
  
  grid_local = grids[[which(siteOrder == '09A1')]]
  scats_local = scatsReferenced_spdf[scatsReferenced$RndBySt == '09A1.Clearing',]
  
  # It looks like that tough part in the left hand side.
  ggplot() + 
    geom_tile(data = grid_local %>% data.frame, aes(x = x, y = y), fill = NA, color = 'black') + 
    geom_text(data = grid_local %>% data.frame, aes(x = x, y = y, label = id), color = 'gray50', size = 3) + 
    geom_point(data = tracks2016_points %>% data.frame %>% filter(RndBySt == '09A1.Clearing'), aes(x = Easting, y = Northing), shape = 1) +
    geom_point(data = scats_local %>% data.frame, aes(x = Easting, y = Northing, color = 'red')) + 
    coord_equal(xlim = c(min(scats_local$Easting) + c(-100,100), max(scats_local$Easting) + c(-100,100)), 
                ylim = c(min(scats_local$Northing) + c(-100,100), max(scats_local$Northing) + c(-100,100)))
  
  # End skip
}

# Get maxI, or total number of grid cells ever visited -------------------------------------------------

# Again, why am I using rleTracks? I need to know ALL of the grid cells I ever crossed, whether we found anything or not, to build y properly.

visitedGridIDs = foreach(i = seq_along(rleTracks), .combine = c) %do% {
  
  rleTracks[[i]]$values
  
}

nGridsSampled = visitedGridIDs %>% unique %>% NROW

visitedGridCells = data.frame(gridID = visitedGridIDs %>% unique)

# Only need visited grids - that's in allGridIDs

allGridInfo = lapply(X = grids, FUN = function(x){x@data}) %>% do.call(what = rbind) %>% select(id, Site) %>% rename(gridID = id)
row.names(allGridInfo) = NULL

visitedGridInfo = visitedGridCells %>% left_join(y = allGridInfo, by = c("gridID" = "gridID"))
visitedGridInfo$y_row = 1:nrow(visitedGridInfo)

rm(visitedGridCells, allGridInfo, visitedGridIDs)


# Get maxT, or total number of rounds -----------------------------------------------------------------------

# Can use the tracks points data frame for this

maxT = {scats2016 %>% pull(RoundNo) %>% unique %>% length} + 1

if(!skip){
  
  # Intuition. The rle$values provides the order of visitation for grid cells. rle$lengths provides the quantity of scats found per grid cell. This is easily adapted to y. 
  rleScats$`12B2.Clearing` %>% do.call(what = cbind.data.frame, args = .)
  
  ID = '09A1.Clearing'
  Site = '09A1'
  grid_local = grids[[which(siteOrder == Site)]]
  scats_local = scatsReferenced_spdf[scatsReferenced$RndBySt == ID,]
  
  over(scats_local, grid_local)
  
  # It looks like that tough part in the left hand side.
  ggplot() + 
    geom_tile(data = grid_local %>% data.frame, aes(x = x, y = y), fill = NA, color = 'black') + 
    geom_text(data = grid_local %>% data.frame, aes(x = x, y = y, label = id), color = 'gray50', size = 3) + 
    geom_point(data = tracks2016_points %>% data.frame %>% filter(RndBySt == ID), aes(x = Easting, y = Northing), shape = 1) +
    geom_point(data = scats_local %>% data.frame, aes(x = Easting, y = Northing, color = 'red')) + 
    coord_equal(xlim = c(min(scats_local$Easting) + c(-100,100), max(scats_local$Easting) + c(-100,100)), 
                ylim = c(min(scats_local$Northing) + c(-100,100), max(scats_local$Northing) + c(-100,100)))
  
  
}

# y matrix -------------------------------------------------------------------------------------------

# This will be total number of scats found at a grid cell. Might end up providing covariates in this data frame later.

# Build skeleton
y = array(data = 0, dim = c(nGridsSampled, maxT, maxV))
# Fill in with data
y = make_y(gridInfo = visitedGridInfo, scatInfo = rleScats, visitInfo = roundVisits, debug = F)


if(!skip){
  
  # Check to see that it was put together correctly
  
  check12B2 = visitedGridInfo %>% filter(Site == '12B2')
  
  # Scats collected per round 
  y[check12B2$y_row, 2:5,] %>% colSums %>% rowSums
  # Is the same as . . . 
  scatsReferenced %>% filter(Site == '12B2') %>% group_by(RndBySt) %>% summarize(n = n())
  # YES!
  
  collectionsByGridRound = apply(X = y, MARGIN = c(1,2), FUN = sum)
  
  collectionsByGridRound = collectionsByGridRound[(collectionsByGridRound %>% rowSums()) > 0,]
  
  Cairo::Cairo(width = 1024*3, height = 760*3, file = 'images/collectionsPerRound.png', dpi = 150*3)
  collectionsByGridRound %>% data.frame %>% rename(Round0 = X1, Round1 = X2, Round2 = X3, Round3 = X4, Round4 = X5) %>% mutate(gridID = 1:nrow(.)) %>% reshape2:::melt.data.frame(id.vars = 'gridID', measure.vars = c(1:5)) %>% filter(variable != "Round0") %>% 
    ggplot() + 
    geom_point(aes(x = variable, y = value)) + 
    geom_line(aes(x = variable, y = value, group = gridID), size = 4, alpha = 0.01) + 
    theme_bw() + theme(axis.title.x = element_blank()) + ylab("Collections Made")
  dev.off()
  
}


# Vis matrix ------------------------------------------------------------------------------------

# Same dimension of y[,,], but it's not just 'where y == 0', because
# there are PLENTY of grid cells that were searched, perhaps multiple times,
# where no scats were found.

# maybe I can adapt make_y for this, just submitting the track points.

# This gets how many track points fall in a particular replicate, on a particular round, on a particular grid cell. As if the track points were scat points. 
vis = make_vis(gridInfo = visitedGridInfo, trackInfo = rleTracks, visitInfo = roundVisits, debug = T, siteToExamine = '12B2')

# Just turn the numbers into 1's and 0's

vis = (vis > 0) + 0 

# Double check - are all vis == 1 where y > 0? 

# TRUE
all(vis[y > 0] == 1)

# days interval matrix ------------------------------------------------------------------------------------

# We need to set N[i,1,v] to be June 1, 2016. Why? Some sites we did not get to visit until later that month! In fact, 12B1 was only first visited on July 4. 

# Therefore, having N[i,1,v] occur on a fixed date allows us to use theta[i,int] to account for the disparity in clearing date. We can assume all N[i,1,v] are Poisson random if we fix it to be the same date, but they will be different (with a left-skew in distribution) if we do not, and assume that N[i,1,v] is simply the deposition immediately before the first visit. 

# This allows more inference into theta; we can attribute *a little* bit of the variance in clearing counts due to the temporal lag between sites; all else being equal, a site sampled at a later date should have more collections on it. 

# So, what this means is that we need to append a 'pretend' round visit to roundVisits that consists of 'visiting' all the sites on June 1, and calculating the diffDays again.

depositionDays = data.frame(Site = roundVisits$Site %>% unique, Date = '2016-06-01', VisitRank = 0, stringsAsFactors = F) %>% mutate(RndBySt = paste0(Site, '.Deposition'))

roundVisits_with_time0 = tracks2016_points@data %>% select(Site, Date, RndBySt) %>% unique %>% arrange(Site, Date)

roundVisits_with_time0$VisitRank = siteVisitRank(sortedSites = roundVisits_with_time0$Site)

roundVisits_with_time0 = roundVisits_with_time0 %>% bind_rows(depositionDays) %>% arrange(Date, VisitRank)


# Get records of intervals between visits. For each site, what are the days between the visits?

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
diffDays = roundVisits_with_time0 %>% group_by(Site) %>% do(out = diff(.$Date %>% as.Date()))
names(diffDays$out) = {roundVisits_with_time0 %>% group_by(Site) %>% attr('labels')}$Site
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

select_groups(roundVisits_with_time0 %>% group_by(Site), group = 12)

# Ok, this is good, it is referenced by site. But, it needs to be referenced by each grid cell now, since the model is structured as days[i,t].
# The rows of y are referenced by visitedGridInfo.

days = make_days(diffDays = diffDays, visitedGridInfo = visitedGridInfo, maxT = maxT)

if(!skip){
  
  # Check to see that it's right. Look at 12B2. It's right.
  
  roundVisits %>% filter(Site == '12B2') # Four visits
  diffDays %>% filter(Site == '12B2') %>% pull(out)
  index = visitedGridInfo %>% filter(Site == '12B2') %>% pull(y_row)
  
  days[index,] %>% unique
  
  
}

# Replicate each row of diffDays_mat by however many grid cells are in that row's site.

data = list(y = y,
            vis = vis,
            days = days,
            nSites = nGridsSampled,
            maxT = maxT,
            maxV = maxV)

save(data, file = 'data_cleaned.Rdata')

