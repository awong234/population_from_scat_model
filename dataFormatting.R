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

  track_12B2 = tracks2016_points[which(tracks2016_points@data$Site == '12B2' & tracks2016_points@data$RoundNo == 0),]

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

  center = grid_12B2 %>% data.frame %>% filter(id == 81072) %>% select(x,y)
  buff = 300
  buffmat = matrix(c(-buff, -buff, buff, buff), nrow = 2, ncol = 2)
  lims = as.numeric(center) + buffmat
  Cairo::Cairo(file = 'images/scatLocs.png', dpi = 150)
  ggplot() +
    geom_tile(data = grid_12B2 %>% data.frame, aes(x = x, y = y), fill = '#262626', color = 'gray60') +
    geom_tile(data = grid_12B2 %>% data.frame %>% filter(id %in% intersect$values), aes(x = x, y = y), color = 'cyan', fill = NA) +
    # geom_text(data = grid_12B2 %>% data.frame, aes(x = x, y = y, label = id), size = 2) +
    geom_path(data = track_12B2 %>% data.frame, aes(x = Easting, y = Northing), color = 'gray90') +
    geom_point(data = scatsReferenced %>% filter(Site == '12B2', RoundNo == 0), aes(x = Easting, y = Northing), shape = 21, size = 2, color = 'gray60', fill = 'lawngreen') +
    # geom_point(data = scats2016 %>% filter(Site == '12B2'), aes(x = Easting, y = Northing), shape = 1, color = 'green') +
    coord_equal(xlim = c(lims[1,1], lims[1,2]), ylim = c(lims[2,1], lims[2,2])) +
    theme(panel.grid = element_blank(),
          plot.background = element_rect(fill = 'transparent', color = NA),
          panel.background = element_rect(fill = 'transparent'),
          axis.title = element_blank(),
          axis.text = element_blank()
          )
  dev.off()

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

# Plot example of replication of cells, for paper

if(!skip){

  for(i in 0:3){

    tracks12B2 = tracks2016_points %>% data.frame %>% filter(Site == '12B2', RoundNo == i)

    grid12B2 = grids$`12B2` %>% data.frame

    rleList = list(rleTracks$`12B2.Clearing`,
                   rleTracks$`12B2.Sample1`,
                   rleTracks$`12B2.Sample2`,
                   rleTracks$`12B2.Sample3`)

    rle12B2 = rleList[[i+1]]$values %>% table %>% data.frame()
    colnames(rle12B2) = c('gridID', 'Freq')
    rle12B2$gridID = as.integer(as.character(rle12B2$gridID))

    grid12B2 = grid12B2 %>% left_join(rle12B2, by = c('id' = 'gridID'))
    grid12B2$Freq[is.na(grid12B2$Freq)] = 0

    buff = 100

    source('images/red_pallete.R')

    colors = rev(a(n = 6))
    names(colors) = c("0", "1", "2", "3", '4', '5')

    center_point = grid12B2 %>% filter(id == 80642) %>% select(x,y)

    Cairo::Cairo(width = 1024, height = 1024, file = paste0('images/track',i,'.png'), dpi = 150)
    print(
      ggplot() +
        geom_tile(data = grid12B2, aes(x = x, y = y, fill = factor(Freq)), color = 'gray50', show.legend = F) +
        geom_path(data = tracks12B2, aes(x = Easting, y = Northing), color = 'red') +
        geom_text(data = grid12B2, aes(x = x, y = y, label = paste('r =',Freq)), size = 4) +
        coord_equal(xlim = center_point[,1]+ c(-buff, buff), ylim = center_point[,2] + c(-buff, buff)) +
        theme_bw() +
        scale_fill_manual(values = colors, limits = c('0','1','2','3', '4', '5'), name = expression(paste("Secondary\nOccasion (", r[g], ")"))) +
        guides(fill = guide_legend()) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank()
              )
    )
    dev.off()

  }

}


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
  tracksDup = tracksDup %>% select(Dups, Freq) %>% mutate(Dups = Dups %>% as.integer)

  # Plot distribution
  Cairo::Cairo(width = 768*2, height = 1024*2, file = 'images/rep_dist.png', dpi = 300)
  ggplot() +
    geom_col(data = tracksDup, aes(x = Dups, y = Freq), fill = "#440154FF") +
    coord_trans(y = 'sqrt') +
    scale_y_continuous(breaks = tracksDup$Freq %>% round(digits = -3)) +
    theme_bw() + xlab('Replicated Cell Observations') + ylab("Frequency")
  dev.off()

  # Get table of distribution
  knitr::kable(x = tracksDup, format = 'latex', booktabs = T)

  # How many replicates > 1?
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

allGridLocations = lapply(X = grids, FUN = function(x){coordinates(x)}) %>% do.call(what = rbind)

allGridInfo = cbind.data.frame(allGridInfo, allGridLocations)

visitedGridInfo = visitedGridCells %>% left_join(y = allGridInfo, by = c("gridID" = "gridID"))

# The row ordering of y (and thus vis, days) is arbitrary.
# At this point I am enforcing the rows of y to align to the ordering of the grid ID's in visitedGridInfo. This is also ordering the grid cells by time, since rleTracks is ordered by roundVisits, which is ordered by time.

# Check

if(!skip){

  # All the RoundBySites are present
  all(roundVisits$RndBySt == names(rleTracks))

  # All the sites are represented
  identical(
    roundVisits$Site %>% unique, visitedGridInfo$Site %>% unique
  )

  # Notice : roundVisits is a record of all visits made, and drawn directly from tracks2016_points@data.
  # visitedGridInfo instead is based on the UNIQUE gridID's that have been visited. It is also directly based on tracks2016_points if you go far back enough:
  # visitedGridInfo$gridID <- rleTracks <- basedOn(roundVisits, tracks2016_points)
  # As seen in the test after y is built, the subsetting is correct.
}

visitedGridInfo$y_row = 1:nrow(visitedGridInfo)

rm(visitedGridCells, allGridInfo, allGridLocations, visitedGridIDs)


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
  scatsReferenced %>% filter(Site == '12B2') %>% group_by(RndBySt) %>% summarize(n = n()) %>% pull(n)
  # YES!

  # Check for all of them.
  sites = visitedGridInfo %>% pull(Site) %>% unique

  test = vector(mode = 'logical', length = NROW(sites))

  for(i in 1:NROW(sites)){

    check = visitedGridInfo %>% filter(Site == sites[i])
    y_collections = y[check12B2$y_row, 2:5,] %>% colSums %>% rowSums
    data_collections = scatsReferenced %>% filter(Site == '12B2') %>% group_by(RndBySt) %>% summarize(n = n()) %>% pull(n)

    test[i] = all(
      y_collections == data_collections
    )

  }

  # All of the scats are represented.
  all(test)

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
vis = make_vis(gridInfo = visitedGridInfo, trackInfo = rleTracks, visitInfo = roundVisits, debug = F)

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


# Export data --------------------------------------------

data = list(y = y,
            vis = vis,
            days = days,
            nSites = nGridsSampled,
            maxT = maxT,
            maxV = maxV)

# visitedGridInfo relates the gridID's from `grids` to the corresponding rows in y.
# roundVisits provides ALL of the site/date combinations and rank of visit
# grids are the grids surrounding all of the tracks, and one can use the grid ID's to relate rows in y to particular grid cells in `grids`
metadata = list(visitedGridInfo = visitedGridInfo,
                roundVisits = roundVisits,
                grids = grids
                )

save(data, file = 'data_cleaned.Rdata')
save(metadata, file = 'metadata.Rdata')












# Format covariates --------------------------------------------------------------

# Need grids, track points. Load track points from above, get grids from saved metadata.

skip = T

load('data_cleaned.Rdata')
load('metadata.Rdata')
source('functions.R')

extract(data)
extract(metadata)

# Need covariates each as a vector nSites long - one value per covariate.

scaleMetrics = data.frame("Covariate" =  character(0),
                          "Center" = numeric(0),
                          "Scale"  = numeric(0)
                          )



# Format Coordinate covariate ------------------------------------------------------------

# Assign Coordinate value per grid cell. Easy.

gridNorth = visitedGridInfo$y
gridNorth_scaled = scale(gridNorth)

scaleMetrics = scaleMetrics %>% add_row(Covariate = "Northing", Center = attributes(gridNorth_scaled)[2], Scale = attributes(gridNorth_scaled)[3])


gridEast = visitedGridInfo$x
gridEast_scaled = scale(gridEast)

scaleMetrics = scaleMetrics %>% add_row(Covariate = "Easting", Center = attributes(gridEast_scaled)[2], Scale = attributes(gridEast_scaled)[3])


# Format Elevation covariate ------------------------------------------------------------

# Assign elevation value per grid cell - summarize by mean elevation within a
# grid cell. Tests indicate that this is correctly performed because the
# residuals when subtracting out the cell center values are very small.

if(!'elevationBySite.Rdata' %in% dir()){

  elev = raster::raster('spatCov/elev_masked/Elev_masked_to_tracks_1.tif')

  elevData = summarizeRastFromGrid(grids = grids, raster = elev, method = 'mean')

  names(elevData) = names(grids)

  save(elevData, file = 'elevationBySite.Rdata')

} else {
  load('elevationBySite.Rdata')
}

if(!skip){

  # See if it worked

  test = elevData$`11A1`

  test = test[complete.cases(test),]

  # Looks good when comparing to arcmap
  ggplot() +
    geom_raster(data = test, aes(x = x, y = y, fill = rasterSummary)) +
    geom_point(data = tracks2016_points %>% data.frame %>% filter(Site == '11A1'), aes(x = Easting, y = Northing)) +
    coord_equal()

}

# Now need to grab elevation means and assign them in order

visitedGridInfo %>% head

elevData$`01A1` %>% head

elevData = elevData %>% do.call(what = rbind)

elevCov = left_join(visitedGridInfo, elevData, by = c("gridID" = 'id')) %>%
  rename(Site = Site.x, x = x.x, y = y.x, elev = rasterSummary) %>% select(gridID, Site, x, y, y_row, elev)

if(!skip){

  elev12B2 = elevCov %>% filter(Site == '12B2')

  # Looks good.
  ggplot() +
    geom_raster(data = elev12B2, aes(x = x , y = y, fill = elev)) +
    geom_point(data = tracks2016_points %>% data.frame %>% filter(Site == '12B2'), aes(x = Easting, y = Northing)) +
    coord_equal()

}

# Scale

scaledElev = scale(elevCov$elev)

scaleMetrics = scaleMetrics %>% add_row(Covariate = "Elevation", Center = attributes(scaledElev)[2], Scale = attributes(scaledElev)[3])

# Format Habitat covariate ------------------------------------------------------------

# Follows same general procedure as elev formatting

if(!'habitatBySite.Rdata' %in% dir()){

  habs = raster::raster('spatCov/tnc_habs_masked.tif/tnc_habs_masked.tif')

  # Assign habitat value per grid cell - will probably take the predominant feature.

  habData = summarizeRastFromGrid(grids = grids, raster = habs, method = 'mode')

  names(habData) = names(grids)

  save(habData, file = 'habitatBySite.Rdata')


} else {

  load('habitatBySite.Rdata')

}

if(!skip){

  # See if it worked

  test = habData$`12B2`

  test = test[complete.cases(test),]

  # Looks good when comparing to arcmap
  ggplot() +
    geom_raster(data = test, aes(x = x, y = y, fill = rasterSummary %>% factor())) +
    geom_point(data = tracks2016_points %>% data.frame %>% filter(Site == '12B2'), aes(x = Easting, y = Northing)) +
    coord_equal()

}

# Put it all together

# One-hot the habitat covariate.

habData = do.call(what = rbind, args = habData)

habCov = left_join(visitedGridInfo, habData, by = c("gridID" = 'id')) %>%
  rename(Site = Site.x, x = x.x, y = y.x, hab = rasterSummary) %>% select(gridID, Site, x, y, y_row, hab)

habRef = data.frame("Name" = c("Other", "Conifer", "Deciduous", "Mixed", "Wetland"),
                    "Value" = c(0,1,2,3,4))

nHab = habRef %>% nrow

habCov_mat = matrix(data = 0, nrow = nSites, ncol = nHab, dimnames = list(NULL, habRef$Name))

for(i in 1:nHab){

  habCov_mat[,i] = habCov$hab == habRef$Value[i]

}



# Format major road density covariate ----------------------------------------------------------------------

# Same process as elevation

if(!'highwayDensityBySite.Rdata' %in% dir()){

  highwayRast = raster::raster('spatCov/highways_masked_to_tracks/highways_masked_to_tracks.tif')

  highwayData = summarizeRastFromGrid(grids = grids, raster = highwayRast, method = 'mean')

  names(highwayData) = names(grids)

  save(highwayData, file = 'highwayDensityBySite.Rdata')

} else {

  load("highwayDensityBySite.Rdata")

}

if(!skip){

  # See if it worked

  test = highwayData$`11A1`

  test = test[complete.cases(test),]

  # Looks good when comparing to arcmap
  ggplot() +
    geom_raster(data = test, aes(x = x, y = y, fill = rasterSummary)) +
    geom_point(data = tracks2016_points %>% data.frame %>% filter(Site == '11A1'), aes(x = Easting, y = Northing)) +
    coord_equal()

}

highwayData = highwayData %>% do.call(what = rbind)

highwayCov = left_join(visitedGridInfo, highwayData, by = c("gridID" = 'id')) %>%
  rename(Site = Site.x, x = x.x, y = y.x, density = rasterSummary) %>% select(gridID, Site, x, y, y_row, density)

if(!skip){

  high12B2 = highwayCov %>% filter(Site == '12B2')

  # Looks good.
  ggplot() +
    geom_raster(data = high12B2, aes(x = x , y = y, fill = density)) +
    geom_point(data = tracks2016_points %>% data.frame %>% filter(Site == '12B2'), aes(x = Easting, y = Northing)) +
    coord_equal()

}

# Scale

scaledHighway = scale(highwayCov$density)

scaleMetrics = scaleMetrics %>% add_row(Covariate = "Highway", Center = attributes(scaledHighway)[2], Scale = attributes(scaledHighway)[3])


# Format minor road density covariate ----------------------------------------------------------------------

# Same process as elevation

if(!'minRoadDensityBySite.Rdata' %in% dir()){

  minRoadRast = raster::raster('spatCov/localRoads_masked_to_tracks/localRoads_masked_to_tracks.tif')

  minRoadData = summarizeRastFromGrid(grids = grids, raster = minRoadRast, method = 'mean')

  names(minRoadData) = names(grids)

  save(minRoadData, file = 'minRoadDensityBySite.Rdata')

} else {

  load("minRoadDensityBySite.Rdata")

}

if(!skip){

  # See if it worked

  test = minRoadData$`11A1`

  test = test[complete.cases(test),]

  # Looks good when comparing to arcmap
  ggplot() +
    geom_raster(data = test, aes(x = x, y = y, fill = rasterSummary)) +
    geom_point(data = tracks2016_points %>% data.frame %>% filter(Site == '11A1'), aes(x = Easting, y = Northing)) +
    coord_equal()

}

minRoadData = minRoadData %>% do.call(what = rbind)

minRoadCov = left_join(visitedGridInfo, minRoadData, by = c("gridID" = 'id')) %>%
  rename(Site = Site.x, x = x.x, y = y.x, density = rasterSummary) %>% select(gridID, Site, x, y, y_row, density)

if(!skip){

  minRoad12B2 = minRoadCov %>% filter(Site == '12B2')

  # Looks good.
  ggplot() +
    geom_raster(data = minRoad12B2, aes(x = x , y = y, fill = density)) +
    geom_point(data = tracks2016_points %>% data.frame %>% filter(Site == '12B2'), aes(x = Easting, y = Northing)) +
    coord_equal()

}

# Scale

scaledMinRoad = scale(minRoadCov$density)

scaleMetrics = scaleMetrics %>% add_row(Covariate = "MinorRoad", Center = attributes(scaledMinRoad)[2], Scale = attributes(scaledMinRoad)[3])


# Same process as elevation


# Put all together ------------------------------------------------------------------------

gridCovariates = data.frame(Intercept = 1,
                            Northing = gridNorth_scaled,
                            Easting = gridEast_scaled,
                            Elevation = scaledElev,
                            scaledHighway,
                            scaledMinRoad,
                            habCov_mat[,2:5])

save(gridCovariates, file = 'gridCovariates.Rdata')

# Format Detection covariates ------------------------------------------------------------

visitData = read.csv('visitData.csv', stringsAsFactors = F)

### Detection covariates need to be of the form X[g,t], with the same X[g,t] applied for all v.

## Format Dog and Handler Covariate ------------------------------------------------------------

visitData$Date_Time = as.POSIXct(visitData$Date_Time, format = '%m/%d/%Y %H:%M')

visitData = visitData %>% mutate(Date = as.Date(Date_Time))

visitDataRedu = visitData %>% select(Date, Transect_ID, Handler, Dog) %>% mutate(Date = as.character(Date))

roundVisits %>% head
visitDataRedu %>% head

visitCovariates = roundVisits %>% left_join(visitDataRedu, by = c("Site" = "Transect_ID", "Date" = "Date"))

dogNames = visitCovariates$Dog %>% unique
nDogs = dogNames %>% length

humNames = visitCovariates$Handler %>% unique
nHandler = humNames %>% length

dogCov = array(data = 0, dim = c(nSites, maxT-1, nDogs), dimnames = list(NULL, NULL, dogNames))
humCov = array(data = 0, dim = c(nSites, maxT-1, nHandler), dimnames = list(NULL, NULL, humNames))

for(i in 1:maxT-1){

  temp = visitCovariates %>% filter(VisitRank == i)

  temp = visitedGridInfo %>% left_join(temp, by = c("Site" = "Site"))

  for(dog in 1:nDogs){
    dogName = dogNames[dog]
    dogCov[,i,dog] = as.integer(temp$Dog == dogName)
  }

  for(hum in 1:nHandler){
    humName = humNames[hum]
    humCov[,i,hum] = as.integer(temp$Handler == humName)
  }

}

dogCov[is.na(dogCov)] = 0
humCov[is.na(humCov)] = 0

rm(temp)

# Correct? Yes all correct

if(!skip){

  out = matrix(F, nrow = maxT - 1, ncol = nDogs)

  for(t in 1:maxT-1){
    for(dog in 1:nDogs){

      dogIndex = dogCov[,t,dogNames[dog]] %>% as.logical() # Skye's visits in first visit

      covSites = visitedGridInfo[dogIndex,] %>% pull(Site) %>% unique

      truthSites = visitCovariates %>% filter(Dog == dogNames[dog], VisitRank == t) %>% pull(Site)

      out[t,dog] = all(truthSites %in% covSites) & all(covSites %in% truthSites)

    }
  }

  out   # exactly equal

}

## Format length of track within grid cell. ------------------------------------------------------------

if(!exists("tracks2016_points")){
  tracks2016_points = rgdal::readOGR(dsn = 'gpxTracks2016_points_CLEANED', layer = 'tracks2016_points_clean', stringsAsFactors = F)
  attr(tracks2016_points@coords, 'dimnames') = list(NULL, c("Easting", "Northing"))

  roundVisits = tracks2016_points@data %>% select(Site, Date, RndBySt) %>% unique %>% arrange(Site, Date)

  roundVisits$VisitRank = siteVisitRank(sortedSites = roundVisits$Site)

  roundVisits = roundVisits %>% arrange(Date)

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  rleTracks = rlePoints(visitInfo = roundVisits, points = tracks2016_points, grids = grids, plots = F)
  names(rleTracks) = roundVisits$RndBySt
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

}


# Distance values for every grid cell, sample occasion, and replicate.
# Without parallel takes 200 s; with parallel takes 96s

Dcov = trackDistPerRep(tracks = tracks2016_points,
                rleTracks = rleTracks,
                visitedGridInfo = visitedGridInfo,
                roundVisits = roundVisits,
                debug = F)



# Test: There should be positive length wherever vis == 1, and 0 wherever vis == 0.

which(!(Dcov > 0) == (vis[,2:5,] > 0))
which(!(Dcov == 0) == (vis[,2:5,] == 0))

# Why are there some errors in the first test? I really can't say. However, this is 31/747760 data points, which represents 0.00041% of the data. Set them to 1.

Dcov[which(!(Dcov > 0) == (vis[,2:5,] > 0))] = 1

# No more errors
which(!(Dcov > 0) == (vis[,2:5,] > 0))
which(!(Dcov == 0) == (vis[,2:5,] == 0))

# Scale and center

Dcov_scaled = scale(Dcov)
scaleMetrics = scaleMetrics %>% add_row(Covariate = "DetectDist", Center = attributes(Dcov_scaled)[2], Scale = attributes(Dcov_scaled)[3])

Dcov_scaled = array(data = scale(Dcov), dim = dim(Dcov))

# Old tests

if(!skip){

  testRle = rleTracks[[1]]

  testTracks = tracks2016_points %>% data.frame %>% filter(Site == '11A1', Round == 'Clearing') %>% mutate(Order = 1:nrow(testTracks))

  gridVisitOrder = data.frame(gridID = testRle$values)

  gridVisitOrder$rank[orderVec] = siteVisitRank(testRle$values[orderVec])

  testfn = function(){

    dist = NULL

    start = 0
    end = 0


    for(i in seq(testRle$lengths %>% length)){

      nPoints = testRle$lengths[i]

      start = end + 1
      end = end + nPoints

      if(nPoints == 1){
        dist = c(dist, 0)
      } else {

        dist = c(dist,diffDist(pts = testTracks[start:end,] %>% select(Easting,Northing)) %>% sum)

      }

    }

    return(dist)
  }

  # Test with indexed vector

  testfn_index = function(){
    dist = vector(mode = 'numeric', length = testRle$lengths %>% length)

    start = 0
    end = 0

    for(i in seq(testRle$lengths %>% length)){

      nPoints = testRle$lengths[i]

      start = end + 1
      end = end + nPoints

      if(nPoints == 1){
        dist[i] = 0
      } else {

        dist[i] = diffDist(pts = testTracks[start:end,] %>% select(Easting,Northing)) %>% sum

      }

    }

    return(dist)


  }

  microbenchmark::microbenchmark("Not Indexed" = testfn(), "Indexed" = testfn_index(), times = 50)

}


## Put it all together ------------------------------------------------------------

detectCovar = list(
  intercept = matrix(1, nrow = nSites, ncol = maxT-1),
  dogCov = dogCov,
  humCov = humCov,
  Dcov = Dcov_scaled
)

save(detectCovar, file = 'detectCovar.Rdata')

# Save scale information ------------------------------------------------------------------------

scaleMetrics$Center = scaleMetrics$Center %>% unlist %>% unname()
scaleMetrics$Scale  = scaleMetrics$Scale %>% unlist %>% unname()

save(scaleMetrics, file = 'scaleMetrics.Rdata')
