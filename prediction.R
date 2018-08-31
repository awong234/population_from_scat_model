# Take outputs from directory, and using appropriate covariate layers, predict moose abundance

# Note; prediction will take place on all parts of the Adirondacks excluding the following areas:

# Querying the APA dataset, the following masks out the uninhabitable areas.
# LCCD_Def.LCCD_DEF IN ( 'Hamlet', 'Industrial Use', 'Intensive Use', 'Moderate Intensity', 'Water')

# Setup ---------------

library(rgeos)
library(sp)
library(ggplot2)

source('functions.R')

skip = T


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# -------------------- Obtain 1000m grid, mask by uninhabitable -----------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Simple method - reference points spaced x meters apart, and remove points falling within uninhabitable areas

# Complex method - reference grid to uninhabitable areas, and either remove squares intersecting, or clip them, keeping track of square area.

# Load adk boundary

rgdal::ogrListLayers('spatCov/adkbound/adkbound.shp')

adkbound = rgdal::readOGR(dsn = 'spatCov/adkbound/adkbound.shp', layer = 'adkbound')

# Load uninhabitable areas

uhm = rgdal::readOGR(dsn = 'spatCov/uninhabitable_mask/uninhabitable_mask.shp')

# Use complex option.
# Build grid, promote to SpatialGridDataFrame. Clip to adbound, then to uhm, and keep track of area. 
# Use gIntersection from `rgeos`

predict_grid = makegrid(x = adkbound, cellsize = 1000) %>% rename(x = x1, y = x2)
coordinates(predict_grid) = ~x + y
proj4string(predict_grid) = proj4string(adkbound)
# Promote to spatialgriddataframe
gridded(predict_grid) = TRUE
predict_grid = as(predict_grid, "SpatialGrid")

# Convert to raster to mask by polygons
predict_grid = raster::raster(predict_grid)
predict_grid[] = 1
predict_grid = raster::mask(x = predict_grid, mask = adkbound)

predict_grid_raster = raster::mask(x = predict_grid, mask = uhm, inverse = T)

# Convert back to grid

predict_grid = as(predict_grid_raster, 'SpatialPointsDataFrame')
predict_grid = predict_grid[!is.na(predict_grid@data$layer),]
gridded(predict_grid) = T
summary(predict_grid)
plot(predict_grid)

# Obtain covariates summarized to prediction grid --------------------------------------------------------------

# Put centroid value for Northing and Easting in data slot

predict_grid@data$Northing = coordinates(predict_grid)[,2]
predict_grid@data$Easting  = coordinates(predict_grid)[,1]

test = predict_grid %>% data.frame

which(!test$Northing == test$y)
all(test$Easting  == test$x)

# Format elev raster ----------------------------------------------

elev = raster::raster('spatCov/Extract_Elev1/Extract_Elev1.tif')
raster::image(elev)

system.time({
elevSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = elev, method = 'mean')
})

# There are 8 NA's. I could remove those pixels, or inherit data from nearest point.
elevSummary$rasterSummary %>% summary

predict_grid@data$Elevation = elevSummary$rasterSummary

# Inherit from nearest neighbors

predict_grid = imputeMissing(spdf = predict_grid, dataCol = "Elevation")

nadist = fields::rdist(predict_grid@coords[!complete.cases(predict_grid@data),], predict_grid@coords)
nadist_order = nadist %>% apply(MARGIN = 1, FUN = order) %>% t

nineSmall = nadist_order[,2:10]

predict_grid[nineSmall[1,],] %>% data.frame %>% 
  ggplot() + 
  geom_tile(aes(x = x, y = y, fill = Elevation)) + 
  scale_fill_continuous(na.value = 'red')

predict_grid[nineSmall[1,],][['Elevation']] %>% mean(na.rm = T)

predict_grid[nadist_order[1,1],'Elevation']

if(!skip){
  
  ggplot() + 
    geom_tile(data = predict_grid %>% data.frame, aes(x = x, y = y, fill = Elevation)) + 
    scale_fill_gradient(limits = c(0,1600), na.value = 'red') + 
    coord_equal()
  
  # Or use resample - essentially the same thing, but 
  
  test = raster::resample(elev, predict_grid_raster)
  
  test_sp = as(test, "SpatialPixelsDataFrame")
  test_df = test_sp %>% data.frame %>% rename(Elevation = Extract_Elev1)
  
  test_df$Elevation %>% summary
  
  ggplot() + 
    geom_tile(data = test_df, aes(x = x, y = y, fill = Elevation)) + 
    scale_fill_gradient(limits = c(0,1600)) + 
    coord_equal()
  
}

# Format hwy raster ----------------------------------------------

hwy = raster::raster('spatCov/KernelD_Highways/KernelD_Highways.tif')
raster::image(hwy)

system.time({
  hwySummary = summarizeRastFromGrid(grids = list(predict_grid), raster = hwy, method = 'mean')
})

predict_grid@data$Highway = hwySummary$rasterSummary

# Format minRoad raster ----------------------------------------------

minRoad = raster::raster('spatCov/KernelD_localRoads/KernelD_localRoads.tif')
raster::image(minRoad)

system.time({
  minRoadSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = minRoad, method = 'mean')
})

predict_grid@data$MinorRoad = minRoadSummary$rasterSummary

# Format habitat raster ----------------------------------------------

hab = raster::raster('spatCov/tnc_habs_clip_proj_reclass/tnc_habs_clip_proj_reclass.tif')
# Doesn't seem to work; try resample
# system.time({
#   habSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = hab, method = 'mode')
# })

habSummary = raster::resample(hab, predict_grid_raster, method = 'ngb')
habSummary = raster::mask(habSummary, predict_grid_raster)

if(!skip){
  
  # Plot to see that it's okay
  
  raster::getValues(habSummary) %>% str
  
  test_sp = as(habSummary, "SpatialPixelsDataFrame")
  test_df = test_sp %>% data.frame %>% rename(Habitat = tnc_habs_clip_proj_reclass)
  
  test_df$Habitat %>% table
  
  ggplot() + 
    geom_tile(data = test_df, aes(x = x, y = y, fill = Habitat %>% factor)) + 
    coord_equal()
  
  
}

# One-hot habitat raster ----------------------------------------------

# Remove NA values, which are masked out; see that it is the same length as predict_grid@data %>% nrow
(raster::getValues(habSummary)[!is.na(raster::getValues(habSummary))] %>% length) == nrow(predict_grid@data)

# Plot it to see if it's the same - yep looks good. 

if(!skip){
  
  coords = predict_grid@coords
  data   = raster::getValues(habSummary)[!is.na(raster::getValues(habSummary))]
  
  toPlot = data.frame(coords, data)
  
  toPlot %>% 
    ggplot() + 
    geom_tile(aes(x = x, y = y, fill = data %>% factor)) + coord_equal()
  
  
}

habData = raster::getValues(habSummary)[!is.na(raster::getValues(habSummary))]

habDummy = cbind(
  habData == 0 | habData == 5271,
  habData == 1,
  habData == 2,
  habData == 3,
  habData == 4
) * 1

colnames(habDummy) = c("Other", "Conifer", "Deciduous", "Mixed", "Wetland")

habDummy = data.frame(habDummy)

predict_grid@data$Conifer = habDummy$Conifer
predict_grid@data$Deciduous = habDummy$Deciduous
predict_grid@data$Mixed = habDummy$Mixed
predict_grid@data$Wetland = habDummy$Wetland

if(!skip){
  ggplot(predict_grid %>% data.frame) + 
    geom_tile(aes(x = x, y = y, fill = Conifer %>% factor)) + 
    coord_equal() + theme_bw()
}

# Scale covariates by DATA scaling factors - obtained from dataFormatting.R ---------------------------------------------------------

load('scaleMetrics.Rdata')

predict_grid_scaled = predict_grid

for(covariate in scaleMetrics$Covariate[1:5]){
  
  try({
    predict_grid_scaled@data[[covariate]] = scale(predict_grid_scaled@data[[covariate]], 
                                                center = scaleMetrics %>% filter(Covariate == covariate) %>% pull(Center),
                                                scale  = scaleMetrics %>% filter(Covariate == covariate) %>% pull(Scale)
  ) %>% as.numeric
  })
  
}

# Convert back to make sure it's right - yes it's right, within rounding error.

for(covariate in scaleMetrics$Covariate[1:5]){
  print(
  ((((predict_grid_scaled@data[[covariate]] * scaleMetrics %>% filter(Covariate == covariate) %>% pull(Scale)) + 
      scaleMetrics %>% filter(Covariate == covariate) %>% pull(Center)) - predict_grid@data[[covariate]]) < 0.01) %>% all
  )
}

# Save prediction items -------------------------------------------------

save(predict_grid, file = paste0('predict_grid_', raster::res(predict_grid_raster)[1], '.Rdata'))
save(predict_grid_scaled, file = paste0('predict_grid_scaled_', raster::res(predict_grid_raster)[1], '.Rdata'))
save(predict_grid_raster, file = paste0('predict_grid_raster_', raster::res(predict_grid_raster)[1], '.Rdata'))

# Load prediction items -------------------------------------------------

# Start here after setting up prediction grid

load('predict_grid_1000.Rdata')
load('predict_grid_raster_1000.Rdata')
load('predict_grid_scaled_1000.Rdata')


# Get area of cells ----------------------------------------

cellDim = raster::res(predict_grid_raster)[1]

area = cellDim^2

analysisArea = 50*50

# Output folders ---------------------------------------------------

folders = dir(path = 'modelOutputs/', full.names = T)

files_in_folders = lapply(X = folders, FUN = function(x){list.files(path = x, full.names = T)})
files_info = lapply(files_in_folders, FUN = function(x){
  dates = file.info(x)$mtime
  tempDf = data.frame('paths' = x, dates) %>% arrange(dates)
  return(tempDf)
})

# Only care about the latest file. 

latest_files = lapply(files_info, FUN = function(x){tail(x, n = 1)}) %>% do.call(what = rbind.data.frame)
latest_files$modName = c("null", "cont", "crit", "dcov")

### Perform prediction -----------------------------------------------------------------------------------------

# Defecation rates --------------------------------------------------------------------------------------------

defecationRates = data.frame(reference = c(rep("Miquelle", 4),
                                           rep("Joyal&Ricard", 3)),
                             mean = c(10.9, 19, 13, 11.2,
                                      12.3, 13.5, 12.1),
                             se = c(0.5, 0.5, 0.7, 1,
                                    NA,NA,NA),
                             sd = c(NA, NA, NA, NA,
                                    5.8, 6.3, 3.9),
                             N = c(22, 8, 22, 21, 
                                   38,  38,  38),
                             ageClass = c('yearling', 'calf', 'yearling', 'adult',
                                          'calf', 'adult', 'adult'),
                             sex = c('mf', 'm', 'm', 'mf',
                                     'mf', 'f', 'f'),
                             status = c('captive', 'captive', 'captive', 'free-range',
                                        'free-range', 'free-range', 'free-range'),
                             season = c(rep('summer', 4),
                                        rep('winter', 3))
)

miquelleRows = defecationRates$reference == 'Miquelle'

defecationRates$sd[miquelleRows] = with(defecationRates, expr = {se * sqrt(N)})[miquelleRows]


# Predict NULL ------------------------------------------------------------------------------------------------

# Using nimble output for null model because it's the best one
latest_files$paths[1] = 'modelOutputs/null/out_null_nimble_2018-08-19.Rdata'

# Object is named 'samples'
load(latest_files$paths[1] %>% as.character())

samples_mcmc = samples %>% lapply(X = ., FUN = function(x){coda::mcmc(x)}) %>% coda::as.mcmc.list()
samples_mcmc_summ = samples_mcmc %>% summary

meanTheta = samples_mcmc_summ$statistics[3]
ciTheta = samples_mcmc_summ$quantiles[c(3,15)]


# New function 

estimates = list(
  meanTheta = meanTheta,
  ciTheta = ciTheta
)

outList = summarizeOutput(predict_grid = predict_grid_scaled, estimates = estimates)

save(outList, file = 'predictionOutput/null/prediction_NULL.Rdata')

# Predict DCOV ------------------------------------------------------------------------------------------------

# Named 'output'
load(latest_files$paths[4] %>% as.character())

meanTheta = output$summary[1]
ciTheta = c(output$summary[1,3], output$summary[1,7])

totalArea = area * nrow(predict_grid_scaled)

(MooseAbundanceMean = (exp(meanTheta) / mean(defecationRates$mean)) * (totalArea / analysisArea))
MooseAbundanceCI = (exp(ciTheta) / mean(defecationRates$mean)) * (totalArea / analysisArea)

(MooseAbundance = data.frame("Estimate" = MooseAbundanceMean, "Lower 95%" = MooseAbundanceCI[1], "Upper 95%" = MooseAbundanceCI[2]))

# Bootstrap defecation rate 

bootstrapDef = rnorm(n = 1e6, mean = defecationRates$mean %>% mean, sd = defecationRates$mean %>% sd)
bootstrapDef = bootstrapDef[bootstrapDef > 0]

MooseAbundanceMeanBS = (exp(meanTheta) / bootstrapDef) * (totalArea / analysisArea)
MooseAbundanceCIBS = foreach(i = ciTheta, .combine = cbind) %do% (exp(i) / bootstrapDef) * (totalArea / analysisArea)

MooseAbundanceBS = data.frame("Estimate" = MooseAbundanceMeanBS, "Lower 95%" = MooseAbundanceCIBS[,1], "Upper 95%" = MooseAbundanceCIBS[,2])

MooseAbundanceBS %>% reshape2::melt() %>% 
  ggplot() + 
  geom_density(aes(x = value))

(finalQuantiles = MooseAbundanceBS %>% unlist %>% quantile(prob = c(0.025, 0.5, 0.95)))

outList = list(MooseAbundance, MooseAbundanceBS, finalQuantiles)

save(outList, file = 'predictionOutput/dCov/prediction_dcov.Rdata')

# Predict Continuous ------------------------------------------------------------------------------------------

# named 'output'
load(latest_files$paths[2] %>% as.character())

relOutput = output$summary[c(1,4,5,6,7),c(1,3,7)] %>% data.frame

# For every cell in the raster, predict moose abundance.

predicted_grid_theta = exp(relOutput$mean[1] + 
                             relOutput$mean[2] * predict_grid_scaled@data$Elevation +
                             relOutput$mean[3] * predict_grid_scaled@data$Highway   +
                             relOutput$mean[4] * predict_grid_scaled@data$MinorRoad +
                             relOutput$mean[5] * predict_grid_scaled@data$Northing
                           )

predicted_grid_abundance = (predicted_grid_theta / mean(defecationRates$mean)) * (area / analysisArea)

predict_grid_scaled@data$MeanAbundance = predicted_grid_abundance


ggplot() + 
  geom_tile(data = predict_grid_scaled %>% data.frame, aes(x = x, y = y, fill = MeanAbundance)) + 
  scale_fill_continuous(trans = 'log') + 
  coord_equal()
  
(sumAbundance = predict_grid_scaled@data$MeanAbundance %>% sum)

plot(predict_grid@data$Elevation, predict_grid_scaled@data$MeanAbundance)
plot(predict_grid@data$Highway, predict_grid_scaled@data$MeanAbundance)
plot(predict_grid@data$MinorRoad, predict_grid_scaled@data$MeanAbundance)
plot(predict_grid@data$Northing, predict_grid_scaled@data$MeanAbundance)


# Are we representing the covariates' populations?
load('gridCovariates.Rdata')

# Elevation - not upper elevations. Maximum is 760
data = ((gridCovariates$Elevation * scaleMetrics[3,3] ) + scaleMetrics[3,2] )

elevation_df = data.frame(variable = c(rep('data', length(data)), rep('population', nrow(predict_grid))),
                          value    = c(data, predict_grid@data$Elevation))

ggplot() + 
  geom_density(data = elevation_df, aes(x = value, color = variable))

# Highway - YES

data = ((gridCovariates$scaledHighway * scaleMetrics[4,3] ) + scaleMetrics[4,2] )

elevation_df = data.frame(variable = c(rep('data', length(data)), rep('population', nrow(predict_grid))),
                          value    = c(data, predict_grid@data$Highway))

ggplot() + 
  geom_density(data = elevation_df, aes(x = value, color = variable))

# Minor road - mostly

data = ((gridCovariates$scaledMinRoad * scaleMetrics[5,3] ) + scaleMetrics[5,2] )

elevation_df = data.frame(variable = c(rep('data', length(data)), rep('population', nrow(predict_grid))),
                          value    = c(data, predict_grid@data$MinorRoad))

ggplot() + 
  geom_density(data = elevation_df, aes(x = value, color = variable))

# Northing - YES

data = ((gridCovariates$Northing * scaleMetrics[1,3] ) + scaleMetrics[1,2] )

elevation_df = data.frame(variable = c(rep('data', length(data)), rep('population', nrow(predict_grid))),
                          value    = c(data, predict_grid@data$Northing))

ggplot() + 
  geom_density(data = elevation_df, aes(x = value, color = variable))


# What if we lop off the outlier elevation? Max observed is 2.947981 

# This results in an estimate of 390 animals. Including high peaks in the
# prediction implies that we have 1105 - 390 = 715 animals living in the high
# peaks. Not particularly sensible.

# First option - we don't care about the number of moose in the high peaks
# because we can't or won't effectively manage it. Eliminate those areas that
# are greater than the 95th %ile, which are outside our observed data anyway.
# This is a value of 3.16 in scaled elevation.

elevIndex = (predict_grid_scaled@data$Elevation < quantile(predict_grid_scaled@data$Elevation, prob = c(.99))) %>% as.logical()

predict_grid_scaled_test = predict_grid_scaled[elevIndex, ]

predicted_grid_theta = exp(relOutput$mean[1] + 
                             relOutput$mean[2] * predict_grid_scaled_test@data$Elevation +
                             relOutput$mean[3] * predict_grid_scaled_test@data$Highway   +
                             relOutput$mean[4] * predict_grid_scaled_test@data$MinorRoad +
                             relOutput$mean[5] * predict_grid_scaled_test@data$Northing
)

predicted_grid_abundance = (predicted_grid_theta / mean(defecationRates$mean)) * (area / analysisArea)

predict_grid_scaled_test@data$MeanAbundance = predicted_grid_abundance


ggplot() + 
  geom_tile(data = predict_grid_scaled_test %>% data.frame, aes(x = x, y = y, fill = MeanAbundance)) + 
  # scale_fill_continuous(trans = 'log') + 
  geom_density_2d(data = scatsReferenced, aes(x = ScatEasting, y = ScatNorthing), alpha = .25, color = 'green') +
  geom_point(data = tracks_points %>% data.frame %>% group_by(Site) %>% sample_n(size = 1), aes(x = Easting, y = Northing), color = 'red', shape = 1, alpha = 0.5) + 
  # geom_point(data = scatsReferenced, aes(x = ScatEasting, y = ScatNorthing), color = 'white', alpha = 0.01) + 
  coord_equal() + theme_bw() + 
  theme(panel.background = element_rect(fill = 'gray5'),
        plot.background  = element_rect(fill = 'gray5'),
        legend.background = element_rect(fill = 'gray5'),
        legend.text = element_text(color = 'gray50'),
        legend.title = element_text(color = 'gray50'),
        line = element_blank())

plot_ly(data = predict_grid_scaled_test %>% data.frame %>% select(x,y,MeanAbundance)) %>% 
  add_heatmap(x = ~x, y = ~y, z = ~MeanAbundance) %>% 
  add_markers(data = tracks_points %>% data.frame %>% group_by(Site) %>% sample_n(size = 1), x = ~Easting, y = ~Northing)

(sumAbundance = predict_grid_scaled_test@data$MeanAbundance %>% sum)

popElev = (predict_grid_scaled_test@data$Elevation * scaleMetrics[3,3]) + scaleMetrics[3,2]

plot(popElev, predict_grid_scaled_test@data$MeanAbundance)

# If this is the case, we should also recalculate the null model. There are about 528 here, vs. 555. Not a major difference.

meanTheta = samples_mcmc_summ$statistics[3]
ciTheta = samples_mcmc_summ$quantiles[c(3,15)]

totalArea = area * nrow(predict_grid_scaled_test)

(MooseAbundanceMean = (exp(meanTheta) / mean(defecationRates$mean)) * (totalArea / analysisArea))


# Second option
# We impute the mean elevation there; saying that we don't know how many animals are up there, but it is probably an average number unrelated to elevation
# Answer is also more reasonable; 420 moose total, perhaps 30 moose in the high peaks.

predict_grid_scaled_test = predict_grid_scaled

elevIndex = (predict_grid_scaled@data$Elevation < quantile(predict_grid_scaled@data$Elevation, prob = c(.95))) %>% as.logical()

predict_grid_scaled_test@data$Elevation[!elevIndex, ] = 0

predicted_grid_theta = exp(relOutput$mean[1] + 
                             relOutput$mean[2] * predict_grid_scaled_test@data$Elevation +
                             relOutput$mean[3] * predict_grid_scaled_test@data$Highway   +
                             relOutput$mean[4] * predict_grid_scaled_test@data$MinorRoad +
                             relOutput$mean[5] * predict_grid_scaled_test@data$Northing
)

predicted_grid_abundance = (predicted_grid_theta / mean(defecationRates$mean)) * (area / analysisArea)

predict_grid_scaled_test@data$MeanAbundance = predicted_grid_abundance


ggplot() + 
  geom_tile(data = predict_grid_scaled_test %>% data.frame, aes(x = x, y = y, fill = MeanAbundance)) + 
  scale_fill_continuous(trans = 'log') + 
  coord_equal()

(sumAbundance = predict_grid_scaled_test@data$MeanAbundance %>% sum)




# Predict Critical model --------------------------------------------------------------------------------------

# No sensible results yet.


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# -------------------- Obtain 100m grid, mask by uninhabitable -----------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# Simple method - reference points spaced x meters apart, and remove points falling within uninhabitable areas

# Complex method - reference grid to uninhabitable areas, and either remove squares intersecting, or clip them, keeping track of square area.

# Load adk boundary

rgdal::ogrListLayers('spatCov/adkbound/adkbound.shp')

adkbound = rgdal::readOGR(dsn = 'spatCov/adkbound/adkbound.shp', layer = 'adkbound')

# Load uninhabitable areas

uhm = rgdal::readOGR(dsn = 'spatCov/uninhabitable_mask/uninhabitable_mask.shp')

# Use complex option.
# Build grid, promote to SpatialGridDataFrame. Clip to adbound, then to uhm, and keep track of area. 
# Use gIntersection from `rgeos`

predict_grid = makegrid(x = adkbound, cellsize = 100) %>% rename(x = x1, y = x2)
coordinates(predict_grid) = ~x + y
proj4string(predict_grid) = proj4string(adkbound)
# Promote to spatialgriddataframe
gridded(predict_grid) = TRUE
predict_grid = as(predict_grid, "SpatialGrid")

# Convert to raster to mask by polygons
predict_grid = raster::raster(predict_grid)
predict_grid[] = 1
predict_grid = raster::mask(x = predict_grid, mask = adkbound)

predict_grid_raster = raster::mask(x = predict_grid, mask = uhm, inverse = T)

# Convert back to grid

predict_grid = as(predict_grid_raster, 'SpatialPointsDataFrame')
predict_grid = predict_grid[!is.na(predict_grid@data$layer),]
gridded(predict_grid) = T
summary(predict_grid)
plot(predict_grid)

# Obtain covariates summarized to prediction grid --------------------------------------------------------------

# Put centroid value for Northing and Easting in data slot

predict_grid@data$Northing = coordinates(predict_grid)[,2]
predict_grid@data$Easting  = coordinates(predict_grid)[,1]

test = predict_grid %>% data.frame

which(!test$Northing == test$y)
all(test$Easting  == test$x)

# Format elev raster ----------------------------------------------

elev = raster::raster('spatCov/Extract_Elev1/Extract_Elev1.tif')
raster::image(elev)

system.time({
  elevSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = elev, method = 'mean')
})

# There may be NA's. I could remove those pixels, or inherit data from nearest points.
elevSummary$rasterSummary %>% summary

predict_grid@data$Elevation = elevSummary$rasterSummary

# Inherit from nearest neighbors

predict_grid = imputeMissing(spdf = predict_grid, dataCol = "Elevation")

nadist = fields::rdist(predict_grid@coords[!complete.cases(predict_grid@data),], predict_grid@coords)
nadist_order = nadist %>% apply(MARGIN = 1, FUN = order) %>% t

nineSmall = nadist_order[,2:10]

predict_grid[nineSmall[1,],] %>% data.frame %>% 
  ggplot() + 
  geom_tile(aes(x = x, y = y, fill = Elevation)) + 
  scale_fill_continuous(na.value = 'red')

predict_grid[nineSmall[1,],][['Elevation']] %>% mean(na.rm = T)

predict_grid[nadist_order[1,1],'Elevation']

if(!skip){
  
  ggplot() + 
    geom_tile(data = predict_grid %>% data.frame, aes(x = x, y = y, fill = Elevation)) + 
    scale_fill_gradient(limits = c(0,1600), na.value = 'red') + 
    coord_equal()
  
  # Or use resample - essentially the same thing, but 
  
  test = raster::resample(elev, predict_grid_raster)
  
  test_sp = as(test, "SpatialPixelsDataFrame")
  test_df = test_sp %>% data.frame %>% rename(Elevation = Extract_Elev1)
  
  test_df$Elevation %>% summary
  
  ggplot() + 
    geom_tile(data = test_df, aes(x = x, y = y, fill = Elevation)) + 
    scale_fill_gradient(limits = c(0,1600)) + 
    coord_equal()
  
}

# Format hwy raster ----------------------------------------------

hwy = raster::raster('spatCov/KernelD_Highways/KernelD_Highways.tif')
raster::image(hwy)

system.time({
  hwySummary = summarizeRastFromGrid(grids = list(predict_grid), raster = hwy, method = 'mean')
})

predict_grid@data$Highway = hwySummary$rasterSummary

# Format minRoad raster ----------------------------------------------

minRoad = raster::raster('spatCov/KernelD_localRoads/KernelD_localRoads.tif')
raster::image(minRoad)

system.time({
  minRoadSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = minRoad, method = 'mean')
})

predict_grid@data$MinorRoad = minRoadSummary$rasterSummary

# Format habitat raster ----------------------------------------------

hab = raster::raster('spatCov/tnc_habs_clip_proj_reclass/tnc_habs_clip_proj_reclass.tif')
# Doesn't seem to work; try resample
# system.time({
#   habSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = hab, method = 'mode')
# })

habSummary = raster::resample(hab, predict_grid_raster, method = 'ngb')
habSummary = raster::mask(habSummary, predict_grid_raster)

if(!skip){
  
  # Plot to see that it's okay
  
  raster::getValues(habSummary) %>% str
  
  test_sp = as(habSummary, "SpatialPixelsDataFrame")
  test_df = test_sp %>% data.frame %>% rename(Habitat = tnc_habs_clip_proj_reclass)
  
  test_df$Habitat %>% table
  
  ggplot() + 
    geom_tile(data = test_df, aes(x = x, y = y, fill = Habitat %>% factor)) + 
    coord_equal()
  
  
}

# One-hot habitat raster ----------------------------------------------

# Remove NA values, which are masked out; see that it is the same length as predict_grid@data %>% nrow
(raster::getValues(habSummary)[!is.na(raster::getValues(habSummary))] %>% length) == nrow(predict_grid@data)

# Plot it to see if it's the same - yep looks good. 

if(!skip){
  
  coords = predict_grid@coords
  data   = raster::getValues(habSummary)[!is.na(raster::getValues(habSummary))]
  
  toPlot = data.frame(coords, data)
  
  toPlot %>% 
    ggplot() + 
    geom_tile(aes(x = x, y = y, fill = data %>% factor)) + coord_equal()
  
  
}

habData = raster::getValues(habSummary)[!is.na(raster::getValues(habSummary))]

habDummy = cbind(
  habData == 0 | habData == 5271,
  habData == 1,
  habData == 2,
  habData == 3,
  habData == 4
) * 1

colnames(habDummy) = c("Other", "Conifer", "Deciduous", "Mixed", "Wetland")

habDummy = data.frame(habDummy)

predict_grid@data$Conifer = habDummy$Conifer
predict_grid@data$Deciduous = habDummy$Deciduous
predict_grid@data$Mixed = habDummy$Mixed
predict_grid@data$Wetland = habDummy$Wetland

if(!skip){
  ggplot(predict_grid %>% data.frame) + 
    geom_tile(aes(x = x, y = y, fill = Conifer %>% factor)) + 
    coord_equal() + theme_bw()
}

# Scale covariates by DATA scaling factors - obtained from dataFormatting.R ---------------------------------------------------------

load('scaleMetrics.Rdata')

predict_grid_scaled = predict_grid

for(covariate in scaleMetrics$Covariate[1:5]){
  
  try({
    predict_grid_scaled@data[[covariate]] = scale(predict_grid_scaled@data[[covariate]], 
                                                  center = scaleMetrics %>% filter(Covariate == covariate) %>% pull(Center),
                                                  scale  = scaleMetrics %>% filter(Covariate == covariate) %>% pull(Scale)
    ) %>% as.numeric
  })
  
}

# Convert back to make sure it's right - yes it's right, within rounding error.

for(covariate in scaleMetrics$Covariate[1:5]){
  print(
    ((((predict_grid_scaled@data[[covariate]] * scaleMetrics %>% filter(Covariate == covariate) %>% pull(Scale)) + 
         scaleMetrics %>% filter(Covariate == covariate) %>% pull(Center)) - predict_grid@data[[covariate]]) < 0.01) %>% all
  )
}

# Save prediction items -------------------------------------------------

save(predict_grid, file = paste0('predict_grid_', raster::res(predict_grid_raster)[1], '.Rdata'))
save(predict_grid_scaled, file = paste0('predict_grid_scaled_', raster::res(predict_grid_raster)[1], '.Rdata'))
save(predict_grid_raster, file = paste0('predict_grid_raster_', raster::res(predict_grid_raster)[1], '.Rdata'))

# Load prediction items -------------------------------------------------

# Start here after setting up prediction grid

load('predict_grid_100.Rdata')
load('predict_grid_raster_100.Rdata')
load('predict_grid_scaled_100.Rdata')


# Get area of cells ----------------------------------------

cellDim = raster::res(predict_grid_raster)[1]

area = cellDim^2

analysisArea = 50*50

# Output folders ---------------------------------------------------

folders = dir(path = 'modelOutputs/', full.names = T)

files_in_folders = lapply(X = folders, FUN = function(x){list.files(path = x, full.names = T)})
files_info = lapply(files_in_folders, FUN = function(x){
  dates = file.info(x)$mtime
  tempDf = data.frame('paths' = x, dates) %>% arrange(dates)
  return(tempDf)
})

# Only care about the latest file. 

latest_files = lapply(files_info, FUN = function(x){tail(x, n = 1)}) %>% do.call(what = rbind.data.frame)
latest_files$modName = c("null", "cont", "crit", "dcov")