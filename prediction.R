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
# -------------------- Obtain grid, mask by uninhabitable -----------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

cellSize = 100

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

predict_grid = makegrid(x = adkbound, cellsize = cellSize) %>% rename(x = x1, y = x2)
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
dev.off()
raster::image(elev)

# Not completely necessary! Works with elevation raster and prediction grid.
# system.time({
# elevSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = elev, method = 'mean')
# })

elevSummary = raster::extract(elev, predict_grid, fun = 'mean', na.rm = T)

ggplot() + 
  geom_tile(data = test, aes(x = x, y = y, fill = Elevation)) + 
  scale_fill_gradient(limits = c(0,1600)) + 
  coord_equal()


# There are NA's. I could remove those pixels, or inherit data from nearest point.
elevSummary %>% summary

predict_grid@data$Elevation = elevSummary

# Inherit from nearest neighbors

predict_grid = imputeMissing(spdf = predict_grid, dataCol = "Elevation")

predict_grid@data$Elevation %>% is.na %>% any

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
dev.off()
raster::image(hwy)

# system.time({
#   hwySummary = summarizeRastFromGrid(grids = list(predict_grid), raster = hwy, method = 'mean')
# })

hwySummary = raster::extract(hwy, predict_grid, fun = 'mean', na.rm = T)

predict_grid@data$Highway = hwySummary

# Format minRoad raster ----------------------------------------------

minRoad = raster::raster('spatCov/KernelD_localRoads/KernelD_localRoads.tif')
dev.off()
raster::image(minRoad)

# system.time({
#   minRoadSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = minRoad, method = 'mean')
# })

minRoadSummary = raster::extract(minRoad, predict_grid, fun = 'mean', na.rm = T)

predict_grid@data$MinorRoad = minRoadSummary

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

save(predict_grid, file = paste0('predict_grid_', cellSize, '.Rdata'))
save(predict_grid_scaled, file = paste0('predict_grid_scaled_', cellSize, '.Rdata'))
save(predict_grid_raster, file = paste0('predict_grid_raster_', cellSize, '.Rdata'))

# Load prediction items -------------------------------------------------

# Start here after setting up prediction grid

cellSize = 1000

load(paste0('predict_grid_', cellSize, '.Rdata'))
load(paste0('predict_grid_raster_', cellSize, '.Rdata'))
load(paste0('predict_grid_scaled_', cellSize, '.Rdata'))


# Get area of cells ----------------------------------------

area = cellSize^2

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

# For bootstrapping purposes, fit to gamma distribution.
par = c(10,1)

gammafn = function(par){
  nll = dgamma(defecationRates$mean, shape = par[1], rate = par[2], log = T)
  joint.nll = sum(nll)
  return(-joint.nll)
}

gammafn(par)

(fit_optim = optim(par = par, fn = gammafn))

(fit_nlm = nlm(p = par, f = gammafn))

gShape = fit_nlm$estimate[1]
gRate  = fit_nlm$estimate[2]


# Predict NULL ------------------------------------------------------------------------------------------------

# Using nimble output for null model because it's the best one
latest_files$paths[1] = 'modelOutputs/null/out_null_nimble_2018-08-19.Rdata'

# Object is named 'samples'
load(latest_files$paths[1] %>% as.character())

samples_mcmc = samples %>% lapply(X = ., FUN = function(x){coda::mcmc(x)}) %>% coda::as.mcmc.list()
samples_mcmc_summ = samples_mcmc %>% summary

theta = numeric(length = 3)
theta[1] = samples_mcmc_summ$statistics[3]
theta[c(2,3)] = samples_mcmc_summ$quantiles[c(3,15)]


# New function 

outList = summarizeOutput(predict_grid = predict_grid_scaled, theta = theta, covariates = NULL)

save(outList, file = paste0('predictionOutput/null/prediction_NULL_', cellSize, 'm_all_elev.Rdata'))

# Predict DCOV ------------------------------------------------------------------------------------------------

# Named 'output'
load(latest_files$paths[4] %>% as.character())

theta = output$summary[1,c(1,3,7)]

outList = summarizeOutput(predict_grid = predict_grid_scaled, theta = theta, covariates = NULL)

save(outList, file = 'predictionOutput/dCov/prediction_dcov_1000m_all_elev.Rdata')

# Predict Continuous ------------------------------------------------------------------------------------------

# named 'output'
load(latest_files$paths[2] %>% as.character())

relOutput = output$summary[c(1,4,5,6,7),c(1,3,7)] %>% data.frame %>% rename(Lower = `X2.5.`, Upper = `X97.5.`)

# For every cell in the raster, predict moose abundance.

predicted_grid_theta = foreach(col = 1:ncol(relOutput), 
                               .combine = cbind, 
                               .final = function(x){attr(x = x, which = 'dimnames') = list(NULL, names(relOutput)); return(x)}) %do% {
  exp(relOutput[1,col] + 
        relOutput[2,col] * predict_grid_scaled@data$Elevation +
        relOutput[3,col] * predict_grid_scaled@data$Highway   +
        relOutput[4,col] * predict_grid_scaled@data$MinorRoad +
        relOutput[5,col] * predict_grid_scaled@data$Northing
      )
                               }

outList = summarizeOutput(predict_grid = predict_grid_scaled, 
                          theta = predicted_grid_theta, 
                          covariates = c("Northing", "Elevation", "Highway", "MinorRoad"))

save(outList, file = 'predictionOutput/cont/prediction_cont_1000m_all_elev.Rdata')

load('predictionOutput/cont/prediction_cont_1000m_all_elev.Rdata')

extract(outList)

predict_grid_scaled@data$MeanAbundance = gridEstimates[,1]

ggplot() + 
  geom_tile(data = predict_grid_scaled %>% data.frame, aes(x = x, y = y, fill = MeanAbundance)) + 
  scale_fill_continuous(trans = 'log') + 
  coord_equal()

if(!skip){
  
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
  
  
}
  

# What if we lop off the outlier elevation? Max observed is 2.947981 

# This results in an estimate of 390 animals. Including high peaks in the
# prediction implies that we have 1105 - 390 = 715 animals living in the high
# peaks. Not particularly sensible.

# First option - we don't care about the number of moose in the high peaks
# because we can't or won't effectively manage it. Eliminate those areas that
# are greater than the 95th %ile, which are outside our observed data anyway.
# This is a value of 3.16 in scaled elevation.

# 99% option -------------------------------------------------------------------------------

elevIndex = (predict_grid_scaled@data$Elevation < quantile(predict_grid_scaled@data$Elevation, prob = c(.99))) %>% as.logical()

pr_grid_sc_elev = predict_grid_scaled[elevIndex, ]

save(pr_grid_sc_elev, file = 'predict_grid_scaled_1000_elev.Rdata')


predicted_grid_theta = foreach(col = 1:ncol(relOutput), 
                               .combine = cbind, 
                               .final = function(x){attr(x = x, which = 'dimnames') = list(NULL, names(relOutput)); return(x)}) %do% {
                                 exp(relOutput[1,col] + 
                                       relOutput[2,col] * pr_grid_sc_elev@data$Elevation +
                                       relOutput[3,col] * pr_grid_sc_elev@data$Highway   +
                                       relOutput[4,col] * pr_grid_sc_elev@data$MinorRoad +
                                       relOutput[5,col] * pr_grid_sc_elev@data$Northing
                                 )
                               }


outList = summarizeOutput(predict_grid = pr_grid_sc_elev, 
                          theta = predicted_grid_theta, 
                          covariates = c("Northing", "Elevation", "Highway", "MinorRoad"))

save(outList, file = 'predictionOutput/cont/prediction_cont_1000m_elev_trim_99.Rdata')

load('predictionOutput/cont/prediction_cont_1000m_elev_trim.Rdata')

extract(outList)

pr_grid_sc_elev@data$MeanAbundance = gridEstimates[,1]

ggplot() + 
  geom_tile(data = pr_grid_sc_elev %>% data.frame, aes(x = x, y = y, fill = MeanAbundance)) + 
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

plot_ly(data = pr_grid_sc_elev %>% data.frame %>% select(x,y,MeanAbundance)) %>% 
  add_heatmap(x = ~x, y = ~y, z = ~MeanAbundance) %>% 
  add_markers(data = tracks_points %>% data.frame %>% group_by(Site) %>% sample_n(size = 1), x = ~Easting, y = ~Northing)

(sumAbundance = pr_grid_sc_elev@data$MeanAbundance %>% sum)

popElev = (pr_grid_sc_elev@data$Elevation * scaleMetrics[3,3]) + scaleMetrics[3,2]

plot(popElev, pr_grid_sc_elev@data$MeanAbundance)

# If this is the case, we should also recalculate the null model. There are about 528 here, vs. 555. Not a major difference.

meanTheta = samples_mcmc_summ$statistics[3]
ciTheta = samples_mcmc_summ$quantiles[c(3,15)]

theta = c(samples_mcmc_summ$statistics[3], samples_mcmc_summ$quantiles[c(3,15)])

outList = summarizeOutput(predict_grid = pr_grid_sc_elev, theta = theta, covariates = NULL)

save(outList, file = 'predictionOutput/null/prediction_NULL_1000m_elev_trim_99.Rdata')

# Second option
# We impute the mean elevation there; saying that we don't know how many animals are up there, but it is probably an average number unrelated to elevation
# Answer is also more reasonable; 420 moose total, perhaps 30 moose in the high peaks.

pr_grid_sc_elev = predict_grid_scaled

elevIndex = (predict_grid_scaled@data$Elevation < quantile(predict_grid_scaled@data$Elevation, prob = c(.99))) %>% as.logical()

pr_grid_sc_elev@data$Elevation[!elevIndex ] = 0

predicted_grid_theta = foreach(col = 1:ncol(relOutput), 
                               .combine = cbind, 
                               .final = function(x){attr(x = x, which = 'dimnames') = list(NULL, names(relOutput)); return(x)}) %do% {
                                 exp(relOutput[1,col] + 
                                       relOutput[2,col] * pr_grid_sc_elev@data$Elevation +
                                       relOutput[3,col] * pr_grid_sc_elev@data$Highway   +
                                       relOutput[4,col] * pr_grid_sc_elev@data$MinorRoad +
                                       relOutput[5,col] * pr_grid_sc_elev@data$Northing
                                 )
                               }

outList = summarizeOutput(predict_grid = pr_grid_sc_elev, 
                          theta = predicted_grid_theta, 
                          covariates = c("Northing", "Elevation", "Highway", "MinorRoad"))

save(outList, file = 'predictionOutput/cont/prediction_cont_1000m_elev_mean_99.Rdata')

# predicted_grid_abundance = (predicted_grid_theta / mean(defecationRates$mean)) * (area / analysisArea)
# 
# pr_grid_sc_elev@data$MeanAbundance = predicted_grid_abundance

pr_grid_sc_elev@data$MeanAbundance = outList$gridEstimates[,1]

ggplot() + 
  geom_tile(data = pr_grid_sc_elev %>% data.frame, aes(x = x, y = y, fill = MeanAbundance)) + 
  # scale_fill_continuous(trans = 'log') + 
  coord_equal()

(sumAbundance = pr_grid_sc_elev@data$MeanAbundance %>% sum)

# 95% option -------------------------------------------------------------------------------

elevIndex = (predict_grid_scaled@data$Elevation < quantile(predict_grid_scaled@data$Elevation, prob = c(.95))) %>% as.logical()

pr_grid_sc_elev = predict_grid_scaled[elevIndex, ]

save(pr_grid_sc_elev, file = 'predict_grid_scaled_1000_elev.Rdata')


predicted_grid_theta = foreach(col = 1:ncol(relOutput), 
                               .combine = cbind, 
                               .final = function(x){attr(x = x, which = 'dimnames') = list(NULL, names(relOutput)); return(x)}) %do% {
                                 exp(relOutput[1,col] + 
                                       relOutput[2,col] * pr_grid_sc_elev@data$Elevation +
                                       relOutput[3,col] * pr_grid_sc_elev@data$Highway   +
                                       relOutput[4,col] * pr_grid_sc_elev@data$MinorRoad +
                                       relOutput[5,col] * pr_grid_sc_elev@data$Northing
                                 )
                               }


outList = summarizeOutput(predict_grid = pr_grid_sc_elev, 
                          theta = predicted_grid_theta, 
                          covariates = c("Northing", "Elevation", "Highway", "MinorRoad"))

save(outList, file = 'predictionOutput/cont/prediction_cont_1000m_elev_trim_95.Rdata')

extract(outList)

pr_grid_sc_elev@data$MeanAbundance = gridEstimates[,1]

ggplot() + 
  geom_tile(data = pr_grid_sc_elev %>% data.frame, aes(x = x, y = y, fill = MeanAbundance)) + 
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

plot_ly(data = pr_grid_sc_elev %>% data.frame %>% select(x,y,MeanAbundance)) %>% 
  add_heatmap(x = ~x, y = ~y, z = ~MeanAbundance) %>% 
  add_markers(data = tracks_points %>% data.frame %>% group_by(Site) %>% sample_n(size = 1), x = ~Easting, y = ~Northing)

(sumAbundance = pr_grid_sc_elev@data$MeanAbundance %>% sum)

popElev = (pr_grid_sc_elev@data$Elevation * scaleMetrics[3,3]) + scaleMetrics[3,2]

plot(popElev, pr_grid_sc_elev@data$MeanAbundance)

# If this is the case, we should also recalculate the null model. There are about 528 here, vs. 555. Not a major difference.

meanTheta = samples_mcmc_summ$statistics[3]
ciTheta = samples_mcmc_summ$quantiles[c(3,15)]

theta = c(samples_mcmc_summ$statistics[3], samples_mcmc_summ$quantiles[c(3,15)])

outList = summarizeOutput(predict_grid = pr_grid_sc_elev, theta = theta, covariates = NULL)

save(outList, file = 'predictionOutput/null/prediction_NULL_1000m_elev_trim_95.Rdata')

# Second option
# We impute the mean elevation there; saying that we don't know how many animals are up there, but it is probably an average number unrelated to elevation
# Answer is also more reasonable; 420 moose total, perhaps 30 moose in the high peaks.

pr_grid_sc_elev = predict_grid_scaled

elevIndex = (predict_grid_scaled@data$Elevation < quantile(predict_grid_scaled@data$Elevation, prob = c(.95))) %>% as.logical()

pr_grid_sc_elev@data$Elevation[!elevIndex ] = 0

predicted_grid_theta = foreach(col = 1:ncol(relOutput), 
                               .combine = cbind, 
                               .final = function(x){attr(x = x, which = 'dimnames') = list(NULL, names(relOutput)); return(x)}) %do% {
                                 exp(relOutput[1,col] + 
                                       relOutput[2,col] * pr_grid_sc_elev@data$Elevation +
                                       relOutput[3,col] * pr_grid_sc_elev@data$Highway   +
                                       relOutput[4,col] * pr_grid_sc_elev@data$MinorRoad +
                                       relOutput[5,col] * pr_grid_sc_elev@data$Northing
                                 )
                               }

outList = summarizeOutput(predict_grid = pr_grid_sc_elev, 
                          theta = predicted_grid_theta, 
                          covariates = c("Northing", "Elevation", "Highway", "MinorRoad"))

save(outList, file = 'predictionOutput/cont/prediction_cont_1000m_elev_mean_95.Rdata')

# predicted_grid_abundance = (predicted_grid_theta / mean(defecationRates$mean)) * (area / analysisArea)
# 
# pr_grid_sc_elev@data$MeanAbundance = predicted_grid_abundance

pr_grid_sc_elev@data$MeanAbundance = outList$gridEstimates[,1]

ggplot() + 
  geom_tile(data = pr_grid_sc_elev %>% data.frame, aes(x = x, y = y, fill = MeanAbundance)) + 
  # scale_fill_continuous(trans = 'log') + 
  coord_equal()

(sumAbundance = pr_grid_sc_elev@data$MeanAbundance %>% sum)


# Predict Critical model --------------------------------------------------------------------------------------

# No sensible results yet.

