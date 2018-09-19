

# Population estimation from unstructured scat surveys

Author: Alec Wong

Compiled: 09-19-2018




## Statistical Model

This model is intended to be used with broad-scale scat searches, where it is inefficient or impossible to perform clearing of standing crops of scats, and where distance sampling is not easily performed, such as with detection dog searches. 

Three parameters are estimated - <a href="https://www.codecogs.com/eqnedit.php?latex=\lambda" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda" title="\lambda" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\theta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta" title="\theta" /></a>, and <a href="https://www.codecogs.com/eqnedit.php?latex=p" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p" title="p" /></a>. 

## Analysis

### NOTICE: YOU WILL NEED TO INSTALL THE VERSION OF JAGS IN MY FORK IN ORDER TO RUN THIS ANALYSIS.

Run this to install. `devtools::install_github(repo = 'awong234/jagsUI')`.

Included are files to carry out an analysis using data from a survey of the Adirondacks moose population, collected in 2016. The file to run this will be `analysis.R`; it loads the files `data_cleaned.Rdata` and `metadata.Rdata`. The `data` object contains the observation data `y`, the visit array `vis`, the days between visits `days`, the number of grid cells visited `nSites`, the number of primary occasions (including initial deposition) `maxT`, and the maximum number of grid cell replicates `maxV`. 

The `metadata` object contains information on the grid cells visited `visitedGridInfo`, the temporal information of transect visitation `roundVisits`, and the grid cells overlaid upon each transect `grids`. 

See the article referenced for details on these structures.

There are several models, each defined and elaborated upon in the text files within this directory, including:

* `model_cov_full.txt`
* `model_cov_full_tl_shared.txt`
* `model_cov_no_temp.txt`
* `model_cov_no_temp_no_dcov.txt`
* `model_cov_reduced_continuous`
* `model_cov_reduced_continuous_tl_shared.txt`
* `model_cov_reduced_crit.txt`
* `model_cov_reduced_crit_tl_shared.txt`
* `model_cov_reduced_Dcov_only.txt`
* `model_null.txt`

## Data cleaning

The order of operations for cleaning and formatting the data are the following:

1. `trackDataCleaning.R`
1. `dataFormatting.R`

**Not all of the files are available, so these scripts will NOT develop the final data products included.** Much of the files required to run these scripts are larger than the file size limit permits. The scripts are included for transparency reasons alone, and missing files will need to be obtained upon request.

## Simulation

Included are files to perform simulations -- `scatCounting.R` -- testing the viability of the method described herein. The simulation uses the null model; view `model_null.txt`.

In short, the output describes an ability to identify parameters of the model, provided 25/100 grid cells were replicated once; that is, visited twice.

### Parameter estimation

![](https://github.com/awong234/population_from_scat_model/blob/master/images/lamplot.png)
![](https://github.com/awong234/population_from_scat_model/blob/master/images/p00.png)
![](https://github.com/awong234/population_from_scat_model/blob/master/images/thetaplot.png)

Figure 1: Parameter estimates across a gradient of lambda -- initial scat deposition -- and the proportion of sites that were replicated. Probability of detection was held at 50%, and theta was held at 1, representing a mean deposition rate of 1 scat per visit per grid cell. These numbers were deemed reasonable given experience at sites of heavy moose density.

Note that estimation improves dramatically after approximately 25% of the sites are revisited. Following this result, we proceeded with the analysis.
