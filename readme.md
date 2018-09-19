

# Population estimation from unstructured scat surveys

Author: Alec Wong

Compiled: 09-19-2018




## Statistical Model

This model is intended to be used with broad-scale scat searches, where it is inefficient or impossible to perform clearing of standing crops of scats, and where distance sampling is not easily performed, such as with detection dog searches. 

Three parameters are estimated - 
<a href="https://www.codecogs.com/eqnedit.php?latex=\lambda" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda" title="\lambda" /></a>,
<a href="https://www.codecogs.com/eqnedit.php?latex=\theta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta" title="\theta" /></a>, and
<a href="https://www.codecogs.com/eqnedit.php?latex=p" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p" title="p" /></a>.

For the moment, consider these parameters with respect to a single grid cell among many visited over the course of a survey.

<a href="https://www.codecogs.com/eqnedit.php?latex=\lambda" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda" title="\lambda" /></a> 
describes the mean initial deposition of scats before any sampling is performed.

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{gathered}&space;\Delta_0&space;\sim&space;\text{Poisson}(\lambda)\\&space;\lambda&space;=&space;\exp(\boldsymbol{\beta}*\textbf{X})&space;\end{gathered}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{gathered}&space;\Delta_0&space;\sim&space;\text{Poisson}(\lambda)\\&space;\lambda&space;=&space;\exp(\boldsymbol{\beta}*\textbf{X})&space;\end{gathered}" title="\begin{gathered} \Delta_0 \sim \text{Poisson}(\lambda)\\ \lambda = \exp(\boldsymbol{\beta}*\textbf{X}) \end{gathered}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\theta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta" title="\theta" /></a>
represents the mean daily accumulation rate of scats, such that:

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{gathered}&space;\Delta_t&space;\sim&space;\text{Poisson}(\theta*d_t)\\&space;\theta&space;=&space;\exp(\boldsymbol{\beta}*\textbf{X})&space;\end{gathered}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{gathered}&space;\Delta_t&space;\sim&space;\text{Poisson}(\theta*d_t)\\&space;\theta&space;=&space;\exp(\boldsymbol{\beta}*\textbf{X})&space;\end{gathered}" title="\begin{gathered} \Delta_t \sim \text{Poisson}(\theta*d_t)\\ \theta = \exp(\boldsymbol{\beta}*\textbf{X}) \end{gathered}" /></a>

where, <a href="https://www.codecogs.com/eqnedit.php?latex=d_t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?d_t" title="d_t" /></a> is the intervening days between visits on occasions t and t+1. 

Finally, p is derived from visitation of a grid cell multiple times within a single occasion t. That is to say, if the grid cell is observed more than once, we can construct a probability statement conditional on the latent variables 
<a href="https://www.codecogs.com/eqnedit.php?latex=\lambda" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda" title="\lambda" /></a> and 
<a href="https://www.codecogs.com/eqnedit.php?latex=\theta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta" title="\theta" /></a>. 

<a href="https://www.codecogs.com/eqnedit.php?latex=y_{t,r}&space;\sim&space;\text{Binomial}(N_{t,r},p)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{t,r}&space;\sim&space;\text{Binomial}(N_{t,r},p)" title="y_{t,r} \sim \text{Binomial}(N_{t,r},p)" /></a>

Here, the observations 
<a href="https://www.codecogs.com/eqnedit.php?latex=y_{t,r}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{t,r}" title="y_{t,r}" /></a> 
are indexed by occasion t and replicate observation r. The observations are conditional on the population of scats available to be sampled 
<a href="https://www.codecogs.com/eqnedit.php?latex=N_{t,r}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N_{t,r}" title="N_{t,r}" /></a>

The population of scats has the following model:

<a href="https://www.codecogs.com/eqnedit.php?latex=N_{t,r|t=1}&space;=&space;\Delta_{0}\\&space;N_{t,r|t>1,r=1}&space;=&space;N_{t-1,r_{max}}&space;-&space;y_{t-1,r_{max}}&space;&plus;&space;\Delta_{t-1}\\&space;N_{t,r|t>1,r>1}&space;=&space;N_{t,r-1}&space;-&space;y_{t,r-1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N_{t,r|t=1}&space;=&space;\Delta_{0}\\&space;N_{t,r|t>1,r=1}&space;=&space;N_{t-1,r_{max}}&space;-&space;y_{t-1,r_{max}}&space;&plus;&space;\Delta_{t-1}\\&space;N_{t,r|t>1,r>1}&space;=&space;N_{t,r-1}&space;-&space;y_{t,r-1}" title="N_{t,r|t=1} = \Delta_{0}\\ N_{t,r|t>1,r=1} = N_{t-1,r_{max}} - y_{t-1,r_{max}} + \Delta_{t-1}\\ N_{t,r|t>1,r>1} = N_{t,r-1} - y_{t,r-1}" /></a>

where 
<a href="https://www.codecogs.com/eqnedit.php?latex=N_{t,r}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N_{t,r}" title="N_{t,r}" /></a>
is deterministically reduced by observations 
<a href="https://www.codecogs.com/eqnedit.php?latex=y_{t,r}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{t,r}" title="y_{t,r}" /></a> 
and increased by accumulation
<a href="https://www.codecogs.com/eqnedit.php?latex=\Delta_t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Delta_t" title="\Delta_t" /></a>.

This model is applied to all grid cells observed, assuming independence between them.


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
