# To do

## Track Data

* [ ] Clean 2017
    * [x] Separate out multi-track files
        * Edits done in ArcMap, see log
    * [ ] Clean scat data
        * [x] Check to see that scat sites aren't mislabeled
* [x] Clean 2016
    * [x] Clean transect data
        * [x] Separate out multi-export tracks
    * [x] Clean scat data
        * [x] Some scats found to be mislabeled. Fix. 
* [x] Format for analysis
    * [x] Get grid overlay over tracks, count up visits & times of visit if possible. 
    * [x] Get interval between visits as data. **IMPORTANT**, as we're adjusting theta to be site-specific *daily deposition rate*
    * [x] Get index of transects as data, for transect-level effect. Not sure the utility for prediction, though, since other parts of ADK are unreferenced to transect

## Scat Data

* [x] Make scats inherit location from nearest point in corresponding track. 
* [x] Reference scats to grid cells, indicate whether first, second, third, etc. replicate.

## Covariate Data

* [ ] Format spatial habitat covariate. 
    * Reduce resolution to 50m, pull out relevant habitat types. Do this during runtime for null model, and model with transect effect
* [ ] Format distance to road covariate
    * Options include:
        * Distance to highway and distance to minor road as separate covariates
            * Extra parameters to estimate, potentially covarying
        * Distance to any road
            * May weaken the signal if response is different between road classes
        * Distance to any road, weighted by highway class
            * How to do this? 
* [ ] Format distance to water body covariate
    * Should **definitely** weight by water body type or separate. Lakes are different from rivers are different from ephemeral pools are different from wetlands.
        * Start with lakes, wetlands, due to cooling + forage. Moose & rivers not known to be associated.
* [ ] Format distance to human settlement
    * Weight by population size? Potentially - 
        * $(max(\text{Dist}, .01) * \text{Population})$ = Weighted population distance
        * Closer distance only matters if population is high. 
        * One method is to calculate for all towns, and take the maximum value.
* [ ] Format temperature maximum prediction
* [ ] Format elevation (easy)
* [ ] Potential covariate: mobile data availability, from the NYS Office of Information Technology Services- GIS Program Office. Latest data is 2014, so fairly applicable.
    * Rationale is that 4G will be available in all areas where people demand it (*i.e.* tend to spend more time), 2G/3G where people spend less time, and no signal where people spend no time. It is thus a proxy of time spenty as well as population size.

## Grid Overlay

* [ ] For spatial prediction, will need a 50m x 50m grid over all prediction areas. This means masking out the towns and water bodies. Do this in ArcMap. 

## Analysis

* [ ] Make model work according to commentary in `model.txt` file
* [x] Format data. Need:
    * [x] Index of transect ID as data
    * [x] Interval between visits as data
* [ ] Summarize model output

## Writing

* [ ] Revise intro based on Angela's comments
* [x] Methods
* [ ] Results
* [ ] Discussion
    * [ ] Don't forget commentary about model extensions