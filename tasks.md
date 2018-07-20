# To do

## Track Data

* [x] Finish cleaning 2017
    * [x] Separate out multi-track files
        * Edits done in ArcMap, see log
* [ ] Clean 2016
* [ ] Format for analysis
    * [ ] Get interval between visits as data. **IMPORTANT**, as we're adjusting theta to be site-specific *daily deposition rate*
    * [ ] Get index of transects as data, for transect-level effect. Not sure the utility for prediction, though, since other parts of ADK are unreferenced to transect

## Scat Data

* [ ] Make scats inherit location from nearest point in corresponding track. 

## Covariate Data

* [ ] Format spatial habitat covariate. 
    * Reduce resolution to 50m, pull out relevant habitat types. Do this during runtime for null model, and model with transect effect

## Analysis

* [ ] Make model work according to commentary in `model.txt` file
* [ ] Format data. Need:
    * [ ] Index of transect ID as data
    * [ ] Interval between visits as data
* [ ] Summarize model output

## Writing

* [ ] Import intro from deprecated Ch 1. 
* [ ] Methods
* [ ] Results
* [ ] Discussion
    * [ ] Don't forget commentary about model extensions