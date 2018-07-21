# To do

## Track Data

* [x] Clean 2017
    * [x] Separate out multi-track files
        * Edits done in ArcMap, see log
    * [x] Clean scat data
        * [x] Check to see that scat sites aren't mislabeled
* [ ] Clean 2016
    * [x] Clean transect data
        * [x] Separate out multi-export tracks
    * [x] Clean scat data
        * [x] Some scats found to be mislabeled. Fix. 
* [ ] Format for analysis
    * [ ] Get grid overlay over tracks, count up visits & times of visit if possible. 
        *   **MAY NOT BE POSSIBLE FOR RECONSTRUCTED VISITS** be careful.
    * [ ] Get interval between visits as data. **IMPORTANT**, as we're adjusting theta to be site-specific *daily deposition rate*
    * [ ] Get index of transects as data, for transect-level effect. Not sure the utility for prediction, though, since other parts of ADK are unreferenced to transect

## Scat Data

* [ ] Make scats inherit location from nearest point in corresponding track. 
* [ ] Reference scats to grid cells, indicate whether first, second, third, etc. replicate.

## Covariate Data

* [ ] Format spatial habitat covariate. 
    * Reduce resolution to 50m, pull out relevant habitat types. Do this during runtime for null model, and model with transect effect

## Grid Overlay

* [ ] For spatial prediction, will need a 50m x 50m grid over all prediction areas. This means masking out the towns and water bodies. Do this in ArcMap. 

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