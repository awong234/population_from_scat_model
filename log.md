# 2018-07-20

Viewing 2017 data, most scats are in the vicinity of dog tracks. Some are far away, certainly outside the 50m box and won't be counted in a 'visited' cell. Task is to have all the scat locations inherit the location of the **nearest** point on the corresponding dog track line.

# 2018-07-19

GPS errors in 2017 data fixed, some of the export errors made by JL fixed. Next step is to get 2016 data together.

## Compile all of the gpx tracks from 2017 into a single shapefile. 

Workflow: 

1. Fix inconsistencies in naming, split out files exported incorrectly.
1. Assemble point data together in R with the same metadata
2. Trim off GPS errors in ArcMap
3. Re-export for final use

* Work on cleaning 2016 tracks, and do the same.

# 2018-07-18

## MISLABELS

* 06C1-2017-08-27 should be 06C2. 
* 06C2-2017-08-27 should be 06C1.
* 08B3-2017-08-15 should be 08B4.
* 08B4-2017-08-15 should be 08B3.
* 10B5-2017-08-11 should be 10B3.
* 12A6-2017-07-09 should be 12A4.

Plan is to export all of the .gpx tracks to ESRI shapefile, perform spatial join with a distance of ~ 100m, then find which ones are in error. Fix names.

Now fixed, explicit in R code.

## GPS error

There are some GPS errors. Will fix manually where they occur and re-export.

# 2018-03-02

need to move bbox into simScats.

Repo created 2018-02-17.

The script is intended to verify whether counting scats using a modified Jolly-Seber model is feasible. 

In general, the steps will be to:

1. Simulate scat deposition with a fixed rate \theta ('recruitment')
2. Simulate scat removals (i.e. encounters + 'survival') with p(detect) ~ track in grid cell `s`
3. Model scat encounters as a modified spatial Jolly-Seber model.
3a. Need to implement a per-capita recruitment rate, as data augmentation implies a varying recruitment rate. 

---------

Interesting note:

Time between fixes is almost always under a minute. The only times it exceeds a minute are in 12B2 round 1. Most fixes are exactly 1s apart. 

An index for search effort will need to include not only distance traveled within a block, but also time spent traveling it; one might imagine that for a given distance, moving quicker results in lower detection than moving slowly. 

I would need to - for each grid cell - tabulate the sum of the distances traveled from one fix to another, and then the sum of the time taken to travel it, and multiply them for an index of effort.