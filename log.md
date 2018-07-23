# 2018-07-23

## ArcMap Activity

* [x] More errors found in track points, but not in lines. Fixing.
* [x] Errors in scat location, some not referenced properly to nearest site. Some scat collections labeled as collected in 07A1 were actually collected in 07A3.
* Not fixable : Four scats located >100m from nearest track point, but not in the vicinity of another site. Acceptable losses.

## Data Cleaning COMPLETE

* Scats referenced to nearest track point, saved into `scatsData.Rdata` as the object `scatsReferenced`. 
    * Note that distMoved is the 'residual' distance between the original collection location and the moved location.

| Percentile|   DistMoved|
|----------:|-----------:|
|       0.00|   0.0000083|
|       0.10|   0.6086150|
|       0.20|   1.1026085|
|       0.30|   1.5471718|
|       0.40|   2.0019971|
|       0.50|   2.5559814|
|       0.60|   3.3894581|
|       0.70|   4.2767710|
|       0.80|   5.7106574|
|       0.90|   8.2930854|
|       0.91|   8.6496837|
|       0.92|   9.1293406|
|       0.93|   9.7858546|
|       0.94|  11.3028923|
|       0.95|  12.1658219|
|       0.96|  13.1328228|
|       0.97|  17.1034756|
|       0.98|  23.2711577|
|       0.99|  46.2693123|
|       1.00| 923.7595729|

# 2018-07-22

The track lines layer as managed did not retain time information, which is critical for aligning scats to tracks. The scats need to be aligned to the nearest track point in time. 

I make the assumption of aligned time in scat records and track records. This ought to be mostly kept due to the fact that time was not manually recorded, but automatically by the device taken. The strong assumption is, then, that the devices are aligned in time such that matching the scat to the nearest track point in time does **not** strongly deviate its position in space. This residual can be compared.

## ArcMap activity

Cleaning will be undertaken again in ArcMap due to above loss of time information in the lines layers. The activity on 2017-07-21 will be replicated. This means the lines feature **should** be deprecated, but will be kept in case.

* [x] `07A3.Clearing` Renamed to `07A2.Clearing`.
* [ ] Error found in tracks, cleaning now; `02B3.Clearing` ought to be `02A2.Clearing`. This was actually noted on 07-20.


# 2018-07-21

## Formatting data for use in script.

## Scat cleaning 2016

* [x] Validate that scats sites are all labeled as their nearest transect; not done for 2017 yet.
* [x] Eliminated ~30 scat records without a corresponding track.

## ArcMap activity

### Track cleaning 2016

* [x] Error in 07A1.Clearing, this was the day Justin cleared a whole transect adjacent to the correct one. Correct.
* [x] One feature in 07A2 is a multi-track line. Split and name correctly.

# 2018-07-20

Viewing 2017 data, most scats are in the vicinity of dog tracks. Some are far away, certainly outside the 50m box and won't be counted in a 'visited' cell. Task is to have all the scat locations inherit the location of the **nearest** point on the corresponding dog track line.

## Track cleaning Issues 2017

* [x] 08B3 and 08B4 still mislabeled. Looks like two sites each that weren't picked up by the GIS.
    * Maybe renaming was unnecessary. Testing without rename.
* [x] One site 12A4 mislabeled as 12A6. 
    * Part of 12A6 (one of the multitracks exports) belongs to 12A4. This is from `12A6_07.10.17_JL.gpx`. 
    * Rename to `12A4_07.10.17_JL.gpx`.
* [x] Part of `12A6_07.09.17_SM.gpx` should be 12A4. Need to separate out.

## ArcMap activity

* A feature that consisted of 12A6 and part of 12A4 on July 9 2017 was split into two. 
    * The portion within 12A6 was renamed to match; 12A4.Sample1.
    * The portion within 12A4 was correctly identified and not altered.
* Renamed '12A6.Sample1' to '12A4.Sample1', which is the single-line portion of 12A4 sampled on 
* Made a number of edits to layer to fix GPS error, and points obviously out of search area. Original kept. 
* Renamed '08B4.Sample4' and '08B4.Sample5' to 08B3, each.
* Renamed '08B3.Sample4' and '08B3.Sample5' to 08B4, each.
* **Exported cleaned data, copy exists in Analysis.gdb**

## 2016 tracks

Justin had entered his date as a slightly different format than the rest. All fixed now. 

Now need to export to shapefile to examine in ArcMap.

### Issues

* One layer of 02A2 looks out of place
* One layer of 02B3 looks out of place, potentially swapped with the previous.

## 2016 scat data

# 2018-07-19

GPS errors in 2017 data fixed, some of the export errors made by JL fixed. Next step is to get 2016 data together.

## Compile all of the gpx tracks from 2017 into a single shapefile. 

Workflow: 

1. Fix inconsistencies in naming, split out files exported incorrectly.
2. Assemble point data together in R with the same metadata
3. Trim off GPS errors in ArcMap
4. Re-export for final use

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