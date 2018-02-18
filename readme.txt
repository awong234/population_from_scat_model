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