Repo created 2018-02-17.

The script is intended to verify whether counting scats using a modified Jolly-Seber model is feasible. 

In general, the steps will be to:

1. Simulate scat deposition with a fixed rate \theta ('recruitment')
2. Simulate scat removals (i.e. encounters + 'survival') with p(detect) ~ track in grid cell `s`
3. Model scat encounters as a modified spatial Jolly-Seber model.
3a. Need to implement a per-capita recruitment rate, as data augmentation implies a varying recruitment rate. 

