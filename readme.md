# Scat Counting Simulation
Alec Wong  
February 19, 2018  




# Overview

The major components of this simulation are as follows:

* Simulate a population of scats with a rate of deposition, to estimate the population based on simulated collections.
* Simulate collection of scats by probabilistically drawing Bernoulli samples based upon dog tracks within the grid cell coincident with the scat.
* Estimate scat deposition rate per visit using a modified spatial Jolly-Seber model that incorporates scat 'recruitment' and 'survival', where 'recruitment' is a rate of deposition, and 'survival' is fixed to 0, since all scats encountered are immediately removed.

# Simulating the data

## Dog tracks

I obtained two example dog tracks from 2017 to examine and test under this framework. The sites referenced are 12B2, and 15A4, during three consecutive visits within mid-July to August. 

The site 12B2 exhibited extremely high moose density, and thus the dog movement is no doubt affected heavily by this. In the case that this movement -- being unexplained in the model we're using -- affects model performance, the other site 15A4 was selected to be tested as well. This site exhibited no moose collections throughout the summer, and the dog movement is taken to be the most 'natural' search pattern under normal conditions. 15A4 is also relatively easy terrain, so the dogs' movement was not inhibited by thorns, slope, water, or the like.

### Examples



![](readme_files/dogTrack.gif)

In the figure above, we see two triangular transects, and notice that multiple 'Rounds' or visits are made to each site, depicted in the image. They are centered such that the centroid of the transects are on 0, but independently scaled so that the transects are the same size relative to each other. The grid was first generated to be 50m x 50m, and then scaled by the same amount. To easily verify the grid size, count the squares covering the vertical portion of 12B2; they number 20, indicating that the length of the side is approximately 1km, which is true.

I scaled and centered the dog tracks so that they may sample the same population of scats, to observe any differences in estimates. Notice that 12B2 is more 'wiggly', owing to the frequency of moose scat encounters.

# Simulation of scat deposition and encounters

The deposition of scats is to be done in a Poisson random fashion, with additional 'recruitment' being added with some rate $\theta$. 

Validation requires knowledge of where the scats were generated, where dog tracks intersect the scats' grid cells, and whether they are being removed properly. 

## Scat simulation

This part is easy enough. I test using an initial deposition of 500 uniformly distributed scat piles. 

![](readme_files/figure-html/scatPlot-1.png)<!-- -->

In the above plot, we have scat locations (red '+'), and grid cells with numbers marking how many scat piles exist within the grid. The function written to count any given point layer within grid cells is working properly - it is zoomed in to demonstrate this, but is correct for the wider grid on the whole. 

## Dog track points within grid

I also need to verify that I'm detecting any dog track within the grids. Below, I test a function generating a probability of detection based on track length within a grid, or based on some baseline probability of detection, but for validation I am using detection == 1.



![](readme_files/scatDeposPlot.gif)

In the plot above, we see 500 initial scats deposited, followed by a random number of scats deposited afterward modeled as Poisson with mean = 200. Since they are not reproducing in the normal sense, the recruitment rate is independent of the scat population size[^1]. 

[^1]: However, the recruitment rate might be related to the initial size deposited, since they both indicate more moose on the transect.

We observe a series of 500, 201, 198 scats deposited in the initial sample, after the first sample, and after the second sample, respectively. Of course, there is no sampling after the third sample, so that recruitment is not simulated.

In the following plot, pay particular attention to the highlighted areas:

![](readme_files/scatDeposPlot2.gif)

Those scats are changing their status from 'not removed', to 'removed', demonstrating the function's proper operation. They are deposited amongst the initial sample, and on the first visit the dog track does not intersect their grid cell, and so they are not removed in the first round. In the second round, notice that the dog track intersects their grid cell, and so they are removed[^2].

[^2]: The simulation has p(detect) = 1 if there is any track within the grid whatsoever. Later, we can adjust to make it a function of distance traveled, area covered, time, or a combination of these.

## Dataset obtained

The simulated dataset is obtained by filtering out only those scats that were removed (since we would not have information about those not removed), as below:



<table>
 <thead>
  <tr>
   <th style="text-align:center;"> ID </th>
   <th style="text-align:center;"> x </th>
   <th style="text-align:center;"> y </th>
   <th style="text-align:center;"> RoundDeposited </th>
   <th style="text-align:center;"> pEnc </th>
   <th style="text-align:center;"> Removed </th>
   <th style="text-align:center;"> RoundRemoved </th>
   <th style="text-align:center;"> gridID </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 2.07 </td>
   <td style="text-align:center;"> -0.39 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 317 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 16 </td>
   <td style="text-align:center;"> -0.01 </td>
   <td style="text-align:center;"> -1.24 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 160 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 1.29 </td>
   <td style="text-align:center;"> 0.55 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 487 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 37 </td>
   <td style="text-align:center;"> 1.37 </td>
   <td style="text-align:center;"> -0.52 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 284 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 41 </td>
   <td style="text-align:center;"> 1.49 </td>
   <td style="text-align:center;"> -0.41 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 314 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 43 </td>
   <td style="text-align:center;"> 1.31 </td>
   <td style="text-align:center;"> 0.61 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 487 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 45 </td>
   <td style="text-align:center;"> 0.14 </td>
   <td style="text-align:center;"> -1.20 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 161 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 50 </td>
   <td style="text-align:center;"> 0.90 </td>
   <td style="text-align:center;"> 0.93 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 542 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 57 </td>
   <td style="text-align:center;"> -0.85 </td>
   <td style="text-align:center;"> -1.61 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 68 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 60 </td>
   <td style="text-align:center;"> -0.43 </td>
   <td style="text-align:center;"> 1.31 </td>
   <td style="text-align:center;"> 0 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 621 </td>
  </tr>
</tbody>
</table>

<br />


Let's look at individual 41. The scat was deposited in round 0, but wasn't removed until round 2 - we shall see if this is correct.



![](readme_files/indScatGIF.gif)

In the above plot, notice that this particular scat pile is not encountered in the first occasion, due to the lack of any track points in its grid cell. In the second occasion, it is encountered due to the change in the dog track pattern directing it into the grid cell of the scat.

## Incorporating detection probability

### Indicator model

Simulation of detection probability in the 'indicator' fashion proceeds as follows:

$$
p(\text{detect}) = p_0 * I(\text{track in grid})
$$

where $I$ is the indicator function evaluating to 1 if there exists any track within the grid, and 0 if there does not. 

A Bernoulli trial is applied to each scat, simulating its encounter (and removal). Since dogs appear to be extremely discerning when it comes to detecting scats, I test using a detection probability of 0.8, but any value can be simulated.

### Length model

In future revisions, a scat's detection probability will be determined by the length of track in the grid cell that the scat occupies, in this fashion:

$$
p(\text{detect}) = p_0 * \sum_{i=1}^{I-1} \text{dist}(x_{i+1},x_i)
$$

where, $\text{dist}(x_i,x_j)$ is the euclidean distance between points $x_i$ and $x_j$. If there exists no dog track points in the grid cell, the probability of detection is zero. If by chance there is a singular point represented in a grid cell, instead of assigning it an arbitrary distance I will consider it with distance 0, and as such equivalent to no track length within the grid.

## Observation of detection process

### Indicator process

I identify a series of scats that were deposited in the first round and removed in the last for observation. 


  ID       x       y  RoundDeposited    pEnc  Removed   RoundRemoved    gridID
----  ------  ------  ---------------  -----  --------  -------------  -------
  66   -1.12    0.78  0                  0.8  1         3                  501
  99    1.44   -0.48  0                  0.8  1         3                  285
 157    0.14   -1.22  0                  0.8  1         3                  161
 346   -1.07   -0.26  0                  0.8  1         3                  328
 419    0.91    0.81  0                  0.8  1         3                  542


#### Example: Individual 66

Observe individual 66:


                     ID           x           y  RoundDeposited    pEnc  Removed   RoundRemoved    gridID
------------------  ---  ----------  ----------  ---------------  -----  --------  -------------  -------
Round 0  Snapshot    66   -1.124147   0.7804769  0                  0.0  0         NA                 501
Round 1  Snapshot    66   -1.124147   0.7804769  0                  0.0  0         NA                 501
Round 2  Snapshot    66   -1.124147   0.7804769  0                  0.0  0         NA                 501
Round 3  Snapshot    66   -1.124147   0.7804769  0                  0.8  1         3                  501

Evidently, the probability of encounter was 0 until the final occasion, meaning that we ought to observe no track until round 3. 



![](readme_files/ind66.gif)

In the above plot, we do indeed see no track until round 3. 

Another way to look at the data is to see whether -- in any round -- given track in grid, about 80% of those scats should be removed. Let's observe:


```r
# Of those with p(encounter) > 0, how many removed?
table1 = scatSim$ScatRecords$`Round 1` %>% filter(pEnc > 0) %>% {summary(.$Removed)} 
names(table1) = c("Not Removed", "Removed")
print(table1)
```

```
## Not Removed     Removed 
##          15          50
```

We see that 50 individuals of 65 individuals are removed in Round 1: this is approximately 80%, or more specifically, 0.7692308 %.

# Analysis of data

The data obtained are a record of collections per grid, per visit. Below, I randomly sample a few records from each visit to show the counts.


```
## # A tibble: 15 x 3
## # Groups:   RoundRemoved [3]
##    RoundRemoved gridID     n
##          <fctr> <fctr> <int>
##  1            1    328     1
##  2            1    222     1
##  3            1    621     2
##  4            1    159     1
##  5            1    317     1
##  6            2    346     1
##  7            2    487     1
##  8            2    566     1
##  9            2     67     1
## 10            2    502     1
## 11            3    569     2
## 12            3    299     2
## 13            3    568     1
## 14            3    501     2
## 15            3     96     1
```

# Final update notes

At this time, the simulation is complete, and the dataset is obtained. Next tasks are to incorporate detection probability into the simulation, and provide a model with which JAGS can estimate the appropriate parameters, with the ultimate goal of a population estimate of scats. 

Options for simulating detection probability:

* Total length of path in grid cell
* Total time spent in grid cell
* Total length of path in grid cell times total time spent in grid cell (my preferred option)
* Fraction of grid cell area covered, assuming a fixed detection radius around each GPS fix point.

### Footnotes
