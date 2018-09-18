

<style>
img.parimg {
  width: 300px;
}
div.caption { 
  padding: 25px 100px 25px 100px;
}
</style>


# Population estimation from unstructured scat surveys

Author: Alec Wong
Compiled: 09-18-2018




## Simulation

Included are files to perform simulations -- `scatCounting.R` -- testing the viability of the method described herein. The simulation uses the null model; view `model_null.txt`.

In short, the output describes an ability to identify parameters of the model, provided 25/100 grid cells were replicated once; that is, visited twice.

### Parameter estimation

<div>
<img class='parimg' src = 'images/old/lamplot.png'>
<img class='parimg' src = 'images/old/thetaplot.png'>
<img class='parimg' src = 'images/old/p00.png'>
</div>

<div class='caption'>
<center>
Figure 1: Parameter estimates across a gradient of lambda -- initial scat deposition -- probability of detection, and the proportion of sites that were replicated.
</center>
</div>



Note that estimation improves dramatically after approximately 25% of the sites are revisited.
