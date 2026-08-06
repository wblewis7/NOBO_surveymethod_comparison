# NOBO_surveymethod_comparison
Formatted data and R/NIMBLE code for performing Bayesian hierarchical analyses estimating fall density of northern bobwhite (Colinus virginianus) in central Florida, USA, via point counts, multi-observer SCR quadrat surveys, and trapping SCR
---
Authors: Robin L. Schmitt, William B. Lewis, Amanda S. Cramer, Geoff Beane, and James A. Martin
<br />
Manuscript Title: A comparison of northern bobwhite density estimators
<br />
Contact: William B. Lewis, University of Georgia, william.lewis27@uga.edu
---

---

# Metadata
Northern bobwhite abundance was surveyed during the falls (October - November) of 2021/2022 at two study properties in central Florida, USA. Properties were divided into 6 total survey courses (4 at Property 1, 2 at Property 2).
Data were surveyed at 19 500x500m sampling units (2 - 4 sampling units per survey course) via three different surveys methods. First, covey abundance was surveyed from the center of each sampling unit via 500-m fixed-radius, distance-sampling
point counts. Second, covey abundance was surveyed at the sampling units via the multi-observer quadrat method. Four observers stood at the mid-points of the edges of the sampling unit and recorded the location and time of 
detection of each unique covey detected. Observers reconciled observations after surveys, identifying the number of unique coveys detected and which coveys were detected by which observer. Both point counts and quadrat surveys occurred
from 45-minutes prior to sunrise until sunrise and were repeated twice per year. Lastly, individual abundance was surveyed via trapping spatial capture-recapture (SCR). We deployed baited walk-in traps in a 4x4 array, spaced 150m apart, at each sampling unit. Each trap
was run for ~ 7.5 (1 - 9) occasions. Upon capture, birds were banded with a uniquely-numbered aluminum leg band and released. Supplemental feeding at the properties was generally restricted ~10 days prior to trapping, though two courses had feed withheld
~20 days prior to trapping. Each of the three datasets was analyzed via a Bayesian framework via NIMBLE in R. 

<br>
<br>
<br>

----

# Point Count
We model point count data via the hierarchical distance-sampling models of Royle et al. (2004, Ecology) with modifications to account for repeat surveys (Lewis et al. 2025. Ecological Solutions and Evidence). We model expected abundance via effects of survey course and year, and model the distance of detections based on a half-normal detection function. Not that this model estimates covey abundance at sampling units, which we convert to individual abundance via an average covey size of 12 birds/covey.
<br>
## NOBO_CFL_modeldata_pointcount.gzip
Formatted data for estimating abundance from the point count dataset are provided in "NOBO_CFL_modeldata_pointcount.gzip".
### Npoints
The number of survey points across properties (19).
### Nyrs
The total number of survey years (2).
### Nvisits
The number of repeat surveys performed at each survey point in each year (2).
### y 
A 3-D array giving the total number of coveys detected at each point (1st dimension), repeat surveys (2nd dimension), and year (3rd dimension).
### ydb
A 4-D array giving the number of covey detections at each point (1st dimension), in each distance bin (2nd dimension), repeat surveys (3rd dimension), and year (4th dimension).
### nB
The total number of distance bins for distance sampling (5).
### db
The boundaries of the distance bins (m).
### pix
The proportional area of each distance bin.
### point.area
The area (square meters) surveyed by point counts.
### ncourse
The total number of survey courses (6).
### course
A vector giving the course number for each point in Npoints.

<br />

## Analysis_NOBO_CFL_pointcount.R
R and NIMBLE code for running the Bayesian hierarchical point count model.


<br>
<br>

----

# Quadrat SCR
We estimate density from reconciled quadrat surveys using a hierarchical SCR model (Royle & Young. 2008. Ecology) with modifications to account for multiple observers. We take a data augmentation approach, modeling whether or not each augmented covey at each point in each year is a true covey (z) based on the inclusion parameter (psi). We calculate expected density (lambda) via course and year effects, then convert to psi by dividing lambda by the total number of possible individuals (M). We also model the spatial locations of each covey bounded within the limits of the 25ha study area. We decompose the detection process into availability (probability of a covey calling during a survey) and conditional detection probability (probability of observers detecting a calling covey). For each survey visit, we model availability of coveys based on z and property-specific probability of availability. Note that this formulation fixes abundance across repeat survey visits at the same sampling unit in each year, but allows for variable availability and spatial locations across visits. For each of the 4 observers, we then calculate detection probability for each covey based on a half-normal detection function and the distance between the covey and the observer. We then model detections for each potential covey by each observer, for each visit in each sampling point in each year, based on a Bernoulli distribution. As with point counts, this model estimates covey abundance and can be converted to individual abundance via an average covey size.
<br>
## NOBO_CFL_modeldata_quadratSCR.gzip
Formatted data for estimating abundance from the point count dataset are provided in "NOBO_CFL_modeldata_quadratSCR.gzip".
### Nproperty 
The total number of study properties (2).
### Ncourse 
The total number of survey courses across the 2 properties (6).
### Npoints 
The total number of survey quadrats (19).
### Nyear
The total number of survey years (2).
### course
A vector giving the course number for each survey point in Npoints.
### Nvisits
The number of repeat surveys performed at each survey point in each year (2).
### property
A vector giving the property for each study course in Npoints.
### M
The maximum number of potential coveys for data augmentation (100).
### Xlim
A matrix giving the lower/upper (columns) x-coordinates of the 25-ha survey quadrats at each survey point in Npoints (rows).
### Ylim
A matrix giving the lower/upper (columns) y-coordinates of the 25-ha survey quadrats at each survey point in Npoints (rows).
### J
The number of observers for each quadrat survey (4).
### X
A 3-D array giving the x (1st element of 3rd dimension) and y (2nd element of 3rd dimension) of each observer in J (1st dimension) during quadrat surveys at each survey point in Npoints (2nd dimension).
### z
An M x Npoints x Nyear array fixing the z-state to 1 for detected coveys during data augmentation, with NA for z-states to be estimated
### y.Q
A 5-D array of observed covey detections denoting whether or not each covey in M (1st dimension) was detected at each survey point in Npoints (2nd detection), by each observer in J (3rd dimension), for each year in nyrs (4th dimension) and each repeat survey in Nvisits (5th dimension).
### constraint.lambda
A Npoints x Nyear matrix of 1s for maximizing the contraint likelihood to ensure that expected abundance (lambda) is not greater than M.

<br />

## Analysis_NOBO_CFL_quadratSCR.R
R and NIMBLE code for running the Bayesian hierarchical quadrat SCR model.


<br>
<br>

----

# Trapping SCR
We analyze trapping data using a SCR model with data augmentation. As with quadrats, we calculate lambda via course and year effects and convert to psi by dividing by M. This model estimates activity centers of captured and augmented individuals, and estimates detection probability at each trap based on the distance between the trap and the activity center. We used the 'nimbleSCR' package (Bischof et al. 2020. Wildlife Biology) to improve computational efficiency. Notably, we grid the trapping area into 50m cells and constrain activity centers to be within upland areas based on a mask of wetland areas developed based on FWC/FNAI 2023. Additionally, these functions only calculate detection probability for traps which could feasibly catch an individual. The furthest distance within a year between captures in our data was 422m, so we model traps >844m from activity centers as being unable to capture individuals. These functions also scale the x/y axes such that the origin is in the bottom-left at 0,0.
<br />
We provide formatted point count, quadrat, and trapping data, as well as R and NIMBLE code to analyze the data.

