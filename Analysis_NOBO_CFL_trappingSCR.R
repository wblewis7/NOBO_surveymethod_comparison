# R and nimble code for estimating fall covey abundance of northern bobwhite (Colinus
#   virginianus) via trapping-based spatial capture-recapture in Central Florida, USA 
#   2021-2022.
# From Schmitt, R. L., A. S. Cramer, G. Beane, W. B. Lewis, and J. A. Martin. A 
#   comparison of northern bobwhite density estimators.
# Contact: Will Lewis, University of Georgia, william.lewis27@uga.edu



# Fall bobwhite coveys were trapped in the falls of 2021/2022 in Central Florida, USA.
#   trappings occurred at 2 properties, each of which contained 2 - 4 trapping courses.
# We trapped bobwhite in early-mid October in each year using baited walk-in funnel
#   traps. We deployed 32-64 traps at each course, approximately 150 m apart. Each course
#   was trapped for up to 9 occasions. Traps were checked twice daily.
# Models could not be run for property 2 in the second year (2022) because there were
#   no spatial recaptures.
# We analyze capture data using a Bayesian spatial capture-recapture model with data
#   augmentation. This model estimates course-level abundance by estimating activity centers
#   of captured and augmented individuals, and estimates detection probability at each trap
#   based on the distance between the trap and the activity center.
# We used the 'nimbleSCR' package (Bischof et al. 2020. Wildlife Biology) to improve computational
#   efficiency. Notably, we grid the trapping area into 50x50m cells and constrain activity centers to be 
#   within upland areas based on a mask of wetland areas based on FWC/FNAI 2023. 
#   Additionally, these functions only calculate detection probability for traps which could
#   feasibly catch an individual. The furthest distance within a year between captures in our
#   data was 422m, so we model traps >844m from activity centers as being unable to capture
#   individuals. These functions also scale the x/y axes such that the origin is in the
#   bottom-left at 0,0.





library(nimbleSCR)


# Loading in data
load("NOBO_CFL_modeldata_trappingSCR.gzip")
str(data.trappingSCR)

# Nproperty is the total number of study properties (2).
# Ncourse is a vector of length Nproperty giving the total number of trapping courses at each property.
# Nyear is a vector of length Nproperty giving the number of trapping years at each property.
# M is a 3-D array giving the potential number of birds within the population at
#   at each course (2nd dimension) at each property (3rd dimension) in each year
#   (1st dimension). Includes captured and augmented individuals.
# x.max is a 3-D array giving the maximum x coordinates (in 50x50m grid cells) of the trapping area at each 
#   course (2nd dimension) at each property (3rd dimension) in each year
#   (1st dimension).
# y.max is a 3-D array giving the maximum y coordinates (in 50x50m grid cells) of the trapping area at each 
#   course (2nd dimension) at each property (3rd dimension) in each year
#   (1st dimension).
# ones is a 4-D array specifying values of 1 to maximize the likelihood of the 
#   'dHabitatMask' function, which constrains activity center locations to be
#   within upland based on the wetland habitat mask (habitat.mx). Dimensions correspond to maximum
#   M across courses/years (1st), courses (2nd), years (3rd), and properties (4th).
#   Values of NA correspond to potential individuals above the value of M for the 
#   given course/year/property and don't affect the likelihood.
# lowerCoords gives he x and y coordinates (units are grid cells) for the origin of plots of
#   trapping area (0,0 in bottom left). This is the same across all courses/years.
# upperCoords is a 4-D array giving the upper x and y (1st and 2nd, respectively,
#   elements of 2nd dimension) for the trapping areas (units are grid cells) at each course (1st dimension),
#   in each year (3rd dimension), and property (4th dimension).
# habitat.mx is a 5-D array denoting whether each 50x50m cell in the gridded trapping 
#   area could be an activity center. Values of 1 denote grid cells which are upland
#   and within 500m of traps. Dimensions correspond to maximum y-coordinates from y.max
#   (1st dimension) across courses/years, maximum x-coordinates from x.max (2nd dimension)
#   across courses/years, courses (3rd dimension), years (4th dimension), and properties
#   (5th dimension). Values of NA correspond to grids outside of the trapping area for the 
#   given course/year/property and don't affect the likelihood.
# z is a 4-D array giving the known state (real/fake) for each individual in M (1st dimension)
#   at each course (2nd dimension), for each year (3rd dimension), and property (4th dimension).
#   Values of 1 denote individuals which were captured (known alive) while z-state will be
#   estimated for augmented birds (NA values).
# nMaxDetectors is a 3-D array giving the maximum number of traps that an individual was captured in
#   at each course (2nd dimension) at each property (3rd dimension) and year (1st dimension).
# y is a 5-D array giving the total number of captures of each individual in M (1st dimension)
#   in different traps (2nd dimension) at each course (3rd dimension), in each year (3th dimension),
#   and property (5th dimension). Note that columns are organized by the number of unique traps
#   and individual was captured in (1st, 2nd, etc.) rather than the trap number. Values of -1 indicate
#   no captures, with augmented birds being all -1. These are the data for the 'dbinomLocal_formal' function.
# nbDetections is a 4-D array giving the total number of traps each individual in M (1st dimension) was 
#   captured in for each course (2nd dimension), year (3rd dimension), and property (4th dimension).
#   Values of NA correspond to values above M for the given course/year/property and don't affect the likelihood.
# yDets is a 5-D array giving the trap ID in which each individual in M (1st dimension) was captured.
#   2nd dimension is the unique traps (max 3), 3rd dimension is course, 4th dimension is year, and
#   5th dimension is property. Values of -1 indicate no captures, with augmented birds being all -1.
#   This is the same as y, but gives the trap ID of captures rather than the number of captures per trap.
# n.detectors is a 3-D array giving the total number of traps run at each course (2nd dimension),
#   at each property (3rd dimension), in each year (1st dimension).
# trials is a 4-D array giving the total number of trapping occasions for each trap. Dimension correspond
#   to the maximum number of traps across courses/years (1st dimension), courses (2nd dimension), years
#   (3rd dimension), and properties (4th dimension).
# sigindex is a 3-D array indexing the detection process at courses (2nd dimension), years (1st dimension),
#   and properties (3rd dimension). Supplemental food was generally restricted ~10 days before trapping, but
#   this was increased to ~20 days before trapping for Courses 1 and 3 at the first property. Modeling separate
#   effects on the detection process (g0 and sigma) based on whether or not food was restricted 20 (1) or 10 (2)
#   days before trapping.
# detector.xy is a 5-D array giving the x and y coordinates (1st and 2nd elements, respectively, of 2nd dimension)
#   in 50x50m grid cells of traps (1st dimension). 3rd dimension is course, 4th dimension is year, and 5th dimension
#   is property. NA values correspond to traps beyond the number of traps (n.detectors) at a given course/year.
# n.cells is a 3-D array giving the total number of 50x50m grid cells in which activity centers may be located for
#   each course (2nd dimension) at each property (3rd dimension) and year (1st dimension). Corresponds to the sum
#   of habitat.mx for each course/property/year.
# maxNBDets is a 3-D array giving the maximum number of traps for which an individual could be captured in. This is 
#   calculated based on locations of traps, potential locations of activity centers, and maximum possible trapping distance 
#   (844m). Dimensions refer to years (1st dimension), courses (2nd dimension), and properties (3rd dimension).
# detectorIndex is a 5-D array. For each 50x50m grid cells (1st dimension), gives the IDs of traps in which individuals
#   whose activity centers are within that grid cell could potentially be captured (2nd dimension). 3rd dimension
#   is course, 4th dimension is year, and 5th dimension is property.
# nDetectors is a 4-D array giving, for each 50x50m grid cell in n.cells (1st dimension) at each course (2nd dimension)
#   in each year (3rd dimension) on each property (4th dimension), the total number of traps in which individuals
#   whose activity centers are within that grid cell could potentially be captured. Corresponds to the number of 
#   non-zero columns in detectorIndex.
# habitatIDDet is a 5-D array giving the grid cell ID for each 50x50m grid cell in the gridded trapping 
#   area. Values of 0 denote a grid cell which cannot host activity centers (classified as wetland or too far
#   from traps to feasibly be captured). Dimensions correspond to maximum y-coordinates from y.max
#   (1st dimension) across courses/properties/years, maximum x-coordinates from x.max (2nd dimension)
#   across courses/properties/years, courses (3rd dimension), years (4th dimension), and properties
#   (5th dimension). Values of NA correspond to grids outside of the trapping area for the 
#   given course/year/property and don't affect the likelihood. This is the same as habitat.mx except
#   giving the grid cell ID rather than a binary variable.
# area is a 3-D array giving the area surveyed by the trapping grid (ha) at each course (2nd dimension),
#   year (1st dimension), and property (3rd dimension). This corresponds to the sum of habitat.mx per
#   course/year/property, multiplied by the cell size (0.25ha).
# s.I is a 5-D array giving values of activity center coordinates for initializing the MCMC sampler.
#   Dimensions correspond to the number of individuals in M (1st dimension), the x and y coordinates of
#   activity centers (1st and 2nd elements, respectively, of 2nd dimension), courses (3rd dimension),
#   years (4th dimension), and property (5th dimension). These values ensure that activity centers 
#   are within the trapping area and correspond to upland grid cells (habitat.mx=1). For captured birds,
#   ensures that activity centers are within 844m of all traps in which the birds were captured.
# z.I is a 4-D array giving the values of z for initializing the MCMC sampler. Dimensions correspond 
#   to the number of individuals in M (1st dimension), courses (2nd dimension), years (3rd dimension),
#   and property (4th dimension). Only need initial values for augmented individuals (set to NA for
#   individuals with z=1).
# constraint.lambda is a 3-D array of 1s for maximizing the constraint likelihood to ensure that expected
#   abundance (lambda) is not greater than M.  Dimensions correspond to years (1st dimension), courses
#   (2nd dimension), and properties (3rd dimension). Values of NA denote years/courses/properties without
#   data and do not affect the likelihood.


# Inputting nDetectors in the likelihood function 'dbinomLocal_normal' doesn't work with subset arrays.
#   Need to do for each course/year.


mod_data <- list(lowerCoords = data.trappingSCR$lowerCoords,
                 upperCoords = data.trappingSCR$upperCoords,
                 y = data.trappingSCR$y,
                 z = data.trappingSCR$z,
                 detector.xy = data.trappingSCR$detector.xy,
                 habitat.mx = data.trappingSCR$habitat.mx,
                 ones = data.trappingSCR$ones,
                 trials = data.trappingSCR$trials,
                 detectorIndex = data.trappingSCR$detectorIndex,
                 nDetectors.p1c1y1 = data.trappingSCR$nDetectors[,1,1,1],
                 nDetectors.p1c2y1 = data.trappingSCR$nDetectors[,2,1,1],
                 nDetectors.p1c3y1 = data.trappingSCR$nDetectors[,3,1,1],
                 nDetectors.p1c4y1 = data.trappingSCR$nDetectors[,4,1,1],
                 nDetectors.p1c1y2 = data.trappingSCR$nDetectors[,1,2,1],
                 nDetectors.p1c2y2 = data.trappingSCR$nDetectors[,2,2,1],
                 nDetectors.p1c3y2 = data.trappingSCR$nDetectors[,3,2,1],
                 nDetectors.p1c4y2 = data.trappingSCR$nDetectors[,4,2,1],
                 nDetectors.p2c1y1 = data.trappingSCR$nDetectors[,1,1,2],
                 nDetectors.p2c2y1 = data.trappingSCR$nDetectors[,2,1,2],
                 habitatIDDet = data.trappingSCR$habitatIDDet,
                 yDets = data.trappingSCR$yDets,
                 nbDetections = data.trappingSCR$nbDetections,
                 constraint.lambda = data.trappingSCR$constraint.lambda)

mod_constants <- list(Nproperty = data.trappingSCR$Nproperty,
                      Ncourse = data.trappingSCR$Ncourse,
                      Nyear = data.trappingSCR$Nyear,
                      M = data.trappingSCR$M,
                      sigindex = data.trappingSCR$sigindex,
                      n.detectors = data.trappingSCR$n.detectors,
                      y.max = data.trappingSCR$y.max,
                      x.max = data.trappingSCR$x.max,
                      n.cells = data.trappingSCR$n.cells,
                      maxNBDets = data.trappingSCR$maxNBDets,
                      nMaxDetectors = data.trappingSCR$nMaxDetectors,
                      area = data.trappingSCR$area)

mod_inits <- list(s = data.trappingSCR$s.I,
                  z = data.trappingSCR$z.I,
                  p0 = rep(0.1, times=2),
                  a0 = log(t(apply(apply(data.trappingSCR$z, c(2,3,4), sum, na.rm=T), c(1,3), max))+1),
                  a1 = 0,
                  sigma = rep(2, times=2))


Trapping_code <- nimbleCode({
  
  ## Priors for detection. Varying based on timing of supplemental feed restriction prior to trapping (2 levels)
  for(h in 1:2){
    sigma[h] ~ dunif(0, 50) # Distance-detection function for half-normal detection function
    p0[h] ~ dunif(0, 1) # probability of trapping bird at 0 distance, not assuming to be 1
  } # end h
  a1 ~ dnorm(0,0.001)
  
  ## Likelihood
  for(p in 1:Nproperty){
    for(c in 1:Ncourse[p]){
      
      a0[p,c] ~ dnorm(0,0.001)
      
      for(e in 1:Nyear[p]){
        
        lambda[e,c,p] <- exp(a0[p,c] + a1 * (e-1)) # Expected abundance
        constraint.lambda[e,c,p] ~ dconstraint(lambda[e,c,p] <= M[e,c,p]) # Abundance can't be higher than M
        psi[e,c,p] <- lambda[e,c,p] / M[e,c,p] # Inclusion probability for data augmentation
        
        # Looping over individuals
        for(d in 1:M[e,c,p]){
          
          # Activity center coordinates
          s[d,1,c,e,p] ~ dunif(0, x.max[e,c,p])
          s[d,2,c,e,p] ~ dunif(0, y.max[e,c,p])
          
          # Habitat constraint using nimbleSCR function. Based on the ones trick.
          ones[d,c,e,p] ~ dHabitatMask(s = s[d,1:2,c,e,p],
                                       xmin = lowerCoords[1],
                                       xmax = upperCoords[c,1,e,p],
                                       ymin = lowerCoords[2],
                                       ymax = upperCoords[c,2,e,p],
                                       habitat = habitat.mx[1:y.max[e,c,p], 1:x.max[e,c,p], c, e, p]) # Flipping x/y due to nature of the matrix construction
          
          # Latent inclusion alive state
          z[d,c,e,p] ~ dbern(psi[e,c,p])
        } #end d
        
        N[e,c,p] <- sum(z[1:M[e,c,p],c,e,p])    # Realized abundance
        EN[e,c,p] <- psi[e,c,p] * M[e,c,p]      # Expected abundance
        D[e,c,p] <- N[e,c,p] / area[e,c,p]      # Realized density
        ED[e,c,p] <- EN[e,c,p] / area[e,c,p]    # Expected density
        
      } #end e
    } # end c
  } # end p
  
  # Likelihood function not working with arrays, so need to do individually
  
  # Property 1, course 1, year 1
  for(a in 1:M[1,1,1]){
    y[a,1:nMaxDetectors[1,1,1],1,1,1] ~ dbinomLocal_normal(detNums = nbDetections[a,1,1,1],
                                                           detIndices = yDets[a,1:nMaxDetectors[1,1,1],1,1,1],
                                                           size = trials[1:n.detectors[1,1,1],1,1,1],
                                                           p0 = p0[sigindex[1,1,1]],
                                                           s = s[a,1:2,1,1,1],
                                                           sigma = sigma[sigindex[1,1,1]],
                                                           trapCoords = detector.xy[1:n.detectors[1,1,1],1:2,1,1,1],
                                                           localTrapsIndices = detectorIndex[1:n.cells[1,1,1],1:maxNBDets[1,1,1],1,1,1],
                                                           localTrapsNum = nDetectors.p1c1y1[1:n.cells[1,1,1]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[1,1,1],1:x.max[1,1,1],1,1,1],
                                                           indicator = z[a,1,1,1])
  } # end a
  
  # Property 1, course 2, year 1
  for(b in 1:M[1,2,1]){
    y[b,1:nMaxDetectors[1,2,1],2,1,1] ~ dbinomLocal_normal(detNums = nbDetections[b,2,1,1],
                                                           detIndices = yDets[b,1:nMaxDetectors[1,2,1],2,1,1],
                                                           size = trials[1:n.detectors[1,2,1],2,1,1],
                                                           p0 = p0[sigindex[1,2,1]],
                                                           s = s[b,1:2,2,1,1],
                                                           sigma = sigma[sigindex[1,2,1]],
                                                           trapCoords = detector.xy[1:n.detectors[1,2,1],1:2,2,1,1],
                                                           localTrapsIndices = detectorIndex[1:n.cells[1,2,1],1:maxNBDets[1,2,1],2,1,1],
                                                           localTrapsNum = nDetectors.p1c2y1[1:n.cells[1,2,1]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[1,2,1],1:x.max[1,2,1],2,1,1],
                                                           indicator = z[b,2,1,1])
  } # end b
  
  # Property 1, course 3, year 1
  for(f in 1:M[1,3,1]){
    y[f,1:nMaxDetectors[1,3,1],3,1,1] ~ dbinomLocal_normal(detNums = nbDetections[f,3,1,1],
                                                           detIndices = yDets[f,1:nMaxDetectors[1,3,1],3,1,1],
                                                           size = trials[1:n.detectors[1,3,1],3,1,1],
                                                           p0 = p0[sigindex[1,3,1]],
                                                           s = s[f,1:2,3,1,1],
                                                           sigma = sigma[sigindex[1,3,1]],
                                                           trapCoords = detector.xy[1:n.detectors[1,3,1],1:2,3,1,1],
                                                           localTrapsIndices = detectorIndex[1:n.cells[1,3,1],1:maxNBDets[1,3,1],3,1,1],
                                                           localTrapsNum = nDetectors.p1c3y1[1:n.cells[1,3,1]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[1,3,1],1:x.max[1,3,1],3,1,1],
                                                           indicator = z[f,3,1,1])
  } # end f
  
  # Property 1, course 4, year 1
  for(g in 1:M[1,4,1]){
    y[g,1:nMaxDetectors[1,4,1],4,1,1] ~ dbinomLocal_normal(detNums = nbDetections[g,4,1,1],
                                                           detIndices = yDets[g,1:nMaxDetectors[1,4,1],4,1,1],
                                                           size = trials[1:n.detectors[1,4,1],4,1,1],
                                                           p0 = p0[sigindex[1,4,1]],
                                                           s = s[g,1:2,4,1,1],
                                                           sigma = sigma[sigindex[1,4,1]],
                                                           trapCoords = detector.xy[1:n.detectors[1,4,1],1:2,4,1,1],
                                                           localTrapsIndices = detectorIndex[1:n.cells[1,4,1],1:maxNBDets[1,4,1],4,1,1],
                                                           localTrapsNum = nDetectors.p1c4y1[1:n.cells[1,4,1]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[1,4,1],1:x.max[1,4,1],4,1,1],
                                                           indicator = z[g,4,1,1])
  } # end g
  
  # Property 1, course 1, year 2
  for(i in 1:M[2,1,1]){
    y[i,1:nMaxDetectors[2,1,1],1,2,1] ~ dbinomLocal_normal(detNums = nbDetections[i,1,2,1],
                                                           detIndices = yDets[i,1:nMaxDetectors[2,1,1],1,2,1],
                                                           size = trials[1:n.detectors[2,1,1],1,2,1],
                                                           p0 = p0[sigindex[2,1,1]],
                                                           s = s[i,1:2,1,2,1],
                                                           sigma = sigma[sigindex[2,1,1]],
                                                           trapCoords = detector.xy[1:n.detectors[2,1,1],1:2,1,2,1],
                                                           localTrapsIndices = detectorIndex[1:n.cells[2,1,1],1:maxNBDets[2,1,1],1,2,1],
                                                           localTrapsNum = nDetectors.p1c1y2[1:n.cells[2,1,1]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[2,1,1],1:x.max[2,1,1],1,2,1],
                                                           indicator = z[i,1,2,1])
  } # end i
  
  # Property 1, course 2, year 2
  for(j in 1:M[2,2,1]){
    y[j,1:nMaxDetectors[2,2,1],2,2,1] ~ dbinomLocal_normal(detNums = nbDetections[j,2,2,1],
                                                           detIndices = yDets[j,1:nMaxDetectors[2,2,1],2,2,1],
                                                           size = trials[1:n.detectors[2,2,1],2,2,1],
                                                           p0 = p0[sigindex[2,2,1]],
                                                           s = s[j,1:2,2,2,1],
                                                           sigma = sigma[sigindex[2,2,1]],
                                                           trapCoords = detector.xy[1:n.detectors[2,2,1],1:2,2,2,1],
                                                           localTrapsIndices = detectorIndex[1:n.cells[2,2,1],1:maxNBDets[2,2,1],2,2,1],
                                                           localTrapsNum = nDetectors.p1c2y2[1:n.cells[2,2,1]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[2,2,1],1:x.max[2,2,1],2,2,1],
                                                           indicator = z[j,2,2,1])
  } # end j
  
  # Property 1, course 3, year 2
  for(k in 1:M[2,3,1]){
    y[k,1:nMaxDetectors[2,3,1],3,2,1] ~ dbinomLocal_normal(detNums = nbDetections[k,3,2,1],
                                                           detIndices = yDets[k,1:nMaxDetectors[2,3,1],3,2,1],
                                                           size = trials[1:n.detectors[2,3,1],3,2,1],
                                                           p0 = p0[sigindex[2,3,1]],
                                                           s = s[k,1:2,3,2,1],
                                                           sigma = sigma[sigindex[2,3,1]],
                                                           trapCoords = detector.xy[1:n.detectors[2,3,1],1:2,3,2,1],
                                                           localTrapsIndices = detectorIndex[1:n.cells[2,3,1],1:maxNBDets[2,3,1],3,2,1],
                                                           localTrapsNum = nDetectors.p1c3y2[1:n.cells[2,3,1]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[2,3,1],1:x.max[2,3,1],3,2,1],
                                                           indicator = z[k,3,2,1])
  } # end k
  
  # Property 1, course 4, year 2
  for(l in 1:M[2,4,1]){
    y[l,1:nMaxDetectors[2,4,1],4,2,1] ~ dbinomLocal_normal(detNums = nbDetections[l,4,2,1],
                                                           detIndices = yDets[l,1:nMaxDetectors[2,4,1],4,2,1],
                                                           size = trials[1:n.detectors[2,4,1],4,2,1],
                                                           p0 = p0[sigindex[2,4,1]],
                                                           s = s[l,1:2,4,2,1],
                                                           sigma = sigma[sigindex[2,4,1]],
                                                           trapCoords = detector.xy[1:n.detectors[2,4,1],1:2,4,2,1],
                                                           localTrapsIndices = detectorIndex[1:n.cells[2,4,1],1:maxNBDets[2,4,1],4,2,1],
                                                           localTrapsNum = nDetectors.p1c4y2[1:n.cells[2,4,1]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[2,4,1],1:x.max[2,4,1],4,2,1],
                                                           indicator = z[l,4,2,1])
  } # end l
  
  # Property 2, course 1, year 1
  for(n in 1:M[1,1,2]){
    y[n,1:nMaxDetectors[1,1,2],1,1,2] ~ dbinomLocal_normal(detNums = nbDetections[n,1,1,2],
                                                           detIndices = yDets[n,1:nMaxDetectors[1,1,2],1,1,2],
                                                           size = trials[1:n.detectors[1,1,2],1,1,2],
                                                           p0 = p0[sigindex[1,1,2]],
                                                           s = s[n,1:2,1,1,2],
                                                           sigma = sigma[sigindex[1,1,2]],
                                                           trapCoords = detector.xy[1:n.detectors[1,1,2],1:2,1,1,2],
                                                           localTrapsIndices = detectorIndex[1:n.cells[1,1,2],1:maxNBDets[1,1,2],1,1,2],
                                                           localTrapsNum = nDetectors.p2c1y1[1:n.cells[1,1,2]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[1,1,2],1:x.max[1,1,2],1,1,2],
                                                           indicator = z[n,1,1,2])
  } # end n
  
  # Property 2, course 2, year 1
  for(q in 1:M[1,2,2]){
    y[q,1:nMaxDetectors[1,2,2],2,1,2] ~ dbinomLocal_normal(detNums = nbDetections[q,2,1,2],
                                                           detIndices = yDets[q,1:nMaxDetectors[1,2,2],2,1,2],
                                                           size = trials[1:n.detectors[1,2,2],2,1,2],
                                                           p0 = p0[sigindex[1,2,2]],
                                                           s = s[q,1:2,2,1,2],
                                                           sigma = sigma[sigindex[1,2,2]],
                                                           trapCoords = detector.xy[1:n.detectors[1,2,2],1:2,2,1,2],
                                                           localTrapsIndices = detectorIndex[1:n.cells[1,2,2],1:maxNBDets[1,2,2],2,1,2],
                                                           localTrapsNum = nDetectors.p2c2y1[1:n.cells[1,2,2]],
                                                           resizeFactor = 1,
                                                           habitatGrid = habitatIDDet[1:y.max[1,2,2],1:x.max[1,2,2],2,1,2],
                                                           indicator = z[q,2,1,2])
  } # end q
})


monitors <- c("psi", "N", "EN", "D", "ED", "p0", "sigma", "a0", "a1")


Rmodel <- nimbleModel(code=Trapping_code, constants=mod_constants, data=mod_data, inits=mod_inits)

# Customize the mcmc sampler for s. By default, applies a Univaiate sampler to each dimension of
#   s (x separately from y), assigning a joint sampler
conf <- configureMCMC(Rmodel, monitors = monitors, print = FALSE)
conf$removeSamplers(c("s"))
ACnodes <- c(paste0("s[",1:mod_constants$M[1,1,1],",1:2,1,1,1]"),
             paste0("s[",1:mod_constants$M[1,2,1],",1:2,2,1,1]"),
             paste0("s[",1:mod_constants$M[1,3,1],",1:2,3,1,1]"),
             paste0("s[",1:mod_constants$M[1,4,1],",1:2,4,1,1]"),
             paste0("s[",1:mod_constants$M[2,1,1],",1:2,1,2,1]"),
             paste0("s[",1:mod_constants$M[2,2,1],",1:2,2,2,1]"),
             paste0("s[",1:mod_constants$M[2,3,1],",1:2,3,2,1]"),
             paste0("s[",1:mod_constants$M[2,4,1],",1:2,4,2,1]"),
             paste0("s[",1:mod_constants$M[1,1,2],",1:2,1,1,2]"),
             paste0("s[",1:mod_constants$M[1,2,2],",1:2,2,1,2]"))
for(node in ACnodes) {
  conf$addSampler(target = node,
                  type = "RW_block",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

trappingout <- runMCMC(Cmcmc, niter = 30000, nchain=3, nburnin=3000, thin=5, samplesAsCodaMCMC=TRUE)
