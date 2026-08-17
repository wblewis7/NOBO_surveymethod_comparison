# R and NIMBLE code for estimating fall covey abundance of northern bobwhite (Colinus
#   virginianus) via spatial capture-recapture quadrat surveys in Central Florida, USA 
#   2021-2022.
# From Schmitt, R. L., A. S. Cramer, G. Beane, W. B. Lewis, and J. A. Martin. A 
#   comparison of northern bobwhite density estimators.
# Contact: Will Lewis, University of Georgia, william.lewis27@uga.edu




# Fall bobwhite coveys were surveyed in the falls of 2021/2022 in Central Florida, USA.
#   Surveys occurred at 2 properties, each of which contained 2 - 4 survey courses.
#   Each course has 2 - 4 survey locations, and each was surveyed twice in a given fall.
# Quadrat surveys consisted of 4 trained observers simultaneously counting bobwhite coveys
#   within a 500x500 m square survey area. Each observer stood at the midpoint on one of the
#  four sides and recorded the location of all detected coveys within the 25ha area. After
#  surveying, the four observers reconciled their observations into a single record of how
#  many coveys were detected and which observers detected each covey.
# We estimate density from reconciled surveys using a hierarchical spatially-explicit capture-recapture
#   model (Royle & Young. 2008. Ecology) with modifications to account for multiple observers.
#   We take a data augmentation approach, modeling whether or not each augmented covey at each sampling unit in 
#   each year is a true covey (z) based on the inclusion parameter (psi). We model expected abundance at the 
#   sampling unit (lambda) based on survey course and year, then convert to psi by dividing lambda by the total number
#   of possible individuals (M). We also model the spatial locations of each covey bounded within the limits of the
#   25ha survey area. We decompose the detection process into availability (probability of a covey calling during
#   a survey) and conditional detection probability (probability of observers detecting a calling covey). For each survey visit, 
#   we model availability of coveys based on z and property-specific probability of availability. This formulation fixes
#   abundance across repeat survey visits at the same sampling unit in each year, but allows for variable availability and spatial locations across visits.
#   For each of the 4 observers, we then calculate detection probability for each covey based on a half-normal 
#   detection function and the distance between the covey and the observer. We then model detections for each
#   potential covey by each observer, for each visit at each sampling unit in each year, based on a Bernoulli distribution.



library(nimble)


# Loading in data
load("NOBO_CFL_modeldata_quadratSCR.gzip")
str(data.quadratSCR)

# Nproperty is the total number of study properties (2).
# Ncourse is the total number of survey courses across the 2 properties (6).
# Npoints is the total number of sampling units (19).
# Nyear is the total number of survey years (2).
# course gives the course number for each sampling unit in Npoints.
# Nvisits gives the number of repeat surveys performed at each study unit in each year (2).
# property is a vector giving the property for each study course in Npoints.
# M is the maximum number of potential coveys at each sampling unit in each year for data augmentation (100).
# Xlim is a matrix giving the lower/upper (columns) x-coordinates of the 25ha survey quadrats
#   at each sampling unit in Npoints (rows).
# Ylim is a matrix giving the lower/upper (columns) y-coordinates of the 25ha survey quadrats
#   at each sampling unit in Npoints (rows).
# J is the number of observers for each quadrat survey (4).
# X is a 3-D array giving the x (1st element of 3rd dimension) and y (2nd element of 3rd dimension)
#   of each observer in J (1st dimension) during quadrat surveys at each sampling unit in
#   Npoints (2nd dimension).
# z is an M x Npoints x Nyear array fixing the z-state to 1 for detected coveys during data augmentation,
#   with NA for z-states to be estimated (augmented indidviduals).
# y.Q is a 5-D array of observed covey detections denoting whether or not each covey in M (1st dimension)
#   was detected at each sampling unit in Npoints (2nd detection), by each observer in J (3rd dimension),
#   for each year in Nyear (4th dimension) and each repeat survey in Nvisits (5th dimension).
# constraint.lambda is a Npoints x Nyear matrix of 1s for maximizing the contraint likelihood to ensure
#   that expected abundance (lambda) is not greater than M.


mod_data <- list(y.Q = data.quadratSCR$y.Q,
                 z = data.quadratSCR$z,
                 constraint.lambda = data.quadratSCR$constraint.lambda)

mod_constants <- data.quadratSCR[!(names(data.quadratSCR) %in% c("y.Q","z","constraint.lambda"))]

# Setting up initial values for z-state, availability, and s
z.init <- array(rbinom(prod(dim(data.quadratSCR$z)),1,0.1), dim=dim(data.quadratSCR$z))
for(h in 1:ncol(z.init)){
  for(s in 1:dim(z.init)[3]){
    z.init[1:sum(!is.na(data.quadratSCR$z[,h,s])),h,s] <- NA
  }
}
z.A.init <- array(rbinom(prod(dim(data.quadratSCR$z))*data.quadratSCR$Nvisits,1,0.9), dim=c(dim(data.quadratSCR$z),data.quadratSCR$Nvisits))
for(d in 1:dim(z.A.init)[4]){
  z.A.init[,,,d] <- z.A.init[,,,d] * z.init
}
z.A.init[is.na(z.A.init)] <- 1
s.init <- array(NA, dim=c(data.quadratSCR$M, data.quadratSCR$Npoints, data.quadratSCR$Nyear, data.quadratSCR$Nvisits, 2))
for(g in 1:data.quadratSCR$Npoints){
  s.init[,g,,,1] <- runif(data.quadratSCR$M * data.quadratSCR$Nyear * data.quadratSCR$Nvisits, data.quadratSCR$Xlim[g,1], data.quadratSCR$Xlim[g,2])
  s.init[,g,,,2] <- runif(data.quadratSCR$M * data.quadratSCR$Nyear * data.quadratSCR$Nvisits, data.quadratSCR$Ylim[g,1], data.quadratSCR$Ylim[g,2])
}
mod_inits <- list(a0 = runif(data.quadratSCR$Ncourse, qlogis(0.05), qlogis(0.15)),
                  a1 = runif(1, -0.5, 0.5),
                  A = runif(data.quadratSCR$Nproperty, 0.1, 0.9),
                  z = z.init,
                  z.A = z.A.init,
                  sigma = runif(1, 100, 400),
                  s = s.init)



Quadrat_code <- nimbleCode({
  
  # Priors
  for(a in 1:Ncourse){
    a0[a] ~ dnorm(0, 0.001)
  } # end a
  a1 ~ dnorm(0, 0.001)
  for(o in 1:Nproperty){
    A[o] ~ dunif(0, 1)
  } # end o
  sigma ~ dunif(0,500) # Effect of distance (m) on detection prob (p)
  
  
  for(g in 1:Npoints){
    for(r in 1:Nyear){
      
      # Inclusion parameter for data augmentation modeled with study course and year effects
      lambda[g,r] <- exp(a0[course[g]] + a1 * (r-1)) # Expected abundance
      constraint.lambda[g,r] ~ dconstraint(lambda[g,r] <= M) # Abundance can't be higher than M
      psi[g,r] <- lambda[g,r] / M # Inclusion parameter for data augmentation
      
      # Derived abundance
      N[g,r] <- sum(z[,g,r])
      
      # Looping through potential coveys 1:M
      for(l in 1:M){
        
        z[l,g,r] ~ dbern(psi[g,r])
        
        # Looping through survey visits
        for(q in 1:Nvisits){
          
          s[l,g,r,q,1] ~ dunif(Xlim[g,1], Xlim[g,2]) # x-coordinate of true covey location
          s[l,g,r,q,2] ~ dunif(Ylim[g,1], Ylim[g,2]) # y-coordinate of true covey location
          
          z.A[l,g,r,q] ~ dbern(A[property[g]] * z[l,g,r])
          
          # Looping through each observer
          for(j in 1:J){
            
            d.1[l,g,j,r,q] <- sqrt((s[l,g,r,q,1] - X[j,g,1])^2 +
                                     (s[l,g,r,q,2] - X[j,g,2])^2) # distance between covey location and observer
            p[l,g,j,r,q] <- exp(-d.1[l,g,j,r,q]^2/(2*sigma^2)) # detection probability of covey l by observer j
            
            # Data Likelihood
            y.Q[l,g,j,r,q] ~ dbern(p[l,g,j,r,q] * z.A[l,g,r,q])
            
          } # end j
        } # end q
      } # end l
    } # end r
  } # end g
})


monitors <- c("sigma", "a0", "a1", "N", "A")


quadmodel <- nimbleModel(code = Quadrat_code,
                     data = mod_data,
                     constants = mod_constants,
                     inits = mod_inits)
quadout <- nimbleMCMC(model = quadmodel,
                       niter = 30000, nchains = 3, nburnin = 3000, thin=5,
                       samplesAsCodaMCMC=TRUE, monitor=monitors)
