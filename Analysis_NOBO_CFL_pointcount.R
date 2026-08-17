# R and NIMBLE code for estimating fall covey abundance of northern bobwhite (Colinus
#   virginianus) via distance-sampling point counts in Central Florida, USA 
#   2021-2022.
# From Schmitt, R. L., A. S. Cramer, G. Beane, W. B. Lewis, and J. A. Martin. A 
#   comparison of northern bobwhite density estimators.
# Contact: Will Lewis, University of Georgia, william.lewis27@uga.edu




# Fall bobwhite coveys were surveyed in the falls of 2021/2022 in Central Florida, USA.
#   Surveys occurred at 2 properties, each of which contained 2 - 4 survey courses.
#   Each course has 2 - 4 sampling units, and each was surveyed twice in a given fall.
# On point counts, trained researchers recorded every calling bobwhite within 500m of the
#   sampling unit center and the distance to the covey. Surveys were conducted from 45-minutes
#   before sunrise to sunrise.
# We estimate density from surveys using the distance sampling models of Royle et al. 2004. 
#   Ecology, with modifications to account for repeat surveys described in Lewis et al. 2025.
#   Ecological Solutions and Evidence.
#  We model year effects on abundance and course-level intercepts.



library(nimble)


# Loading in data
load("NOBO_CFL_modeldata_pointcount.gzip")
str(data.pointcount)

# Npoints is the number of sampling units across properties (19).
# Nyrs is the total number of survey years (2).
# Nvisits gives the number of repeat surveys performed at each sampling unit in each year (2).
# y is a 3-D array giving the total number of coveys detected at each repeat survey (2nd dimension)
#   for each sampling unit (1st dimension) in each year (3rd dimension)
# ydb is a 4-D array giving the number of covey detections in each distance bin (2nd dimension)
#   for each repeat survey (3rd dimension), at each sampling unit (1st dimension), in each year
#   (4th dimension).
# nB gives the total number of distance bins for distance sampling (5).
# db gives the boundaries of the distance bins (m).
# pix gives the proportional area of each distance bin.
# point.area gives the area (square meters) surveyed by each point count.
# ncourse gives the total number of survey courses (6).
# course gives the course number for each sampling unit in Npoints.



mod_data <- list(y = data.pointcount$y,
                 ydb = data.pointcount$ydb)

mod_constants <- data.pointcount[!(names(data.pointcount) %in% c("y","ydb"))]

mod_inits <- list(a0 = runif(data.pointcount$Ncourse, -1, 1),
                  a1 = runif(1, -0.5, 0.5),
                  beta = runif(1, log(150), log(250)),
                  N = apply(data.pointcount$y, c(1,3), max)+1)



PC_code <- nimbleCode({
  
  # Priors
  for(c in 1:Ncourse){
    a0[c] ~ dnorm(0, 0.001)
  } # end c
  a1 ~ dnorm(0, 0.001)
  beta ~ dunif(0, 10)
  log(sigma) <- beta
  
  
  # Conditional detection probability
  for(b in 1:nB){
    p.bin.1[b] <- (sigma^2*(1-exp(-db[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-db[b]^2/(2*sigma^2))))*2*3.1416/(point.area*pix[b]) # Half-normal detection function
    p.bin[b] <- p.bin.1[b]*pix[b]
  } # end b
  pd <- sum(p.bin[1:nB])
  
  # Point count detection model
  for(i in 1:Npoints){
    for(t in 1:Nyear){
      
      # Abundance state model
      lambda[i,t] <- exp(a0[course[i]] + a1 * (t-1))
      N[i,t] ~ dpois(lambda[i,t])
      
      # Observation process
      for(v in 1:Nvisits){
        
        # Likelihood for total number of detections
        y[i,v,t] ~ dbin(pd, N[i,t])
        
        # Likelihood for distance bin of detections
        ydb[i,1:nB,v,t] ~ dmulti(p.bin[1:nB], y[i,v,t])
      } # end v
    } # end t
  } # end i
})


monitors <- c("sigma","a0","a1", "N", "pd")


pcmodel <- nimbleModel(code = PC_code,
                     data = mod_data,
                     constants = mod_constants,
                     inits = mod_inits)
pcout <- nimbleMCMC(model = pcmodel,
                       niter = 30000, nchains = 3, nburnin = 3000, thin=5,
                       samplesAsCodaMCMC=TRUE, monitor=monitors)
