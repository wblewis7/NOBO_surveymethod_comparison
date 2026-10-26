# R and NIMBLE code for estimating fall hunting encounter rates fornorthern bobwhite
#   (Colinus virginianus) coveys in Central Florida, USA 2021-2022.
# From Schmitt, R. L., A. S. Cramer, G. Beane, W. B. Lewis, and J. A. Martin. A 
#   comparison of northern bobwhite density estimators.
# Contact: Will Lewis, University of Georgia, william.lewis27@uga.edu




# Detailed hunting records were recorded on a private property in Central Florida, USA, during
#   the overwinter hunting seasons starting in late 2021 and 2022. The property was divided
#   into 3 survey courses. The hunting party consisted of two shooters and two dog handlers,
#   plus an independent observer. The observer recorded the number of coveys encountered, the
#   hunt duration (not including time spent changing dogs), and the average temperature (C)
#   during the hunt.
# We model the number of coveys encountered per hunt via a Poisson generalized linear model.
#   We include effects of hunting season (2021/2022 vs. 2022/2023) and temperature (standardized)
#   on the expected number of coveys encountered, plus an offset of log-transformed hunt duration
#   to account for differing survey duration. We do not model course effects on covey encounters
#   due to insufficient sample sizes (i.e., <10 hunts on some courses). We model effects of
#   temperature because previous research (Sisson & Stribling. 2009. Tall Timbers Research 
#   Station and Land Conservancy, Wellendorf et al. 2012. Proceedings of the National Quail
#   Symposium) has found covey activity and encounter rates to be lower in cooler temperatures.
# Hunting methods were largely consistent across years, including using the same independent
#   observer. As such, we interpret year effects on covey encounter rates as reflective of
#   changes in abundance between years.



library(nimble)


# Loading in data
load("NOBO_CFL_modeldata_huntingencounters.gzip")
str(data.huntingencounters)

# Encounters gives the number of coveys encountered on each hunt.
# Nhunt gives the total number of hunts (48).
# year gives the year of the hunting season associated with each hunt in Nhunt. Values
#  of 0 correspond to the 2021/2022 hunting season while values of 1 correspond to the
#   2022/2023 season.
# temp gives the average temperature (standardized) during each hunt in Nhunt.
# loghours gives the hunt duration in log-transformed hours of each hunt in Nhunt.



mod_data <- list(H = data.huntingencounters$Encounters)

mod_constants <- data.huntingencounters[!(names(data.huntingencounters) %in% c("Encounters"))]

mod_inits <- list(a0 = runif(1, -1, 1),
                  a1 = runif(1, -0.5, 0.5),
                  a2 = runif(1, -0.5, 0.5))



huntingencounters_code <- nimbleCode({
  
  # Priors
  a0 ~ dnorm(0, 0.001)
  a1 ~ dnorm(0, 0.001)
  a2 ~ dnorm(0, 0.001)
  
  
  # Covey encounter model
  for(w in 1:Nhunt){
    
    lambda[w] <- exp(a0 + a1 * year[w] + a2 * temp[w] + loghours[w])
    H[w] ~ dpois(lambda[w])
    
  } # end w
})


monitors <- c("a0", "a1", "a2")


huntingmodel <- nimbleModel(code = huntingencounters_code,
                            data = mod_data,
                            constants = mod_constants,
                            inits = mod_inits)
huntingout <- nimbleMCMC(model = huntingmodel,
                         niter = 30000, nchains = 3, nburnin = 3000, thin=5,
                         samplesAsCodaMCMC=TRUE, monitor=monitors)
