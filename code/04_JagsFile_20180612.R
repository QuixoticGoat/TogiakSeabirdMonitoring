
################################################################################
# Creates a .jag file that defined the JAGS model for environmental effects on #
# seabird and nest abundance                                                   #
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov>                                    #  
# Date modified: 6/27/2018                                                     #
################################################################################


# Jags model file
cat(file = "./code/jags/mod1.jags", "
    model {
    
    ## Priors
    sd.loglambda ~ dnorm(0, 0.5)                 # Variance on lambda for pop process model
    tau.loglambda <- pow(sd.loglambda, -2)
    
    for(i in 1:P){                               # For each plot (P),
        mu.loglambda[i] ~ dnorm(0, 0.01)         # Avg. log lambda is normally dist. around zero
        N[i,1] ~ dcat(psi[ ,i])                  # Abundance at time 1 is a categorical dist, drawn from probs in psi
    } #i

    for(t in 1:T){                               # For each year (T),
        p[t] ~ dunif(0, 1)                       # Observation prob. is uniformly distributed between 0 and 1 
    } #t

    ## Population process model
    for(i in 1:P){                               # For each plot,
        for(t in 2:T){                             # For years 2 to 28,
            loglambda[i,t-1] <- mu.loglambda[i] + eps.loglambda[i,t-1]
            eps.loglambda[i, t-1] ~ dnorm(0, tau.loglambda)
            lambda[i,t-1] <- exp(loglambda[i,t-1])
            N[i,t] ~ dpois(N[i,t-1] * lambda[i,t-1])
        } #t
    } #i

    ## Observation process model
    for(i in 1:P){                                   # For each plot,
        for(t in 1:T){                                # For each year,
            for(j in 1:J){                             # For each rep,
                C[i,j,t] ~ dbin(p[t],N[i,t])
            } #j
        } #t
    } #i
    
    ## Derived parameters
    for(i in 1:P){                                   # For each plot,
        mu.lambda[i] <- exp(mu.loglambda[i])            # lambda is the exp. of mu.loglambda
    } #i
    }
    ")

#-------------------------------------------------------------------------------

# Load packages:
library(tidyverse)
library(jagsUI)

## Bundle data
# Array of counts:
C <- split(df.jags[, 4:15], as.factor(df.jags$Year))
C <- array(as.numeric(unlist(C)), dim = c(nrow(C[[1]]), ncol(C[[1]]), length(C)))

nplots <- length(unique(df.jags$PlotID))
nreps <- ncol(df.jags[4:15])
nyears <- length(unique(df.jags$Year))

# psi (for dcat distribution) based on the max value (count) in each row (plot) for each dimension (year)
maxCyear1 <- apply(C, c(1,3), max, na.rm = T)[,1]  # Max count in year 1 -- how to deal with NAs (years when plots were not surveyed?)
psi <- matrix(0, 400, nplots)                     # Empty 600x600 matrix (why 600?)
for(i in 1:nplots){                               # Fill matrix with 
  psi[1:(maxCyear1[i]-1), i] <- 0
  psi[maxCyear1[i]:400, i] <- 1/(401-maxCyear1[i])
}

dat <- list(C = C, T = nyears, J = nreps, P = nplots, psi = psi)                                

## Initial values
inits <- function() list(mu.loglambda = runif(30, -0.1, 0.1),
                         sd.loglambda = runif(1, 0, 1),
                         N = matrix(c(round(runif(30, maxCyear1, 400)), 
                                      rep(NA, 30*27)), 30, 28),
                         p = runif(28, 0.2, 1)
                         )

## Parameters monitored
pars <- c("N", "mu.loglambda", "sd.loglambda", "mu.lambda", "lambda", "p")

## MCMC settings
ni <- 15000; nt <- 20; nb <- 5000; nc <- 1

## Call JAGS from R
mod <- jags(dat, inits = NULL, pars, "./code/jags/mod1.jags", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


## Summarize posteriors
print(mod, digits = 3)
plot(mod)

