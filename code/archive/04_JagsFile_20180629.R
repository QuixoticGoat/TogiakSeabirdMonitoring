
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
    sigma.loglambda ~ dunif(0,1)               # Variance on log lambda for pop process model
    sigma2.loglambda <- pow(sigma.loglambda, 2)
    tau.loglambda <- pow(sigma.loglambda, -2)

    for(i in 1:P){                              # For each plot (P),
        mu.loglambda[i] ~ dunif(-100,100)         # Avg. log lambda is normally dist. around zero
        N[i,1] ~ dunif(0,1)                  # Abundance at time 1 is a categorical dist, drawn from probs in psi
    } #i

    for(t in 1:T){                               # For each year (T),
        p[t] ~ dunif(0,1)                       # Observation prob. is uniformly distributed between 0 and 1 
    } #t

    ## Population process model
    for(i in 1:P){                               # For each plot,
        for(t in 2:T){                             # For years 2 to 28,
            loglambda[i,t-1] <- mu.loglambda[i]
            lambda[i,t-1] <- exp(loglambda[i,t-1])
            N[i,t] ~ dpois(N[i,t-1])
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
        mu.lambda[i] <- exp(mu.loglambda[i])         # Derive mean lambda
    } #i
    }
    ")

#-------------------------------------------------------------------------------

# Load packages:
library(tidyverse)
library(jagsUI)

df.jags <- df.list$df.jags

## Bundle data
# Array of counts:
C <- split(df.jags[, 4:15], as.factor(df.jags$Year))
C <- array(as.numeric(unlist(C)), dim = c(nrow(C[[1]]), ncol(C[[1]]), length(C)))

nplots <- length(unique(df.jags$PlotID))
nreps <- ncol(df.jags[4:15])
nyears <- length(unique(df.jags$Year))

# psi (for dcat distribution) based on the max value (count) in each row (plot) for each dimension (year)
maxCyear1 <- apply(C, c(1,3), max, na.rm = T)[,1]  # Max count in year 1
pi <- matrix(0, 400, nplots)                     # Empty 600x600 matrix
for(i in 1:nplots){                               # Fill matrix with 
  pi[1:(maxCyear1[i]-1), i] <- 0
  pi[maxCyear1[i]:400, i] <- 1/(401-maxCyear1[i])
}

dat <- list(C = C, T = nyears, J = nreps, P = nplots, pi = pi)                                

## Initial values
inits <- function() list(mu.loglambda = rep(0, 30),
                         sigma.loglambda = runif(1, 0, 10),
                         N = matrix(c(maxCyear1, 
                                      rep(NA, 30*27)), 30, 28),
                         p = rep(0.5, 28)
                         )

## Parameters monitored
pars <- c("N", "lambda", "p")

## MCMC settings
ni <- 15000; nt <- 20; nb <- 5000; nc <- 1

## Call JAGS from R
mod <- jags(dat, inits = inits, pars, "./code/jags/mod1.jags", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


## Summarize posteriors
print(mod, digits = 3)
plot(mod)

