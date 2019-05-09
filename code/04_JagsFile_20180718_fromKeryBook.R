


#-------------------------------------------------------------------------------
# Open-population binomial-mixture model



# Get the data and put them into 3D array
# bdat <- read.table("./InputData/fritillary.txt", header = TRUE)
# y <- array(NA, dim = c(95, 2, 7))	# 95 sites, 2 reps, 7 days
# 
# for(k in 1:7){
#   sel.rows <- bdat$day == k
#   y[,,k] <- as.matrix(bdat)[sel.rows, 3:4]
# }
# View(y)                          # Look at data set in 3D layout
# str(y)
# 
# # Have a look at raw data
# day.max <- apply(y, c(1, 3), max, na.rm = TRUE)  # Max count each site and day
# day.max
# site.max <- apply(day.max, 1, max, na.rm = TRUE) # Max count each site
# site.max
# table(site.max)         # Frequency distribution of max counts
# plot(table(site.max))
# table(site.max>0)       # Observed occupancy is only 56%
# 
# # Sum of observed max as conventional estimator of total abundance
# max1 <- apply(y, c(1, 3), max)
# obs.max.sum <- apply(max1, 2, sum, na.rm = TRUE)
# obs.max.sum


#-------------------------------------------------------------------------------
# Simple N-mixture model

# Specify model in BUGS language
sink("Nmix0.jags")
cat("
    model {
    
    # Priors
    for (k in 1:28){
    alpha.lam[k] ~ dnorm(0, 0.01)
    p[k] ~ dunif(0, 1)
    }
    
    # Likelihood
    # Ecological model for true abundance
    for (k in 1:28){                          # Loop over 28 years
    lambda[k] <- exp(alpha.lam[k])
    for (i in 1:R){                       # Loop over R sites (30)
    N[i,k] ~ dpois(lambda[k])          # Abundance
    
    # Observation model for replicated counts
    for (j in 1:T){                    # Loop over temporal reps (12)
    y[i,j,k] ~ dbin(p[k],N[i,k])   # Detection
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic E for observed data
    eval[i,j,k] <- p[k] * N[i,k]   	# Expected values
    E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j,k] ~ dbin(p[k], N[i,k])
    E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
    
    } #j 
    } #i 
    } #k 
    
    # Derived and other quantities
    for (k in 1:28){
    totalN[k] <- sum(N[,k])	# Total pop. size across all sites
    mean.abundance[k] <- exp(alpha.lam[k])
    }
    fit <- sum(E[,,])
    fit.new <- sum(E.new[,,])
    }
    ",fill = TRUE)
sink()

# Bundle data
R = nrow(C)
T = ncol(C)
win.data <- list(y = C, R = R, T = T)

# Initial values
Nst <- apply(C, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(28, -1, 1))}

# Parameters monitored
params <- c("totalN", "mean.abundance", "alpha.lam", "p", "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 8
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 1 min)
out0 <- jags(win.data, inits, params, "Nmix0.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(out0, dig = 3)

# Evaluation of fit
plot(out0$sims.list$fit, out0$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(out0$sims.list$fit.new > out0$sims.list$fit)
mean(out0$mean$fit) / mean(out0$mean$fit.new)



#-------------------------------------------------------------------------------
# Zero-inflated N-mixture model (ZIP binomial-mixture model)

# Specify model in BUGS language
sink("Nmix1.jags")
cat("
    model {
    
    # Priors
    omega ~ dunif(0, 1)
    for (k in 1:28){
    alpha.lam[k] ~ dnorm(0, 0.01)
    p[k] ~ dunif(0, 1)
    }
    
    # Likelihood
    # Ecological model for true abundance
    for (i in 1:R){                          # Loop over R sites (30)
    z[i] ~ dbern(omega)                      # Latent suitability state
    for (k in 1:28){                         # Loop over survey periods (years)
    N[i,k] ~ dpois(lam.eff[i,k])             # Latent abundance state
    lam.eff[i,k] <- z[i] * lambda[i,k]
    log(lambda[i,k]) <- alpha.lam[k]
    
    # Observation model for replicated counts
    for (j in 1:T){                    # Loop over temporal reps (12)
    y[i,j,k] ~ dbin(p[k], N[i,k])   # Detection
    # Assess model fit using Chi-squared discrepancy
    
    # Compute fit statistic for observed data
    eval[i,j,k] <- p[k] * N[i,k]
    E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j,k] ~ dbin(p[k], N[i,k])
    E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
    } #j
    } #k
    } #i
    
    # Derived and other quantities
    for (k in 1:28){
    totalN[k] <- sum(N[,k])	# Estimate total pop. size across all sites
    mean.abundance[k] <- exp(alpha.lam[k])
    }
    fit <- sum(E[,,])
    fit.new <- sum(E.new[,,])
    }
    ",fill = TRUE)
sink()

# Bundle data
R = nrow(C)
T = ncol(C)
win.data <- list(y = C, R = R, T = T)

# Initial values
# It is advisable to give initial values for the latent state z, the best option is to provide a vector of 1
Nst <- apply(C, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(28, -1, 1), z = rep(1, 30))}

# Parameters monitored
params <- c("omega", "totalN", "alpha.lam", "p", "mean.abundance", "fit", "fit.new")

# MCMC settings
ni <- 30000
nt <- 15
nb <- 15000
nc <- 3

# Call JAGS from R (BRT 3 min)
out1 <- jags(win.data, inits, params, "Nmix1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(out1, dig = 3)

# Evaluation of fit
plot(out1$sims.list$fit, out1$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(out1$sims.list$fit.new > out1$sims.list$fit)
mean(out1$mean$fit) / mean(out1$mean$fit.new)


#-------------------------------------------------------------------------------
# N-mixture model with overdispersion in both abundance and detection

# Specify model in BUGS language
sink("Nmix2.jags")
cat("
    model{
    # Priors
    for (k in 1:28){
    alpha.lam[k] ~ dnorm(0, 0.1)
    beta[k] ~ dnorm(0, 0.1)
    }
    
    # Abundance site and detection site-by-day random effects
    for (i in 1:R){
    eps[i] ~ dnorm(0, tau.lam)                    # Abundance noise
    }
    tau.lam <- 1 / (sd.lam * sd.lam)
    sd.lam ~ dunif(0, 3)
    tau.p <- 1 / (sd.p * sd.p)
    sd.p ~ dunif(0, 3)
    
    # Likelihood
    # Ecological model for true abundance
    for (i in 1:R){                                 # Loop over R sites (30)
    for (k in 1:28){                              # Loop over year (28)
    N[i,k] ~ dpois(lambda[i,k])               # Abundance
    log(lambda[i,k]) <- alpha.lam[k] + eps[i]
    
    # Observation model for replicated counts
    for (j in 1:T){                           # Loop over temporal reps (12)
    y[i,j,k] ~ dbin(p[i,j,k], N[i,k])      # Detection
    p[i,j,k] <- 1 / (1 + exp(-lp[i,j,k])) 
    lp[i,j,k] ~ dnorm(beta[k], tau.p) # random delta defined implicitly
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j,k] <- p[i,j,k] * N[i,k]
    E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j,k] ~ dbin(p[i,j,k], N[i,k])
    E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
    } #j
    ik.p[i,k] <- mean(p[i,,k])
    } #k
    } #i
    
    # Derived and other quantities
    for (k in 1:28){
    totalN[k] <- sum(N[,k])   # Estimate total pop. size across all sites
    mean.abundance[k] <- mean(lambda[,k])
    mean.N[k] <- mean(N[,k])
    mean.detection[k] <- mean(ik.p[,k])
    }
    fit <- sum(E[,,])
    fit.new <- sum(E.new[,,])
    }
    ",fill = TRUE)
sink()

# Bundle data
R = nrow(C)
T = ncol(C)
win.data <- list(y = C, R = R, T = T)

# Initial values
Nst <- apply(C, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(28, -3, 3), 
                         beta = runif(28, -3, 3), sd.lam = runif(1, 0, 1), 
                         sd.p = runif(1, 0, 1))}

# Parameters monitored
params <- c("totalN", "alpha.lam", "beta", "sd.lam", "sd.p", "mean.abundance", 
            "mean.N", "mean.detection", "fit", "fit.new")

# MCMC settings
ni <- 350000
nt <- 300
nb <- 50000
nc <- 3

# Call JAGS from R (BRT 215 min)
out2 <- jags(win.data, inits, params, "Nmix2.jags", n.chains = nc, n.thin = nt, 
             n.iter = ni, n.burnin = nb)

# Evaluation of fit
plot(out2$sims.list$fit, out2$sims.list$fit.new, main = "", 
     xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", 
     frame.plot = FALSE, xlim = c(50, 200), ylim = c(50, 200))
abline(0, 1, lwd = 2, col = "black")
mean(out2$fit.new > out2$sims.list$fit)
mean(out2$mean$fit) / mean(out2$mean$fit.new)

# Summarize posteriors
print(out2, dig = 2)

max.day.count <- apply(C, c(1, 3), max, na.rm = TRUE)
max.day.count[max.day.count == "-Inf"] <- NA
mean.max.count <- apply(max.day.count, 2, mean, na.rm = TRUE)
mean.max.count

par(mfrow = c(2, 1))
plot(1:28, mean.max.count, xlab = "Day", ylab = "Mean daily abundance", las = 1, ylim = c(0, 16), type = "b", main = "", frame.plot = FALSE, pch = 16, lwd = 2)
lines(1:28, out2$summary[25:31,5], type = "b", pch = 16, col = "blue", lwd = 2)
segments(1:28, out2$summary[25:31,3], 1:28, out2$summary[25:31,28], col = "blue")

plot(1:28, out2$summary[32:38,1], xlab = "Day", ylab = "Detection probability ", las = 1, ylim = c(0, 1), type = "b", col = "blue", pch = 16, frame.plot = FALSE, lwd = 2)
segments(1:28, out2$summary[32:38,3], 1:28, out2$summary[32:38,28], col = "blue")


