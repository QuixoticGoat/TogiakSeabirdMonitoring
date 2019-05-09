

.libPaths("C:/Users/jlaufenberg/Documents/R/win-library/3.3")


set.seed(12345)
nplots <- 30
nyears <- 20
nreps <- 12
mu.lambda <- runif(nplots, 0.95, 1.05)
mu.loglambda <- log(mu.lambda)
sd.loglambda <- 0.025
eps.loglambda <- matrix(rnorm(nplots * (nyears-1), 0, sd.loglambda), nplots, nyears-1)
lambda <- exp(mu.loglambda + eps.loglambda)
N <- matrix(0, nplots, nyears)
N[,1] <- round(runif(nplots, 200, 400))
for(t in 2:nyears){
    N[,t] <- rpois(nplots, N[,t-1] * lambda[,t-1])
}
p <- runif(nyears, 0.5, 0.9)
C <- array(0, c(nplots, nreps, nyears))
for(i in 1:nplots){
    for(t in 1:nyears){
        C[i,,t] <- rbinom(nreps, N[i,t], p[t])
    }
}

## Bundle data
maxCyear1 <- apply(C,c(1,3),max)[,1]
psi <- matrix(0,600,nplots)
for(i in 1:nplots){
    psi[1:(maxCyear1[i]-1),i] <- 0
    psi[maxCyear1[i]:600,i] <- 1/(601-maxCyear1[i])
}
data <- list(C=C, P=nplots, J=nreps, T=nyears, psi=psi)

## Initial values
Ninit <- N
Ninit[,2:20] <- NA
inits <- function() list(mu.loglambda=mu.loglambda,#runif(30,-0.1,0.1),
                         sd.loglambda=sd.loglambda,#runif(1,0,1),
                         N=Ninit,#matrix(c(round(runif(30,maxCyear1,600)),rep(NA,30*19)),30,20),
                         p=p#runif(20,0.5,1)
                         )

## Parameters monitored
pars <- c("N", "mu.loglambda", "sd.loglambda", "mu.lambda", "lambda", "p")


library(rjags)


cat("   ", format(Sys.time()), "\n\n")
system.time({
        jm <- jags.model(file="model.jags", data=data, inits=inits, n.chains=1, n.adapt=500)
        jc <- coda.samples(jm, pars, n.iter=1000, thin=1)
})


save.image(file = "RawMCMCoutput/CMR.TNRB.Model1.ZINITtest.RData")
plot(TNRBout1.0,ask=TRUE)
TNRB.model1.0.mc <- mcmc.list(sapply(TNRBout1.0, function(x) x))






### JAGS model
sink("model.jags")
cat("
model {
    ## Priors
    sd.loglambda ~ dunif(0,10)
    tau.loglambda <- pow(sd.loglambda,-2)
    for(i in 1:P){
        mu.loglambda[i] ~ dnorm(0,0.01)
        N[i,1] ~ dcat(psi[,i])
    } #i
    for(t in 1:T){
        p[t] ~ dunif(0,1)
    } #t
    ## Population process model
    for(i in 1:P){
        for(t in 2:T){
            loglambda[i,t-1] <- mu.loglambda[i] + eps.loglambda[i,t-1]
            eps.loglambda[i,t-1] ~ dnorm(0,tau.loglambda)
            lambda[i,t-1] <- exp(loglambda[i,t-1])
            N[i,t] ~ dpois(N[i,t-1] * lambda[i,t-1])
        } #t
    } #i
    ## Observation process model
    for(i in 1:P){
        for(t in 1:T){
            for(j in 1:J){
                C[i,j,t] ~ dbin(p[t],N[i,t])
            } #j
        } #t
    } #i
    ## Derived parameters
    for(i in 1:P){
        mu.lambda[i] <- exp(mu.loglambda[i])
    } #i
}
",fill=TRUE)
sink()





