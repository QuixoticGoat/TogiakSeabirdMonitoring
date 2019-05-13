
################################################################################
# Creates a .jag file that defined the JAGS model for environmental effects on #
# seabird and nest abundance                                                   #
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov>                                    #  
# Date modified: 6/7/2018                                                      #
################################################################################


# BUGS model file
cat(file = "./OutputData/Jags/mod1.jags", "
model { 
    
    # Priors and constraints
      lambda ~ dgamma(0.001, 0.001)
      p ~ dunif(0,1)

    # Likelihood
      for (i in 1:M) {                 # State model
        N[i] ~ dpois(lambda)
        for (j in 1:J) {
          C[i,j] ~ dbin(p, N[i])      # Observation model
        }
      }
    }
    ")

#-------------------------------------------------------------------------------

# Load packages:
library(tidyverse)
library(jagsUI)


# Load the data:
load("./data/derived/df.list.Rdata")
df.jags <- df.list$df.jags

# Create an wide-formatted input df for jags:
df.jags <- df %>% 
  filter(!is.na(Birds),             # Remove NAs from Birds
         Species == "BLKI") %>%     # Subset only BLKI
  select(-one_of("Nests", "Replicate", "Species")) %>%             # Remove Nests
  arrange(PlotID, Date) %>%         # Reorder by PlotID and Date
  group_by(PlotID, Year) %>%        # Group by PlotID and Date
  mutate(Counter = row_number()) %>%    # Add a column that counts replicates
  select(-Date) %>%                 # Drop the date column
  spread(Counter, Birds)            # Long to wide format
df.jags <- as.data.frame(df.jags)

Year <- rep(1990:2017, nlevels(df.jags$PlotID))
Year <- sort(Year)
PlotID = rep(levels(df.jags$PlotID), length(1990:2017))
joiner <- data.frame(PlotID, Year)
df.jags <- full_join(df.jags, joiner)
rm(Year, PlotID, joiner)

# Replace NA counts with the annual mean (row mean):
foo <- df.jags[, 3:14]  # Select the count columns
k <- which(is.na(foo), arr.ind=TRUE)   # create an index of col/rows with NAs
foo[k] <- rowMeans(foo, na.rm=TRUE)[k[,1]]  # Fill in NAs with row means
foo <- round(foo,0)  # round everything up
df.jags <- cbind(df.jags[, 1:2], foo)  # add the new counts back into the df
df.jags[ is.na(df.jags)] <- NA    # Replace NaNs with NA
df.jags <- df.jags %>%
  arrange(PlotID, Year)
rm(foo)  # Clean up

# Select the Plot 1 data, reps 1:3
df.jags <- subset(df.jags, Year == 2017)
df.jags <- df.jags[complete.cases(df.jags), ]
# df.jags <- df.jags[, 1:5]

# Bundle data
dat <- list(#year = unique(df.jags$Year),               # Years
            #T = length(unique(df.jags$Year)),          # Number of years
            J = ncol(df.jags[3:14]),                     # Number of rep counts/site
            M = length(unique(df.jags$PlotID)),         # Number of sites
            C = as.matrix(df.jags[, 3:14]))              # Counts
                

# Initial values
Nst <- apply(dat$C, 1, max)
inits <- function(){list(N = Nst)}

# Parameters monitored
parameters <- c("lambda", "p", "N")

# MCMC settings
ni <- 15000; nt <- 20; nb <- 5000; nc <- 1

# Call JAGS from R (BRT 3 min)
mod <- jags(dat, inits = inits, parameters, "./OutputData/Jags/mod1.jags", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


# Summarize posteriors
print(mod, digits = 3)
plot(mod)

