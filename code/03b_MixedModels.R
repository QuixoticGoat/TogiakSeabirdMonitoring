
################################################################################
# Regression the "Seabird Population and Productivity Monitoring               #
# at Cape Peirce and Cape Newenham, Alaska" dataset                            #
#                                                                              #
# ServCat ID: 
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov                                     #
# Last modified: 5/31/2018                                                     #            
################################################################################


#-------------------------------------------------------------------------------
##---- Generalized linear mixed effects model 

require(lme4)
require(tidyverse)

load("./data/derived/df.list.Rdata")
df <- df.list$df %>%
  filter(Species=="BLKI",
         Year==2016) %>%
  droplevels()

mod <- list()

mod[1] <- glmer(Birds ~ Yearn + (1|PlotID), data = df, family = poisson(link="log"))  # Null



#-------------------------------------------------------------------------------
##---- Open n-mixture model

library(unmarked)

# Format the data for unmarked:
df.jags <- df.list$df.jags %>%
  filter(Year==1990) %>%
  droplevels()

df <- df.list$df

selected <- c("10BC", "10D", "10E", "10F", "2", "3C", "5A", "7C", "7F", "7G")
df <-  df[!df$PlotID %in% selected, ]
df <- droplevels(df)
rm(selected)

df <- df %>%
  group_by(Year) %>%
  mutate(JulianS = Julian - min(Julian) + 1)  # Add a column that visit day of that year


df.periods <- df %>% 
  filter(!is.na(Birds),             # Remove NAs from Birds
         Species == "BLKI") %>%     # Subset only BLKI
  select(-one_of("Nests", "Replicate", "Species", "Birds", "Julian")) %>%  # Remove Nests, Replicate, Species and Julian columns
  arrange(PlotID, Date) %>%         # Reorder by PlotID and Date
  group_by(PlotID, Year) %>%        # Group by PlotID and Date
  mutate(Counter = row_number()) %>%    # Add a column that counts replicates
  select(-Date) %>%                 # Drop the date column
  spread(Counter, JulianS) %>%            # Long to wide format
  arrange(Year) %>%                    # Reorder by Year
  filter(Year==1990) %>%
  droplevels()
  
# Counts (12 visits, 30 sites)
y1 <- data.matrix(df.jags[,4:12])

# Max. number of primary periods (visits):
np <- 9

# Matrix of primary period values:
primaryPeriod1 <- data.matrix(df.periods[,4:12])

# Create an unmarked object:
umf1 <- unmarkedFramePCO(y=y1, numPrimary=np, primaryPeriod=primaryPeriod1)


