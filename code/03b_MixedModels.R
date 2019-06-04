
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

