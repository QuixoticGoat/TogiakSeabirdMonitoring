
################################################################################
# Plots the "Seabird Population and Productivity Monitoring                    #
# at Cape Peirce and Cape Newenham, Alaska" dataset                            #
#                                                                              #
# ServCat ID: 
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov                                     #
# Last modified: 5/31/2018                                                     #
################################################################################

# Required packages:
require(ggplot2)
require(gridExtra)
require(dplyr)

# Load the dataset:
source("./code/01_ImportFormat.R")
df.list <- ImportFormat("./data/raw/SeabirdSurveyData.csv",
                        "./data/derived/dat.Rdata")


#-------------------------------------------------------------------------------
## @knitr Exploratory plots

# Line plots
(p1 <- ggplot(df.list$df.mean, aes(y = Birds, x = Year, group = Species)) +
  geom_line() + 
  geom_errorbar(aes(ymin = Birds-SD.birds, ymax = Birds + SD.birds)) +
  facet_wrap(~ Species, ncol = 1, scales = "free") +
  scale_x_continuous("Year", breaks = seq(1990, 2017, by=5)))

(p2 <- ggplot(df.list$df.mean, aes(y = Nests, x = Year, group = Species)) +
  geom_line() + 
  geom_errorbar(aes(ymin = Nests-SD.nests, ymax = Nests+SD.nests)) +
  facet_wrap(~ Species, ncol = 1, scales = "free") +
  scale_x_continuous("Year", breaks = seq(1990, 2017, by=5)))

(p <- gridExtra::grid.arrange(p1, p2, nrow = 1))
ggsave("./Plots/Trend.pdf", p)

## Boxplots
ggplot(df.list$df, aes(y = Birds, x = as.factor(Year))) +
  geom_boxplot() +
  facet_wrap(~ Species, ncol = 1, scales = "free") +
  labs(x="Year", y="Birds")
  
ggplot(df.list$df, aes(y = Nests, x = as.factor(Year))) +
  geom_boxplot(na.rm=T) +
  facet_wrap(~ Species, ncol = 1, scales = "free") +
  labs(x="Year", y="Nests")



df.list$df %>%
  filter(Year==1990:2018) %>%
ggplot(aes(y = Birds, x = as.factor(lubridate::month(as.Date(Date))))) +
  geom_violin(na.rm=T) +
  facet_wrap(~as.factor(Year)) +
  labs(x="Date", y="Nests")
