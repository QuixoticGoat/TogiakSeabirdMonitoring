
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
df.list <- ImportFormat("./data/raw/seabird_data.csv",
                        "./data/derived/df.list.Rdata")


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
ggsave("./output/plots/annual_trend_lineSE.pdf", p)

## Boxplots
(p1 <- ggplot(df.list$df, aes(y = Birds, x = as.factor(Year))) +
  geom_boxplot() +
  facet_wrap(~ Species, ncol = 1, scales = "free") +
  scale_x_continuous("Year", breaks = seq(1990, 2017, by=5)) +
  labs(x="Year", y="Birds"))
p2 <- ggplot(df.list$df, aes(y = Nests, x = as.factor(Year))) +
  geom_boxplot(na.rm=T) +
  facet_wrap(~ Species, ncol = 1, scales = "free") +
  labs(x="Year", y="Nests")
(p <- gridExtra::grid.arrange(p1, p2, nrow = 1))
ggsave("./output/plots/annual_trend_boxplot.pdf", p)

#-------------------------------------------------------------------------------

## Plot total annual bird abundances and nests by Julian date

# Total counts:
df.list$df %>%
  group_by(Julian, Year) %>%
  summarise(Birds = sum(Birds, na.rm=T),
            Nests = sum(Nests, na.rm=T)) %>%
  ggplot() +
  #geom_point(aes(y=Nests, x=Julian), col="blue") +
  geom_point(aes(y=Birds, x=Julian)) +
  stat_smooth(aes(y=Birds, x=Julian), method="lm") +
  labs(title="Total Counts", x="Julian day") +
  scale_x_continuous(breaks=seq(150, 200, 25)) +
  facet_wrap(~Year, scales="free_y")
ggsave("./output/plots/julian_counts_total.pdf")

# Total nests:
df.list$df %>%
  group_by(Julian, Year) %>%
  summarise(Birds = sum(Birds, na.rm=T),
            Nests = sum(Nests, na.rm=T)) %>%
  ggplot() +
  #geom_point(aes(y=Nests, x=Julian), col="blue") +
  geom_point(aes(y=Nests, x=Julian)) +
  stat_smooth(aes(y=Nests, x=Julian), method="lm") +
  labs(title="Total Nests", x="Julian day") +
  scale_x_continuous(breaks=seq(150, 200, 25)) +
  facet_wrap(~Year, scales="free_y")
ggsave("./output/plots/julian_nests_total.pdf")

# BLKI counts:
foo <- df.list$df %>%
  group_by(Julian, Year, Species) %>%
  summarise(Birds = sum(Birds, na.rm=T),
            Nests = sum(Nests, na.rm=T)) %>%
  subset(Species=="BLKI") %>%
  ggplot() +
  geom_point(aes(y=Birds, x=Julian)) +
  stat_smooth(aes(y=Birds, x=Julian), method="lm") +
  labs(title="Black-legged Kittiwake Counts", x="Julian day") +
  scale_x_continuous(breaks=seq(150, 200, 25)) +
  facet_wrap(~Year, scales="free_y")
ggsave("./output/plots/julian_counts_blki.pdf")

# COMU counts:
df.list$df %>%
  group_by(Julian, Year, Species) %>%
  summarise(Birds = sum(Birds, na.rm=T),
            Nests = sum(Nests, na.rm=T)) %>%
  subset(Species=="COMU") %>%
  ggplot() +
  geom_point(aes(y=Birds, x=Julian)) +
  stat_smooth(aes(y=Birds, x=Julian), method="lm") +
  labs(title="Common Murre Counts", x="Julian day") +
  scale_x_continuous(breaks=seq(150, 200, 25)) +
  facet_wrap(~Year, scales="free_y")
ggsave("./output/plots/julian_counts_comu.pdf")

# PECO counts:
df.list$df %>%
  group_by(Julian, Year, Species) %>%
  summarise(Birds = sum(Birds, na.rm=T),
            Nests = sum(Nests, na.rm=T)) %>%
  subset(Species=="PECO") %>%
  ggplot() +
  geom_point(aes(y=Birds, x=Julian)) +
  stat_smooth(aes(y=Birds, x=Julian), method="lm") +
  labs(title="Pacific Cormorant Counts", x="Julian day") +
  scale_x_continuous(breaks=seq(150, 200, 25)) +
  facet_wrap(~Year, scales="free_y")
ggsave("./output/plots/julian_counts_peco.pdf")

# BLKI nests:
df.list$df %>%
  group_by(Julian, Year, Species) %>%
  summarise(Birds = sum(Birds, na.rm=T),
            Nests = sum(Nests, na.rm=T)) %>%
  subset(Species=="BLKI") %>%
  ggplot() +
  geom_point(aes(y=Nests, x=Julian)) +
  stat_smooth(aes(y=Nests, x=Julian), method="lm") +
  labs(title="Black-legged Kittiwake Nests", x="Julian day") +
  scale_x_continuous(breaks=seq(150, 200, 25)) +
  facet_wrap(~Year, scales="free_y")
ggsave("./output/plots/julian_nests_blki.pdf")

# PECO nests:
df.list$df %>%
  group_by(Julian, Year, Species) %>%
  summarise(Birds = sum(Birds, na.rm=T),
            Nests = sum(Nests, na.rm=T)) %>%
  subset(Species=="PECO") %>%
  ggplot() +
  geom_point(aes(y=Nests, x=Julian)) +
  stat_smooth(aes(y=Nests, x=Julian), method="lm") +
  labs(title="Pacific Cormorant Nests", x="Julian day") +
  scale_x_continuous(breaks=seq(150, 200, 25)) +
  facet_wrap(~Year, scales="free_y")
ggsave("./output/plots/julian_nests_peco.pdf")

#-------------------------------------------------------------------------------
## Scatterplots of Birds~Nest


foo <- df.list$df %>%
  group_by(Julian, Year, Species) %>%
  summarise(Birds = sum(Birds, na.rm=T),
            Nests = sum(Nests, na.rm=T)) %>%
  subset(Species != "COMU")

foo %>%
ggplot(aes(x=Nests, y=Birds)) +
  geom_point() +
  stat_smooth(method="lm", se=FALSE) +
  facet_wrap(~Species, scales="free")
ggsave("./output/plots/scatter_birds_nests.pdf")
