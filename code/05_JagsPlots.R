
################################################################################
# Plots the posterior estimates of the JAGS model for environmental effects on #
# seabird and nest abundance                                                   #
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov>                                    #  
# Date modified: 6/7/2018                                                      #
################################################################################


library(ggplot2)

abundance <- mod$mean$mean.abundance
lci <- mod$q2.5$mean.abundance
uci <- mod$q97.5$mean.abundance
df.plot <- data.frame(abundance, lci, uci, "plotID" = unique(df.jags$PlotID))
rm(abundance, lci, uci)

ggplot(df.plot, aes(y = abundance, x = plotID)) +
  geom_point() +
  geom_errorbar(data=df.plot, aes(ymin = lci, ymax = uci)) +
  labs(title = "1991", x = "Est. abundance", y = "Plot")
ggsave("./Plots/Plot1990.pdf")
