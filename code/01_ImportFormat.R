
################################################################################
# Imports and formats the "Seabird Population and Productivity Monitoring      #
# at Cape Peirce and Cape Newenham, Alaska" csv dataset                        #
#                                                                              #
# inpath: directory path of the input csv dataset                              #
# outpath: directory path to save a list containing the formatted dataframe    #
# (df) and the formatted dataframe for jags (df.jags).                         #
#                                                                              #
# ServCat ID:                                                                  #
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov>                                    #
# Last modified: 05-09-2019                                                    #                                                                       #
################################################################################


ImportFormat <- function(inpath, outpath) {
  
  #Required packages
  library(tidyverse)
  
  df <- read.csv(inpath)
  
  # Reformats the data:
  df$Date <- as.Date(df$Date, format = "%m/%d/%Y")
  df$Year <- as.numeric(df$Year)
  df$Replicate <- as.factor(df$Replicate)
  df$Birds <- as.numeric(df$Birds)
  df$Nests <- as.numeric(df$Nests)
  df$Yearn <- as.numeric(as.factor(df$Year))

  # Creates a formatted dataframe for jags:
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
  PlotID <-  rep(levels(df.jags$PlotID), length(1990:2017))
  joiner <- data.frame(PlotID, Year)
  df.jags <- full_join(df.jags, joiner)
  rm(Year, PlotID, joiner)
  
  # Replace NA counts with the annual mean (row mean):
  foo <- df.jags[, 4:15]  # Select the count columns
  k <- which(is.na(foo), arr.ind=TRUE)   # create an index of col/rows with NAs
  foo[k] <- rowMeans(foo, na.rm=TRUE)[k[,1]]  # Fill in NAs with row means
  foo <- round(foo,0)  # round everything up
  df.jags <- cbind(df.jags[, 1:3], foo)  # add the new counts back into the df
  df.jags[ is.na(df.jags)] <- NA    # Replace NaNs with NA
  df.jags <- df.jags %>%
    arrange(PlotID, Year)
  rm(foo, k)  # Clean up
  
  # Removes plots that were not surveyed often:
  selected <- c("10BC", "10D", "10E", "10F")
  df.jags <-  df.jags[!df.jags$PlotID %in% selected, ]
  df.jags <- droplevels(df.jags)
  rm(selected)
  
  # Removes plots that were pretty much all zeros:
  selected <- c("2", "3C", "5A", "7C", "7F", "7G")
  df.jags <-  df.jags[!df.jags$PlotID %in% selected, ]
  df.jags <- droplevels(df.jags)
  rm(selected)
  
  # Creates a dataframe of annual mean values for plotting:
  df.mean <- df %>%
    group_by(Year, Species, Replicate) %>%
    summarise(B = sum(Birds, na.rm = T), 
              N = sum(Nests, na.rm = T)) %>%
    summarise(Birds = mean(B, na.rm = T),
              SD.birds = sd(B, na.rm = T),
              Nests = mean(N, na.rm = T),
              SD.nests = sd(N, na.rm = T))
  df.mean[df.mean == 0] <- NA
  df.mean$Year <- as.numeric(as.character(df.mean$Year))
  
  # Combined dataframes into a list:
  df.list <- list()
  df.list$df <- df
  df.list$df.jags <- df.jags
  df.list$df.mean <- df.mean
  
  # Saves the list:
  save(df.list, file=outpath)

  return(df.list)
}
