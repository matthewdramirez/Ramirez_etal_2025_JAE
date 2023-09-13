
# Robust trophic pathways in carnivorous reef fishes across a major human disturbance gradient

# Authors: Matthew D. Ramirez [1,2,3], Kelton W. McMahon [2], Neil Rooney [4], Rana W. El-Sabaawi [1], and Julia K. Baum [1]
# Institution: [1] Department of Biology, University of Victoria, PO BOX 1700 Station CSC, Victoria, British Columbia, V8W 2Y2, Canada
# [2] Graduate School of Oceanography, University of Rhode Island, 215 South Ferry Road, Narragansett, Rhode Island, 02882, USA
# [3] Department of Biology and Marine Biology, University of North Carolina Wilmington, 601 South College Road, Wilmington, NC, 28403, USA
# [4] School of Environmental Sciences, University of Guelph, 50 Stone Road East, Guelph, Ontario, N1G 2W1, Canada

# Corresponding Author: Matthew D. Ramirez, Email: ramirezmd@uncw.edu

# Script to create Fig. 1 B-D (multi-panel figure plotting relative biomass and abundance of fish trophic groups at each local human disturbance level)

# Code was modified from Magel et al. 2020. Direct and indirect effects of climate change‐amplified pulse heat stress events on coral reef fish communities. Ecological Applications, 30, e02124.


##############################

### Figure 1 B-D -----------------------------------------------------------

## Load necessary packages

require(tidyverse)
require(tidyr)
require(ggpol)
library(ggplot2)
library(cowplot)
library(plyr)
library(grid)
library(gridExtra)

# Make sure that this contains the "ki_fish_CSIA_C.csv" file
#setwd("C:/Users/...") # If on a PC
#setwd("/Users/...") # If on a Mac
setwd("~/Desktop/KI_fish/data")

# Load the data
fish <- read.csv("ki_fish_data_sum.csv") 

fish$dist_cat <- factor(fish$dist_cat, level = c("Very Low", "Medium", "Very High"))


# Run the 'summary' script (be sure to set file path accordingly)
source("../analyses/summarySE.R")


### > Biomass -----------------------------------------------------------

# Gather and subset data
fish_BM <- gather(fish, "BM_total", "BM_herb","BM_inv","BM_plank","BM_pisc",
                       "BM_omn","BM_gen","BM_coral", "BM_det", key = "species", value = "BM_value")

# Calculate mean and SE values for each dist_cat x species combination
BM_se <- summarySE(fish_BM, measurevar = "BM_value", groupvars = c("dist_cat","species"))

# Sum across all of KI and disturbance levels separately
BM_se$BM_sum <- BM_se$BM_value_mean[9] + BM_se$BM_value_mean[18] + BM_se$BM_value_mean[27]
BM_se$BM_sum_site <- c(rep(BM_se$BM_value_mean[9], 9), rep(BM_se$BM_value_mean[18],9), rep(BM_se$BM_value_mean[27],9))

# Calculate percent of total biomass
BM_se$BM_scaled <- (BM_se$BM_value_mean/BM_se$BM_sum) * 100

# Calculate percent of disturbance level-specific biomass
BM_se$BM_scaled_site <- (BM_se$BM_value_mean/BM_se$BM_sum_site) * 100

# Reorganize dataset
BM_se$species2 <- rep(c("Corallivore","Detritivore","Generalist carnivore","Herbivore","Invertivore","Omnivore","Piscivore","Planktivore","Total"),3)
BM_se$species2 <- factor(BM_se$species2, level = c("Herbivore","Omnivore","Corallivore","Planktivore","Invertivore","Detritivore","Generalist carnivore","Piscivore","Total"))
BM_se$type <- "Biomass"
BM_se <- subset(BM_se, species2 != "Total")
colnames(BM_se)[4:5] <- c("Mean","Median")

df_BM <- tibble(
  Biomass = BM_se$BM_scaled,
  DistLevel = BM_se$dist_cat,
  FunGroup = BM_se$species2
)

print(df_BM, n =24)

# Percent of fish biomass (across all KI) (Fig. 1 B-D)
ggplot(df_BM, aes(x = FunGroup, y = Biomass, fill = DistLevel)) +
  geom_bar(stat = "identity", colour = "black", width = 0.75, alpha = 0.75) +
  facet_wrap(~DistLevel) +
  #facet_share(~DistLevel, dir = "h", scales = "free", reverse_num = TRUE) +
  coord_flip() +
  theme_test() +
  labs(y = "Percent of Total Fish Biomass", x = "Functional Group", title = " ") +
  theme(legend.position = 'top') + 
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  scale_fill_manual(values = c("#2A0BD9","#ABF8FF","#A60021")) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black"))


df_BM_site <- tibble(
  Biomass = BM_se$BM_scaled_site,
  DistLevel = BM_se$dist_cat,
  FunGroup = BM_se$species2
)

print(df_BM_site, n =24)

# Percent of fish biomass (by disturbance level)
ggplot(df_BM_site, aes(x = FunGroup, y = Biomass, fill = DistLevel)) +
  geom_bar(stat = "identity", colour = "black") +
  facet_wrap(~DistLevel) +
  #facet_share(~DistLevel, dir = "h", scales = "free", reverse_num = TRUE) +
  coord_flip() +
  theme_test() +
  labs(y = "Percent of Total Fish Biomass (by Disturbance Level)", x = "Functional Group", title = " ") +
  theme(legend.position = 'top') + 
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  scale_fill_manual(values = c("#2A0BD9","#ABF8FF","#A60021")) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black"))


### > Abundance -----------------------------------------------------------

# Gather and subset data
fish_AB <- gather(fish, "AB_total", "AB_herb","AB_inv","AB_plank","AB_pisc",
                  "AB_omn","AB_gen","AB_coral", "AB_det", key = "species", value = "AB_value")

# Calculate mean and SE values for each dist_cat x species combination
AB_se <- summarySE(fish_AB, measurevar = "AB_value", groupvars = c("dist_cat","species"))

# Sum across all of KI and disturbance levels separately
AB_se$AB_sum <- AB_se$AB_value_mean[9] + AB_se$AB_value_mean[18] + AB_se$AB_value_mean[27]
AB_se$AB_sum_site <- c(rep(AB_se$AB_value_mean[9], 9), rep(AB_se$AB_value_mean[18],9), rep(AB_se$AB_value_mean[27],9))

# Calculate percent of total abundance
AB_se$AB_scaled <- (AB_se$AB_value_mean/AB_se$AB_sum) * 100

# Calculate percent of disturbance level-specific abundance
AB_se$AB_scaled_site <- (AB_se$AB_value_mean/AB_se$AB_sum_site) * 100

# Reorganize dataset
AB_se$species2 <- rep(c("Corallivore","Detritivore","Generalist carnivore","Herbivore","Invertivore","Omnivore","Piscivore","Planktivore","Total"),3)
AB_se$species2 <- factor(AB_se$species2, level = c("Herbivore","Omnivore","Corallivore","Planktivore","Invertivore","Detritivore","Generalist carnivore","Piscivore","Total"))
AB_se$type <- "Abundance"
AB_se <- subset(AB_se, species2 != "Total")
colnames(AB_se)[4:5] <- c("Mean","Median")

df_AB <- tibble(
  Abundance = AB_se$AB_scaled,
  DistLevel = AB_se$dist_cat,
  FunGroup = AB_se$species2
)

print(df, n =24)

# Percent of fish abundance (across KI) (Fig. 1 E-G)
ggplot(df_AB, aes(x = FunGroup, y = Abundance, fill = DistLevel)) +
  geom_bar(stat = "identity", colour = "black", width = 0.75, alpha = 0.75) +
  facet_wrap(~DistLevel) +
  #facet_share(~DistLevel, dir = "h", scales = "free", reverse_num = TRUE) +
  coord_flip() +
  theme_test() +
  labs(y = "Percent of Total Fish Abundance", x = "Functional Group", title = " ") +
  theme(legend.position = 'top') + 
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  scale_fill_manual(values = c("#2A0BD9","#ABF8FF","#A60021")) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black"))


df_AB_site <- tibble(
  Abundance = AB_se$AB_scaled_site,
  DistLevel = AB_se$dist_cat,
  FunGroup = AB_se$species2
)

# Percent of fish biomass (by disturbance level)
ggplot(df_AB_site, aes(x = FunGroup, y = Abundance, fill = DistLevel)) +
  geom_bar(stat = "identity", colour = "black") +
  facet_wrap(~DistLevel) +
  #facet_share(~DistLevel, dir = "h", scales = "free", reverse_num = TRUE) +
  coord_flip() +
  theme_test() +
  labs(y = "Percent of Total Fish Abundance (by Disturbance Level)", x = "Functional Group", title = " ") +
  theme(legend.position = 'top') + 
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  scale_fill_manual(values = c("#2A0BD9","#ABF8FF","#A60021")) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black"))

