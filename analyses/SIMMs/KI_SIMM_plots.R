
# Robust trophic pathways in carnivorous reef fishes across a major human disturbance gradient

# Authors: Matthew D. Ramirez [1,2,3], Kelton W. McMahon [2], Neil Rooney [4], Rana W. El-Sabaawi [1], and Julia K. Baum [1]
# Institution: [1] Department of Biology, University of Victoria, PO BOX 1700 Station CSC, Victoria, British Columbia, V8W 2Y2, Canada
# [2] Graduate School of Oceanography, University of Rhode Island, 215 South Ferry Road, Narragansett, Rhode Island, 02882, USA
# [3] Department of Biology and Marine Biology, University of North Carolina Wilmington, 601 South College Road, Wilmington, NC, 28403, USA
# [4] School of Environmental Sciences, University of Guelph, 50 Stone Road East, Guelph, Ontario, N1G 2W1, Canada

# Corresponding Author: Matthew D. Ramirez, Email: ramirezmd@uncw.edu

# Script to collate SIMM results (Supplementary Table 1) and plot SIMM results (Fig. 2, Supplementary Fig. 1).


##############################

### Load necessary packages
require(ggplot2)
require(tidyr)

### Set your working directory to folder containing this script.
#setwd("C:/Users/.../Cephalopholis argus") # If on a PC
#setwd("/Users/.../Cephalopholis argus") # If on a Mac
setwd("~/Desktop/KI_fish/analyses/SIMMs")

### Plot mean proportional contributions of carbon sources for individual fish (Fig. 2) ------------

# load and combine datasets
fish.AF <- read.csv("sum_stat_AF.csv") 
fish.CM <- read.csv("sum_stat_CM.csv") 
fish.CA <- read.csv("sum_stat_CA.csv") 
fish.CU <- read.csv("sum_stat_CU.csv") 
fish.LB <- read.csv("sum_stat_LB.csv") 
fish.LF <- read.csv("sum_stat_LF.csv") 

fish <- rbind(fish.AF,
              fish.CM,
              fish.CA,
              fish.CU,
              fish.LB,
              fish.LF)

# make various variables factors for plotting
fish$dist_cat <- factor(fish$dist_cat, level = c("Very Low", "Medium", "Very High"))
fish$fish_code <- as.factor(fish$fish_code)
fish$species_code <- as.factor(fish$species_code)
fish$sci_name <- factor(fish$sci_name, level = c("Aphareus furca", "Caranx melampygus","Cephalopholis argus",  "Lutjanus bohar", "Cephalopholis urodeta", "Lutjanus fulvus"))
fish$source <- as.factor(fish$source)
#fish$group <- as.factor(fish$group)


# for plotting purposes, adjust std_len_mm by 0.01 or 0.02 for individuals of same species that have exact match 
# AF1 = AF7 = AF9
# CA12 = CA13
# CA2 = CA10
# LB 6 = LB8

fish[fish$fish_code == "AF7",21] <- fish[fish$fish_code == "AF7",21] + 0.01
fish[fish$fish_code == "AF9",21] <- fish[fish$fish_code == "AF9",21] + 0.02
fish[fish$fish_code == "CA13",21] <- fish[fish$fish_code == "CA13",21] + 0.01
fish[fish$fish_code == "CA10",21] <- fish[fish$fish_code == "CA10",21] + 0.01
fish[fish$fish_code == "LB8",21] <- fish[fish$fish_code == "LB8",21] + 0.01


ggplot(fish, aes(reorder(fish_code, std_len_mm), Mean, fill=source)) + 
  geom_bar(position="fill",stat = "identity",colour="black", width = 1, alpha=0.9) +
  theme_test() + 
  scale_fill_manual(values=c("#d01c8b","#a6611a", "#008837", "#2c7bb6")) +
  #theme(axis.text.x = element_text(angle = 45, hjust=1)) +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  #theme(legend.title = element_blank()) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=14)) + 
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(plot.title.position = 'plot') +
  theme(legend.position = 'top') + 
  ylab("Proportion Carbon Source") + xlab("Fish ID") +
  #theme(legend.key.size = unit(0.4,"cm"), legend.position = "top",legend.text=element_text(size=12),legend.title =element_text(size=14)) +
  guides(fill=guide_legend("Carbon Source",nrow=1)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0,0.25,0.5,0.75,1)) +
  facet_wrap(sci_name ~ dist_cat, ncol = 3, scales = "free_x")

### Plot mean proportional contributions of carbon sources for each species (Supplementary Fig. 1) ------------

# load and combine datasets
fish.AF.g <- read.csv("sum_stat_AF_global.csv") 
fish.CM.g <- read.csv("sum_stat_CM_global.csv") 
fish.CA.g <- read.csv("sum_stat_CA_global.csv") 
fish.CU.g <- read.csv("sum_stat_CU_global.csv") 
fish.LB.g <- read.csv("sum_stat_LB_global.csv") 
fish.LF.g <- read.csv("sum_stat_LF_global.csv") 

fish.g <- rbind(fish.AF.g,
              fish.CM.g,
              fish.CA.g,
              fish.CU.g,
              fish.LB.g,
              fish.LF.g)

# make various variables factors for plotting
fish.g$species_code <- as.factor(fish.g$species_code)
fish.g$sci_name <- as.factor(fish.g$sci_name)
fish.g$source <- as.factor(fish.g$source)

ggplot(fish.g, aes(sci_name, ymin = p2.5, lower = p25, middle = p50, upper = p75, ymax = p97.5, fill = source)) + 
  geom_boxplot(stat = "identity", width = 0.8) +
  theme_test() + scale_fill_manual(labels = c("Coral", "Detritus", "Algae", "Plankton"),
                  values = c("#d01c8b","#a6611a", "#008837", "#2c7bb6")) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=14)) + 
  theme(axis.ticks.length=unit(2.5,"mm")) +
  ylab("Proportion Carbon Source") + xlab("Species") +
  theme(legend.position = 'top',legend.title.align=0.5) + 
  guides(fill=guide_legend("Carbon Source",nrow=1,title.position = "top",)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  scale_x_discrete(labels = c("A. furca", "C. melampygus", "C. argus", "L. bohar", "C. urodeta", "L. fulvus")) +
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5), lty="dotted")






