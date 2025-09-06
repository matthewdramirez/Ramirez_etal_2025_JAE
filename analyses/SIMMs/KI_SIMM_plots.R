
# R code for Ramirez et al. 2025 "Carbon pathways and trophic attributes are conserved in carnivorous reef fishes across a major human disturbance gradient" Journal of Animal Ecology

# Script to collate SIMM results (Table S2) and plot SIMM results (Figure 3, Figure S2).


##############################

### Load necessary packages
require(ggplot2)
require(tidyr)

### Set your working directory to folder containing this script.
#setwd("C:/Users/.../Cephalopholis argus") # If on a PC
#setwd("/Users/.../Cephalopholis argus") # If on a Mac
setwd("~/Desktop/KI_fish/analyses/SIMMs")

### Plot mean proportional contributions of carbon sources for individual fish (Figure 3; Table S2) ------------

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


# Load the AA d13C dataset
setwd("~/Desktop/KI_fish/data")

csia.N<-read.csv("ki_csia_n.csv") 

# make various variables factors for plotting
csia.N$dist_cat <- factor(csia.N$dist_cat, level = c("Very Low", "Medium", "Very High"))
csia.N$fish_code <- as.factor(csia.N$fish_code)
csia.N$species_code <- as.factor(csia.N$species_code)
#csia.N$sci_name <- factor(csia.N$sci_name, level = c("Aphareus furca", "Caranx melampygus","Cephalopholis argus",  "Lutjanus bohar", "Cephalopholis urodeta", "Lutjanus fulvus"))
csia.N$sci_name <- factor(csia.N$sci_name, level = c("Lutjanus fulvus", "Cephalopholis urodeta", "Lutjanus bohar", "Cephalopholis argus", "Caranx melampygus", "Aphareus furca"))

csia.N$site <- factor(csia.N$site, level = c("10","11","16","15","19","14","8","35","34","27","30"))
# Note site 11 only in AA-C (Zp) data and site 14 NOT in AA-N data but in AA-C (He) and bulk SIA data
# So, AA-N data has VL = 4 sites, M = 3 sites, and VH = 2 sites

csia.N$pub.name <- factor(csia.N$pub.name, level = c("VL6","VL9","VL11","VL1","VL2","M4","M1","M2","M3","VH1","VH3"))
# Note site VL9 only in AA-C (Zp) data and site M4 NOT in AA-N data but in AA-C (He) and bulk SIA data
# So, AA-N data has VL = 4 sites, M = 3 sites, and VH = 2 sites
setwd("~/Desktop/KI_fish/analyses/SIMMs")


# merge datasets

fish <- merge(csia.N, fish, by ="fish_code")

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

# boxplots by individual

ggplot(fish, aes(reorder(fish_code, std_len_mm), ymin = p2.5, lower = p25, middle = p50, upper = p75, ymax = p97.5, fill = source, colour = source)) + 
  geom_boxplot(stat = "identity", width = 0.8) + 
  geom_point(data = fish, aes(reorder(fish_code, std_len_mm),p50), position=position_dodge(width=0.8), colour = "black", size=0.5) + 
  theme_test() + 
  scale_fill_manual(labels = c("Coral", "Detritus", "Turf algae", "Plankton"),
                    values = c("#d01c8b","#a6611a", "#008837", "#2c7bb6")) +
  scale_colour_manual(labels = c("Coral", "Detritus", "Turf algae", "Plankton"),
                    values = c("#d01c8b","#a6611a", "#008837", "#2c7bb6"), guide = "none") +
  theme(axis.text = element_text(colour="black"), text=element_text(size=10), axis.title.x = element_blank(), axis.ticks.length=unit(2.5,"mm")) +
  ylab("Proportion Carbon Source") + #xlab("Species") +
  theme(legend.position = 'top',legend.title.align=0.5) + 
  guides(fill=guide_legend("Carbon Source",nrow=1,title.position = "top",)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(dist_cat ~ sci_name, nrow = 3, scales = "free_x") +
  theme(strip.background = element_blank(), strip.text.x = element_blank())

# Mean species by disturbance category

write.csv(fish,"fish.csv")



### Plot mean proportional contributions of carbon sources for each species (Figure S2) ------------

# load and combine datasets
#fish.g <- read.csv("sum_stat_global.csv") 
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
  theme_test() + 
  scale_fill_manual(labels = c("Coral", "Detritus", "Turf algae", "Plankton"),
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


### Plot mean proportional contributions of carbon sources for each species by SITE (Figure S3) ------------

# make various variables factors for plotting
fish$species_code <- as.factor(fish$species_code)
fish$sci_name <- as.factor(fish$sci_name)
fish$source <- as.factor(fish$source)

ggplot(fish, aes(sci_name, ymin = p2.5, lower = p25, middle = p50, upper = p75, ymax = p97.5, group = site, fill=site)) + 
  geom_boxplot(stat = "identity", width = 0.8) +
  theme_test() + 
  #scale_fill_manual(labels = c("Coral", "Detritus", "Algae", "Plankton"),
  #                                 values = c("#d01c8b","#a6611a", "#008837", "#2c7bb6")) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=14)) + 
  theme(axis.ticks.length=unit(2.5,"mm")) +
  ylab("Proportion Carbon Source") + xlab("Species") +
  theme(legend.position = 'top',legend.title.align=0.5) + 
  guides(fill=guide_legend("Carbon Source",nrow=1,title.position = "top",)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  scale_x_discrete(labels = c("A. furca", "C. melampygus", "C. argus", "L. bohar", "C. urodeta", "L. fulvus")) +
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5), lty="dotted") + 
  facet_wrap(~source)

unique(fish$pub.name)

ggplot(fish, aes(sci_name, Mean, fill = pub.name, colour=pub.name, shape = pub.name)) + 
  geom_errorbar(stat="summary", 
                fun.min=function(x) {mean(x)-1.96* sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+1.96*sd(x)/sqrt(length(x))},
                position = position_dodge(0.75), colour = "black", width = 0, lwd=1) + 
  stat_summary(fun.data = function(x){return(data.frame(y = 1, label = length(x)))}, 
               geom = "text", aes(group=pub.name),
               position = position_dodge(0.75), 
               col = "black") +
  geom_point(stat="summary", fun="mean", position = position_dodge(0.75),  size=5, alpha = 0.75) + 
  scale_fill_manual(values = c(rep("#2A0BD9",4), rep("#ABF8FF",3), rep("#A60021",2)))+
  scale_colour_manual(values = c(rep("#2A0BD9",4), rep("#ABF8FF",3), rep("#A60021",2)))+
  scale_shape_manual(values = c(0,1,2,5,21:23,7,9))+
  #scale_y_continuous(limits=c(3.2,4.5),breaks = seq(3.2,4.4,0.4)) +
  theme_classic()+theme(text = element_text(size=14))+
  #theme(axis.line = element_line(colour = 'black', size = .5))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  theme_test() +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  scale_x_discrete(labels = c("A. furca", "C. melampygus", "C. argus", "L. bohar", "C. urodeta", "L. fulvus")) +
  xlab("Species") + 
  ylab("Proportion Carbon Source") +  
  theme(legend.position = 'top') + 
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black")) +
  guides(fill = guide_legend(title = "Site", title.position = "top", title.hjust = 0.5)) +
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5), lty="dotted") +
  scale_y_continuous(labels = scales::percent, breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  facet_wrap(~source, nrow=4)



