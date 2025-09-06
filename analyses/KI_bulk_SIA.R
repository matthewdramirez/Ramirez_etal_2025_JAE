# R code for Ramirez et al. 2025 "Carbon pathways and trophic attributes are conserved in carnivorous reef fishes across a major human disturbance gradient" Journal of Animal Ecology

# Script to characterize variation in reef fish isotopic niche sizes and positions (i.e., Standard Ellipse Areas) via bulk muscle stable carbon (d13C) and nitrogen (d15N) isotope data. 

# Code was modified from: 
# -- supporting resources associated with the article: Jackson et al. 2011. Comparing isotopic niche widths among and within communities: SIBER – Stable Isotope Bayesian Ellipses in R. Journal of Animal Ecology, 80, 595-602.
# ---- including, vignettes provided at: https://github.com/AndrewLJackson/SIBER
# -- supporting resources associated with the article: Swanson et al. 2015. A new probabilistic method for quantifying <i>n</i> -dimensional ecological niches and niche overlap. Ecology, 96, 318-324.
# ---- including, vignettes provided at:https://github.com/mlysy/nicheROVER

# All figures were modified in Adobe Illustrator prior to publication.

##############################

## Load necessary packages
require(ggplot2)
require(SIBER)
require(dplyr)
require(tidyr)
require(ggpubr)
require(nicheROVER)

## Set your working directory
# Make sure that this contains the "ki_bulk_sia.csv" file
#setwd("C:/Users/...") # If on a PC
#setwd("/Users/...") # If on a Mac
setwd("~/Desktop/KI_fish/data")


## Load the data
fish.bulk <- read.csv("ki_bulk_sia.csv")

fish.bulk$sci_name <- as.factor(fish.bulk$sci_name)
fish.bulk$site <- as.factor(fish.bulk$site)
fish.bulk$dist_cat <- factor(fish.bulk$dist_cat, level = c("Very Low", "Medium", "Very High"))

fish.bulk$species <- factor(fish.bulk$species, level = c("LU.FULV", "CE.UROD", "LU.BOHA","CE.ARGU","CA.MELA","AP.FURC"))

### Bulk SIA baseline corrections  -----------------------------------------------------------
fish.C.AA <- read.csv("ki_csia_c.csv")
fish.N.AA <- read.csv("ki_csia_n.csv")


## AA-N correction factors

fish.N.AA$fish_code <- as.factor(fish.N.AA$fish_code)
fish.N.AA$sci_name <- as.factor(fish.N.AA$sci_name)

fish.N.AA$dist_cat <- as.factor(fish.N.AA$dist_cat)
fish.N.AA$dist_cat <- factor(fish.N.AA$dist_cat, level = c("Very Low", "Medium", "Very High"))

fish.N.AA$site <- as.factor(fish.N.AA$site)
fish.N.AA$site <- factor(fish.N.AA$site, level = c("19","15","16","11","10","14","8","35","34", "27","30"))

fish.N.AA.Sr <- gather(fish.N.AA, 'Phe', 'Lys', key = "Sr.N", value = "d15N")

#Lutjanus fulvus
source.N.LF <- subset(fish.N.AA.Sr, species_code == "LF") %>% 
    group_by(sci_name, site) %>% summarize(n = n(), 
              mean.N = round(mean(d15N, na.rm = TRUE),1),
              sd.N = round(sd(d15N, na.rm = TRUE),1))
source.N.LF$N_dev_from_mean <- round((mean(source.N.LF$mean.N) - source.N.LF$mean.N),2)

#Cephalolphois urodeta
source.N.CU <- subset(fish.N.AA.Sr, species_code == "CU") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.N = round(mean(d15N, na.rm = TRUE),1),
                                         sd.N = round(sd(d15N, na.rm = TRUE),1))
source.N.CU$N_dev_from_mean <- round((mean(source.N.CU$mean.N) - source.N.CU$mean.N),2)

#Lutjanus bohar
source.N.LB <- subset(fish.N.AA.Sr, species_code == "LB") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.N = round(mean(d15N, na.rm = TRUE),1),
                                         sd.N = round(sd(d15N, na.rm = TRUE),1))
source.N.LB$N_dev_from_mean <- round((mean(source.N.LB$mean.N) - source.N.LB$mean.N),2)

#Cephalolphois argus
source.N.CA <- subset(fish.N.AA.Sr, species_code == "CA") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.N = round(mean(d15N, na.rm = TRUE),1),
                                         sd.N = round(sd(d15N, na.rm = TRUE),1))
source.N.CA$N_dev_from_mean <- round((mean(source.N.CA$mean.N) - source.N.CA$mean.N),2)

#Caranx melampygus
source.N.CM <- subset(fish.N.AA.Sr, species_code == "CM") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.N = round(mean(d15N, na.rm = TRUE),1),
                                         sd.N = round(sd(d15N, na.rm = TRUE),1))
source.N.CM$N_dev_from_mean <- round((mean(source.N.CM$mean.N) - source.N.CM$mean.N),2)

#Aphareus furca
source.N.AF <- subset(fish.N.AA.Sr, species_code == "AF") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.N = round(mean(d15N, na.rm = TRUE),1),
                                         sd.N = round(sd(d15N, na.rm = TRUE),1))
source.N.AF$N_dev_from_mean <- round((mean(source.N.AF$mean.N) - source.N.AF$mean.N),2)

#merge tibbles

source.N <- rbind(source.N.LF,source.N.CU,source.N.LB,source.N.CA,source.N.CM,source.N.AF)

source.N <- as.data.frame(source.N)

# Add correction factors for missing sites, based on closest adjacent sire
source.N[nrow(source.N) + 1,] = c("Cephalopholis argus",  14,  NA, NA, NA, 0.54) # sub site 8 for 14
source.N[nrow(source.N) + 1,] = c("Cephalopholis argus",  30,  NA, NA, NA, -0.46) # sub site 27 for 30
source.N[nrow(source.N) + 1,] = c("Cephalopholis urodeta",  34,  NA, NA, NA, 0.22) # sub site 35 for 34
source.N[nrow(source.N) + 1,] = c("Lutjanus fulvus",  14,  NA, NA, NA, 0.57) # sub site 8 for 14



## AA-C correction factors
fish.C.AA$fish_code <- as.factor(fish.C.AA$fish_code)
fish.C.AA$sci_name <- as.factor(fish.C.AA$sci_name)

fish.C.AA$dist_cat <- as.factor(fish.C.AA$dist_cat)
fish.C.AA$dist_cat <- factor(fish.C.AA$dist_cat, level = c("Very Low", "Medium", "Very High"))

fish.C.AA$site <- as.factor(fish.C.AA$site)
fish.C.AA$site <- factor(fish.C.AA$site, level = c("19","15","16","11","10","14","8","35","34", "27","30"))

fish.C.EAA <- gather(fish.C.AA, 'Leu', 'Ile','Lys','Phe','Thr','Val', key = "EAA", value = "d13C")

#Lutjanus fulvus
fish.C.EAA.LF <- subset(fish.C.EAA, species_code == "LF") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.C = round(mean(d13C, na.rm = TRUE),1),
                                         sd.C = round(sd(d13C, na.rm = TRUE),1))
fish.C.EAA.LF$C_dev_from_mean <- round((mean(fish.C.EAA.LF$mean.C) - fish.C.EAA.LF$mean.C),2)

#Cephalolphois urodeta
fish.C.EAA.CU <- subset(fish.C.EAA, species_code == "CU") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.C = round(mean(d13C, na.rm = TRUE),1),
                                         sd.C = round(sd(d13C, na.rm = TRUE),1))
fish.C.EAA.CU$C_dev_from_mean <- round((mean(fish.C.EAA.CU$mean.C) - fish.C.EAA.CU$mean.C),2)

#Lutjanus bohar
fish.C.EAA.LB <- subset(fish.C.EAA, species_code == "LB") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.C = round(mean(d13C, na.rm = TRUE),1),
                                         sd.C = round(sd(d13C, na.rm = TRUE),1))
fish.C.EAA.LB$C_dev_from_mean <- round((mean(fish.C.EAA.LB$mean.C) - fish.C.EAA.LB$mean.C),2)

#Cephalolphois argus
fish.C.EAA.CA <- subset(fish.C.EAA, species_code == "CA") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.C = round(mean(d13C, na.rm = TRUE),1),
                                         sd.C = round(sd(d13C, na.rm = TRUE),1))
fish.C.EAA.CA$C_dev_from_mean <- round((mean(fish.C.EAA.CA$mean.C) - fish.C.EAA.CA$mean.C),2)

#Caranx melampygus
fish.C.EAA.CM <- subset(fish.C.EAA, species_code == "CM") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.C = round(mean(d13C, na.rm = TRUE),1),
                                         sd.C = round(sd(d13C, na.rm = TRUE),1))
fish.C.EAA.CM$C_dev_from_mean <- round((mean(fish.C.EAA.CM$mean.C) - fish.C.EAA.CM$mean.C),2)

#Aphareus furca
fish.C.EAA.AF <- subset(fish.C.EAA, species_code == "AF") %>% 
  group_by(sci_name, site) %>% summarize(n = n(), 
                                         mean.C = round(mean(d13C, na.rm = TRUE),1),
                                         sd.C = round(sd(d13C, na.rm = TRUE),1))
fish.C.EAA.AF$C_dev_from_mean <- round((mean(fish.C.EAA.AF$mean.C) - fish.C.EAA.AF$mean.C),2)

#merge tibbles

fish.C.EAA <- rbind(fish.C.EAA.LF,fish.C.EAA.CU,fish.C.EAA.LB,fish.C.EAA.CA,fish.C.EAA.CM,fish.C.EAA.AF)

fish.C.EAA <- as.data.frame(fish.C.EAA)

# Add correction factors for missing sites, based on closest adjacent sire
fish.C.EAA[nrow(fish.C.EAA) + 1,] = c("Cephalopholis argus",  14,  NA, NA, NA, 0.84) # sub site 8 for 14
fish.C.EAA[nrow(fish.C.EAA) + 1,] = c("Cephalopholis argus",  30,  NA, NA, NA, -0.36) # sub site 27 for 30
fish.C.EAA[nrow(fish.C.EAA) + 1,] = c("Cephalopholis urodeta",  34,  NA, NA, NA, -0.38) # sub site 35 for 34
fish.C.EAA[nrow(fish.C.EAA) + 1,] = c("Lutjanus fulvus",  14,  NA, NA, NA, -1.13) # sub site 8 for 14


# Merge correction factors with bulk SIA dataset
fish.bulk <- left_join(fish.bulk, source.N, by=c("sci_name","site")) %>%
  select(-n, -mean.N, -sd.N)

fish.bulk <- left_join(fish.bulk, fish.C.EAA, by=c("sci_name","site")) %>%
  select(-n, -mean.C, -sd.C)

fish.bulk$N_dev_from_mean <- as.numeric(fish.bulk$N_dev_from_mean)
fish.bulk$C_dev_from_mean <- as.numeric(fish.bulk$C_dev_from_mean)


# Created corrected d15N and d13C columns
fish.bulk$d15Ncor <- fish.bulk$d15N + fish.bulk$N_dev_from_mean
fish.bulk$d13Ccor <- fish.bulk$d13C + fish.bulk$C_dev_from_mean


### Bulk SIA summary statistics for sample set (Table S1) -----------------------------------------------------------

# Species
fish.bulk %>% 
  group_by(species) %>%
  summarize(n = n(), meanSL = mean(standard_length_.mm.)/10,
            minSL = min(standard_length_.mm.)/10,
            maxSL = max(standard_length_.mm.)/10)

# Species x Disturbance Level
fish.bulk %>% 
  group_by(species, dist_cat) %>%
  summarize(n = n(), meanSL = mean(standard_length_.mm.)/10,
            minSL = min(standard_length_.mm.)/10,
            maxSL = max(standard_length_.mm.)/10)

##############################

### Plot SIA versus size (Figure S5) -----------------------------------------------------------

new_labels <- c("AP.FURC" = "A. furca", 
                "CA.MELA" = "C. melampygus", 
                "CE.ARGU" = "C. argus", 
                "LU.BOHA" = "L. bohar", 
                "CE.UROD" = "C. urodeta", 
                "LU.FULV" = "L. fulvus")


ggplot(fish.bulk, aes(standard_length_.mm./10, d13Ccor)) +
  geom_point(aes(fill=dist_cat), pch = 21, col = "black", size = 3, alpha = 0.75) + 
  scale_fill_manual(values = c("#2A0BD9", "#ABF8FF", "#A60021"))+
  geom_smooth(method="gam", colour="black", alpha=0.25)+
  #stat_regline_equation(label.y = 4.6, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = -7, aes(label = after_stat(rr.label)), size = 3.25) + 
  theme_test() +
  #theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_y_continuous(breaks=c(seq(-19,-7,2)),limits=c(-19,-7)) +
  theme(legend.position = 'top') + 
  xlab("Standard Length (cm)") +
  ylab(expression({delta}^13*C~('\u2030'))) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black")) +
  theme(strip.text = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  facet_wrap(~species, scales = "free_x", nrow=1, labeller=labeller(species = new_labels)) 


ggplot(fish.bulk, aes(standard_length_.mm./10, d15Ncor)) +
  geom_point(aes(fill=dist_cat), pch = 21, col = "black", size = 3, alpha = 0.75) + 
  scale_fill_manual(values = c("#2A0BD9", "#ABF8FF", "#A60021"))+
  geom_smooth(method="gam", colour="black", alpha=0.25)+
  #stat_regline_equation(label.y = 4.6, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 16, aes(label = after_stat(rr.label)), size = 3.25) + 
  theme_test() +
  #theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_y_continuous(breaks=c(seq(8,16,2)),limits=c(8,16)) +
  theme(legend.position = 'top') + 
  xlab("Standard Length (cm)") +
  ylab(expression({delta}^15*N~('\u2030'))) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black")) +
  theme(strip.text = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  facet_wrap(~species, scales = "free_x", nrow=1, labeller=labeller(species = new_labels)) 


### Reshape dataset for SIBER analyses -----------------------------------------------------------

# Note the following numeric values will replace letter names as required for the SIBER analyses

fish.bulk.SEA <- fish.bulk
fish.bulk.SEA$species <- factor(fish.bulk.SEA$species, level = c("LU.FULV", "CE.UROD", "LU.BOHA","CE.ARGU","CA.MELA","AP.FURC"))
levels(fish.bulk.SEA$species) <- c("1", "2", "3","4","5","6")

fish.bulk.SEA$dist_cat <- as.character(fish.bulk.SEA$dist_cat)
fish.bulk.SEA$species <- as.character(fish.bulk.SEA$species)

unique(fish.bulk.SEA$species)

# Assign numeric values to Disturbance Levels
fish.bulk.SEA$dist_cat[fish.bulk.SEA$dist_cat == "Very Low"] <- 1
fish.bulk.SEA$dist_cat[fish.bulk.SEA$dist_cat == "Medium"] <- 2
fish.bulk.SEA$dist_cat[fish.bulk.SEA$dist_cat == "Very High"] <- 3

# Assign numeric values to Species
fish.bulk.SEA$species[fish.bulk.SEA$species == "AP.FURC"] <- 6
fish.bulk.SEA$species[fish.bulk.SEA$species == "CA.MELA"] <- 5
fish.bulk.SEA$species[fish.bulk.SEA$species == "CE.ARGU"] <- 4
fish.bulk.SEA$species[fish.bulk.SEA$species == "LU.BOHA"] <- 3
fish.bulk.SEA$species[fish.bulk.SEA$species == "CE.UROD"] <- 2
fish.bulk.SEA$species[fish.bulk.SEA$species == "LU.FULV"] <- 1

# Subset and reorganize data for Species, Disturbance Level, d13C values, and d15N values
fish.bulk.SEA <- fish.bulk.SEA[,c(2,7,17,16)]
#fish.bulk.SEA <- fish.bulk.SEA[,c(2,7,11,12)]
names(fish.bulk.SEA)[1:4] <- c("community","group","iso1","iso2") #community is species #group is disturbance level
fish.bulk.SEA <- fish.bulk.SEA[,c(3,4,2,1)]
fish.bulk.SEA <- fish.bulk.SEA[with(fish.bulk.SEA, order(community, group)),]

# Create separate datasets for each species
fish.bulk.SEA.AF <- subset(fish.bulk.SEA, community == "6")
fish.bulk.SEA.CM <- subset(fish.bulk.SEA, community == "5")
fish.bulk.SEA.CA <- subset(fish.bulk.SEA, community == "4")
fish.bulk.SEA.LB <- subset(fish.bulk.SEA, community == "3")
fish.bulk.SEA.CU <- subset(fish.bulk.SEA, community == "2")
fish.bulk.SEA.LF <- subset(fish.bulk.SEA, community == "1")



### Univariate isotope comparisons among disturbance levels -----------------------------------------------------------

# Perform Kruskal-wallis tests to evaluate differences in univariate d13C and d15N values among disturbance levels 
# 'group' is disturbance level

kruskal.test(iso1 ~ group, data = fish.bulk.SEA.AF) #p-value = 0.1138
kruskal.test(iso2 ~ group, data = fish.bulk.SEA.AF) #p-value = 0.0049  #DIFFERENCE

kruskal.test(iso1 ~ group, data = fish.bulk.SEA.CM) #p-value = 0.0202  #DIFFERENCE
kruskal.test(iso2 ~ group, data = fish.bulk.SEA.CM) #p-value = 0.1163

kruskal.test(iso1 ~ group, data = fish.bulk.SEA.CA) #p-value = 0.4654
kruskal.test(iso2 ~ group, data = fish.bulk.SEA.CA) #p-value = 0.6959

kruskal.test(iso1 ~ group, data = fish.bulk.SEA.LB) #p-value = 0.6324
kruskal.test(iso2 ~ group, data = fish.bulk.SEA.LB) #p-value = 0.5851

kruskal.test(iso1 ~ group, data = fish.bulk.SEA.CU) #p-value = 0.0560
kruskal.test(iso2 ~ group, data = fish.bulk.SEA.CU) #p-value = 0.0013  #DIFFERENCE

kruskal.test(iso1 ~ group, data = fish.bulk.SEA.LF) #p-value = 0.0165  #DIFFERENCE
kruskal.test(iso2 ~ group, data = fish.bulk.SEA.LF) #p-value = 0.1500  


# Perform pairwise comparisons using Wilcoxon rank sum exact test with bonferroni correction for multiple comparisons

pairwise.wilcox.test(fish.bulk.SEA.AF$iso2, fish.bulk.SEA.AF$group, p.adjust.method = "bonferroni", paired = FALSE) #1 (Very Low) v 2 (Medium) key difference

pairwise.wilcox.test(fish.bulk.SEA.CM$iso2, fish.bulk.SEA.CM$group, p.adjust.method = "bonferroni", paired = FALSE) # No differences

pairwise.wilcox.test(fish.bulk.SEA.CU$iso2, fish.bulk.SEA.CU$group, p.adjust.method = "bonferroni", paired = FALSE) #1 (Very Low) v 2 (Medium) AND 1 (Very Low) vs 3 (Very High) key differences

pairwise.wilcox.test(fish.bulk.SEA.LF$iso1, fish.bulk.SEA.LF$group, p.adjust.method = "bonferroni", paired = FALSE) #1 (Very Low) v 2 (Medium) key difference


### Quantify isotopic niche sizes using SIBER -----------------------------------------------------------

# create SIBER objects for each species (used for plotting + niche overlap analysis)
siber.AF <- createSiberObject(fish.bulk.SEA.AF)
siber.CM <- createSiberObject(fish.bulk.SEA.CM)
siber.CA <- createSiberObject(fish.bulk.SEA.CA)
siber.LB <- createSiberObject(fish.bulk.SEA.LB)
siber.CU <- createSiberObject(fish.bulk.SEA.CU)
siber.LF <- createSiberObject(fish.bulk.SEA.LF)

# Full SIBER model (used for all other analyses)
siber.full <- createSiberObject(fish.bulk.SEA)


# > Biplots w/ SEAc (Figure 4) -----------------------------------------------------------

# Create d13C-d15N biplots for each species. Add sample size-correct maximum likelihood standard ellipse areas (SEAc).

# Set color palette
palette(c("#2A0BD9","#ABF8FF","#A60021"))

par(mfrow=c(2,3),mar=c(5,5,1,1))


min(fish.bulk$d13Ccor)
max(fish.bulk$d13Ccor)

min(fish.bulk$d15Ncor)
max(fish.bulk$d15Ncor)

# Create lists of plotting arguments to be passed onwards to each of the three plotting functions.
#community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
#group.ellipses.args  <- list(n = 100, p.interval = NULL, lty = 1, lwd = 2)
#group.hull.args      <- list(lty = 3, col = "black")

# Species-specific biplots with maximum likelihood standard ellipses (SEAc), which account for c. 40% of the data. 
# Standard ellipses were manually removed for disturbance levels with N < 3.
LF <- {
  plotSiberObject(siber.LF,
                  ax.pad = 2, 
                  hulls = F, community.hulls.args, 
                  ellipses = F, group.ellipses.args,
                  group.hulls = F, group.hull.args,
                  bty = "o",
                  iso.order = c(1,2),
                  xlab = expression({delta}^13*C~('\u2030')),
                  ylab = expression({delta}^15*N~('\u2030')),
                  y.limits = c(9,16),
                  x.limits = c(-18,-7))
  points(fish.bulk.SEA.LF$iso1, fish.bulk.SEA.LF$iso2, bg=fish.bulk.SEA.LF$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.LF, m= c(24,14,11), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-18, 16, expression(paste(bold("(A)"), italic(" L. fulvus"), " (n = 49)")), adj = c(0,1))
  #legend("bottomleft",c("Very High", "Medium", "Very Low"), pch=c(21,21,21), 
  #       col = c("black","black","black"), pt.cex=1,
  #       pt.bg = alpha(c("#A60021","#ABF8FF","#2A0BD9"),0.75))
}

CU <- {
  plotSiberObject(siber.CU,
                  ax.pad = 2, 
                  hulls = F, community.hulls.args, 
                  ellipses = F, group.ellipses.args,
                  group.hulls = F, group.hull.args,
                  bty = "o",
                  iso.order = c(1,2),
                  xlab = expression({delta}^13*C~('\u2030')),
                  ylab = expression({delta}^15*N~('\u2030')),
                  y.limits = c(9,16),
                  x.limits = c(-18,-13))
  points(fish.bulk.SEA.CU$iso1, fish.bulk.SEA.CU$iso2, bg=fish.bulk.SEA.CU$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  text(-18, 16, expression(paste(bold("(B)"), italic(" C. urodeta"), " (n = 64)")), adj = c(0,1))
  plotGroupEllipses(siber.CU, m= c(14,4,6), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  #legend("bottomright",c("Very High", "Medium", "Very Low"), pch=c(21,21,21), 
  #       col = c("black","black","black"), pt.cex=1,
  #       pt.bg = alpha(c("#A60021","#ABF8FF","#2A0BD9"),0.75))
}


LB <- {
  plotSiberObject(siber.LB,
                  ax.pad = 2, 
                  hulls = F, community.hulls.args, 
                  ellipses = F, group.ellipses.args,
                  group.hulls = F, group.hull.args,
                  bty = "o",
                  iso.order = c(1,2),
                  xlab = expression({delta}^13*C~('\u2030')),
                  ylab = expression({delta}^15*N~('\u2030')),
                  y.limits = c(9,16),
                  x.limits = c(-18,-13))
  points(fish.bulk.SEA.LB$iso1, fish.bulk.SEA.LB$iso2, bg=fish.bulk.SEA.LB$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.LB, m= c(20,21,23), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-18, 16, expression(paste(bold("(C)"), italic(" L. bohar"), " (n = 24)")), adj = c(0,1))
  #legend("bottomright",c("Very High", "Medium", "Very Low"), pch=c(21,21,21), 
  #       col = c("black","black","black"), pt.cex=1,
  #       pt.bg = alpha(c("#A60021","#ABF8FF","#2A0BD9"),0.75))
}


CA <- {
  plotSiberObject(siber.CA,
                  ax.pad = 2, 
                  hulls = F, community.hulls.args, 
                  ellipses = F, group.ellipses.args,
                  group.hulls = F, group.hull.args,
                  bty = "o",
                  iso.order = c(1,2),
                  xlab = expression({delta}^13*C~('\u2030')),
                  ylab = expression({delta}^15*N~('\u2030')),
                  y.limits = c(9,16),
                  x.limits = c(-18,-13))
  points(fish.bulk.SEA.CA$iso1, fish.bulk.SEA.CA$iso2, bg=fish.bulk.SEA.CA$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.CA, m= c(29,28,11), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-18, 16, expression(paste(bold("(D)"), italic(" C. argus"), " (n = 68)")), adj = c(0,1))
  #legend("bottomright",c("Very Low", "Medium", "Very High"), pch=c(21,21,21), 
  #     col = c("black","black","black"), pt.cex=1,
  #     pt.bg = alpha(c("#2A0BD9","#ABF8FF","#A60021"),0.75))
}


CM <- {
  plotSiberObject(siber.CM,
                  ax.pad = 2, 
                  hulls = F, community.hulls.args, 
                  ellipses = F, group.ellipses.args,
                  group.hulls = F, group.hull.args,
                  bty = "o",
                  iso.order = c(1,2),
                  xlab = expression({delta}^13*C~('\u2030')),
                  ylab = expression({delta}^15*N~('\u2030')),
                  y.limits = c(9,16),
                  x.limits = c(-18,-13))
  points(fish.bulk.SEA.CM$iso1, fish.bulk.SEA.CM$iso2, bg=fish.bulk.SEA.CM$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.CM, m= c(6,3,2), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-18, 16, expression(paste(bold("(E)"), italic(" C.  melampygus"), " (n = 11)")), adj = c(0,1))
  #legend("bottomright",c("Very Low", "Medium", "Very High"), pch=c(21,21,21), 
  #       col = c("black","black","black"), pt.cex=1,
  #       pt.bg = alpha(c("#2A0BD9","#ABF8FF","#A60021"),0.75))
}


AF <- {
  plotSiberObject(siber.AF,
                  ax.pad = 2, 
                  hulls = F, community.hulls.args, 
                  ellipses = F, group.ellipses.args,
                  group.hulls = F, group.hull.args,
                  bty = "o",
                  iso.order = c(1,2),
                  xlab = expression({delta}^13*C~('\u2030')),
                  ylab = expression({delta}^15*N~('\u2030')),
                  y.limits = c(9,16),
                  x.limits = c(-18,-13))
  points(fish.bulk.SEA.AF$iso1, fish.bulk.SEA.AF$iso2, bg=fish.bulk.SEA.AF$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.AF, m= c(9,4,3), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-18, 16, expression(paste(bold("(F)"), italic(" A. furca"), " (n = 16)")), adj = c(0,1))
  legend("bottomright",c("Very Low", "Medium", "Very High"), pch=c(21,21,21), 
         col = c("black","black","black"), pt.cex=1,
         pt.bg = alpha(c("#2A0BD9","#ABF8FF","#A60021"),0.75))
}


# Calculate summary statistics for each group: TA, SEA and SEAc
# TA: convex hull area
# SEA: maximum likelihood standard ellipse area
# SEAc: sample size-corrected maximum likelihood standard ellipse area

group.ML <- groupMetricsML(siber.full)
print(group.ML)


# > Fit Bayesian models to the data -----------------------------------------------------------

# The initial step for comparing compare isotopic niche width among groups is to fit Bayesian multivariate normal distributions to each group in the dataset. 

# set parameters defining how the sampling algorithm is to run
parms <- list()
parms$n.iter <- 2 * 10^4  # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3  # discard the first set of values
parms$n.thin <- 10  # thin the posterior by this many
parms$n.chains <- 2  # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.full, parms, priors)


# calculate the SEA on the posterior distribution of covariance matrix for each group (i.e, Bayesian SEA or SEA-B) (Table S3).
SEA.B <- siberEllipses(ellipses.posterior)

# calculate some credible intervals 
cr.p <- c(0.95,0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
(SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p,2))

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
(SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T))


# > Plot SEA-B data (Figure S4) -------------------------------------------------------
palette(c("#2A0BD9","#ABF8FF","#A60021"))

par(mfrow=c(1,1),mar=c(5,5,1,1))

siberDensityPlot(SEA.B, xticklabels = c(rep(c("VL","M", "VH"),6)), 
                 xlab = c("Disturbance Level"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 ct="median",
                 bty = "o",
                 las = 1, 
                 probs = c(95, 75, 50),
                 ylim = c(0,20))

# add red colored points for SEAc values calculated earlier
points(1:ncol(SEA.B), group.ML[3,], col=c("#2A0BD9","#ABF8FF","#A60021"), pch = c(16,16,16), lwd = 2)

# add vertical lines between species and species labels 
abline(v = c(3.5,6.5,9.5,12.5,15.5), col="black", lwd=1, lty=2)

text(2, 20, expression(italic("L. fulvus")))
text(5, 20, expression(italic("C. urodeta")))
text(8, 20, expression(italic("L. bohar")))
text(11, 20, expression(italic("C. argus")))
text(14, 20, expression(italic("C. melampygus")))
text(17, 20, expression(italic("A. furca")))



### Compare isotopic niche size (Table S4) -------------------------------------------------------

# In order to test whether one group's ellipse is smaller or larger than another, we can simply calculate the probability 
# that its posterior distribution is smaller (or larger). This is achieved by comparing each pair of posterior draws for 
# both groups, and determining which is smaller in magnitude. We then find the proportion of draws that are smaller, and 
# this is a direct proxy for the probability that one group's posterior distribution (of ellipse size in this case) is 
# smaller than the other.

# calculate the proportion, and hence probability, of the SEA.B for group 1 being smaller than the SEA.B for group 2 (e.g., AF.1.2)
# calculate the proportion, and hence probability, of the SEA.B for group 1 being smaller than the SEA.B for group 3 (e.g., AF.1.3)
# calculate the proportion, and hence probability, of the SEA.B for group 2 being smaller than the SEA.B for group 2 (e.g., AF.2.3)

## L. fulvus
(LF.1.2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B))
(LF.1.3 <- sum( SEA.B[,1] < SEA.B[,3] ) / nrow(SEA.B))
(LF.2.3 <- sum( SEA.B[,2] < SEA.B[,3] ) / nrow(SEA.B))

## C. urodeta
(CU.1.2 <- sum( SEA.B[,4] < SEA.B[,5] ) / nrow(SEA.B))
(CU.1.3 <- sum( SEA.B[,4] < SEA.B[,6] ) / nrow(SEA.B))
(CU.2.3 <- sum( SEA.B[,5] < SEA.B[,6] ) / nrow(SEA.B))

## L. bohar
(LB.1.2 <- sum( SEA.B[,7] < SEA.B[,8] ) / nrow(SEA.B))
(LB.1.3 <- sum( SEA.B[,7] < SEA.B[,9] ) / nrow(SEA.B))
(LB.2.3 <- sum( SEA.B[,8] < SEA.B[,9] ) / nrow(SEA.B))

## C. argus
(CA.1.2 <- sum( SEA.B[,10] < SEA.B[,11] ) / nrow(SEA.B))
(CA.1.3 <- sum( SEA.B[,10] < SEA.B[,12] ) / nrow(SEA.B))
(CA.2.3 <- sum( SEA.B[,11] < SEA.B[,12] ) / nrow(SEA.B))

## C. melampygus
(CM.1.2 <- sum( SEA.B[,13] < SEA.B[,14] ) / nrow(SEA.B))

## A. furca
(AF.1.2 <- sum( SEA.B[,16] < SEA.B[,17] ) / nrow(SEA.B))
(AF.1.3 <- sum( SEA.B[,16] < SEA.B[,18] ) / nrow(SEA.B))
(AF.2.3 <- sum( SEA.B[,17] < SEA.B[,18] ) / nrow(SEA.B))




### Quantify isotopic niche overlap (Table S5) -------------------------------------------------------

# Overlap calculation uses nsamples = nprob = 10000 (1e4) for higher accuracy.
clrs <- c("#2A0BD9","#ABF8FF","#A60021") # colors for each species
nsamples <- 1e4


## A. furcus
fish.AF.par <- tapply(1:nrow(fish.bulk.AF), fish.bulk.AF$group,
                      function(ii) niw.post(nsamples = nsamples, X = fish.bulk.AF[ii,1:2]))
AF.over.stat <- overlap(fish.AF.par, nreps = nsamples, nprob = nsamples, alpha = .95)

#The mean AF.overlap metrics calculated across iterations for both niche region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
AF.over.mean <- apply(AF.over.stat, c(1:2), mean)*100
round(AF.over.mean, 1)
AF.over.cred <- apply(AF.over.stat*100, c(1:2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(AF.over.cred,2) # display alpha = .95 niche region

#In the returned plot, Species $A$ is along the rows and Species $B$ is along columns. The plots represent the posterior probability that an individual from the species indicated by the row will be found within the niche of the species indicated by the column header. Before you plot, you must decide upon your $\alpha$-level, and make sure the variable `AF.over.stat` reflects this choice of $\alpha$.
# AF.overlap plot.Before you run this, make sure that you have chosen your alpha level.
overlap.plot(AF.over.stat, col = clrs, mean.cred.col = "black", equal.axis = TRUE,
             xlab = expression(paste(italic("Aphareus furcus"))),
             species.names = c("VL", "M", "VH"))



## C. melampygus
# VH Could not be evaluated because must have more observations than niche dimensions when prior parameter Psi is missing.

fish.bulk.CM.2 <-subset(fish.bulk.CM, group != 3) # remove group with N < 3

fish.CM.par <- tapply(1:nrow(fish.bulk.CM.2), fish.bulk.CM.2$group,
                      function(ii) niw.post(nsamples = nsamples, X = fish.bulk.CM.2[ii,1:2]))
CM.over.stat <- overlap(fish.CM.par, nreps = nsamples, nprob = 1e4, alpha = .95)

#The mean CM.overlap metrics CMlculated across iterations for both niche region sizes (alpha = .95 and alpha = .99) CMn be CMlculated and displayed in an array.
CM.over.mean <- apply(CM.over.stat, c(1:2), mean)*100
round(CM.over.mean, 1)
CM.over.cred <- apply(CM.over.stat*100, c(1:2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(CM.over.cred,2) # display alpha = .95 niche region

#In the returned plot, Species $A$ is along the rows and Species $B$ is along columns. The plots represent the posterior probability that an individual from the species indiCMted by the row will be found within the niche of the species indiCMted by the column header. Before you plot, you must decide upon your $\alpha$-level, and make sure the variable `CM.over.stat` reflects this choice of $\alpha$.
# CM.overlap plot.Before you run this, make sure that you have chosen your alpha level.
overlap.plot(CM.over.stat, col = clrs, mean.cred.col = "black", equal.axis = TRUE,
             xlab = expression(paste(italic("Caranx melampygus"))),
             species.names = c("VL", "M"))



## C. argus
fish.CA.par <- tapply(1:nrow(fish.bulk.CA), fish.bulk.CA$group,
                      function(ii) niw.post(nsamples = nsamples, X = fish.bulk.CA[ii,1:2]))
CA.over.stat <- overlap(fish.CA.par, nreps = nsamples, nprob = nsamples, alpha = .95)

#The mean CA.overlap metrics calculated across iterations for both niche region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
CA.over.mean <- apply(CA.over.stat, c(1:2), mean)*100
round(CA.over.mean, 1)
CA.over.cred <- apply(CA.over.stat*100, c(1:2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(CA.over.cred,2) # display alpha = .95 niche region

#In the returned plot, Species $A$ is along the rows and Species $B$ is along columns. The plots represent the posterior probability that an individual from the species indicated by the row will be found within the niche of the species indicated by the column header. Before you plot, you must decide upon your $\alpha$-level, and make sure the variable `CA.over.stat` reflects this choice of $\alpha$.
# CA.overlap plot.Before you run this, make sure that you have chosen your alpha level.
overlap.plot(CA.over.stat, col = clrs, mean.cred.col = "black", equal.axis = TRUE,
             xlab = expression(paste(italic("Cephalopholis argus"))),
             species.names = c("VL", "M", "VH"))


## L. bohar
fish.LB.par <- tapply(1:nrow(fish.bulk.LB), fish.bulk.LB$group,
                      function(ii) niw.post(nsamples = nsamples, X = fish.bulk.LB[ii,1:2]))
LB.over.stat <- overlap(fish.LB.par, nreps = nsamples, nprob = nsamples, alpha = .95)

#The mean LB.overlap metrics LBlculated across iterations for both niche region sizes (alpha = .95 and alpha = .99) LBn be LBlculated and displayed in an array.
LB.over.mean <- apply(LB.over.stat, c(1:2), mean)*100
round(LB.over.mean, 1)
LB.over.cred <- apply(LB.over.stat*100, c(1:2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(LB.over.cred,2) # display alpha = .95 niche region

#In the returned plot, Species $A$ is along the rows and Species $B$ is along columns. The plots represent the posterior probability that an individual from the species indiLBted by the row will be found within the niche of the species indiLBted by the column header. Before you plot, you must decide upon your $\alpha$-level, and make sure the variable `LB.over.stat` reflects this choice of $\alpha$.
# LB.overlap plot.Before you run this, make sure that you have chosen your alpha level.
overlap.plot(LB.over.stat, col = clrs, mean.cred.col = "black", equal.axis = TRUE,
             xlab = expression(paste(italic("Lutjanus bohar"))),
             species.names = c("VL", "M", "VH"))



## C. urodeta
fish.CU.par <- tapply(1:nrow(fish.bulk.CU), fish.bulk.CU$group,
                      function(ii) niw.post(nsamples = nsamples, X = fish.bulk.CU[ii,1:2]))
CU.over.stat <- overlap(fish.CU.par, nreps = nsamples, nprob = nsamples, alpha = .95)

#The mean CU.overlap metrics CUlculated across iterations for both niche region sizes (alpha = .95 and alpha = .99) CUn be CUlculated and displayed in an array.
CU.over.mean <- apply(CU.over.stat, c(1:2), mean)*100
round(CU.over.mean, 1)
CU.over.cred <- apply(CU.over.stat*100, c(1:2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(CU.over.cred,2) # display alpha = .95 niche region

#In the returned plot, Species $A$ is along the rows and Species $B$ is along columns. The plots represent the posterior probability that an individual from the species indiCUted by the row will be found within the niche of the species indiCUted by the column header. Before you plot, you must decide upon your $\alpha$-level, and make sure the variable `CU.over.stat` reflects this choice of $\alpha$.
# CU.overlap plot.Before you run this, make sure that you have chosen your alpha level.
overlap.plot(CU.over.stat, col = clrs, mean.cred.col = "black", equal.axis = TRUE,
             xlab = expression(paste(italic("Cephalopholis urodeta"))),
             species.names = c("VL", "M", "VH"))



## L. fulvus
fish.LF.par <- tapply(1:nrow(fish.bulk.LF), fish.bulk.LF$group,
                   function(ii) niw.post(nsamples = nsamples, X = fish.bulk.LF[ii,1:2]))
LF.over.stat <- overlap(fish.LF.par, nreps = nsamples, nprob = nsamples, alpha = .95)

#The mean LF.overlap metrics calculated across iterations for both niche region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
LF.over.mean <- apply(LF.over.stat, c(1:2), mean)*100
round(LF.over.mean, 1)
LF.over.cred <- apply(LF.over.stat*100, c(1:2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(LF.over.cred,2) # display alpha = .95 niche region

#In the returned plot, Species $A$ is along the rows and Species $B$ is along columns. The plots represent the posterior probability that an individual from the species indicated by the row will be found within the niche of the species indicated by the column header. Before you plot, you must decide upon your $\alpha$-level, and make sure the variable `LF.over.stat` reflects this choice of $\alpha$.
# LF.overlap plot.Before you run this, make sure that you have chosen your alpha level.
overlap.plot(LF.over.stat, col = clrs, mean.cred.col = "black", equal.axis = TRUE,
             xlab = expression(paste(italic("Lutjanus fulvus"))),
             species.names = c("VL", "M", "VH"))



