
# Robust trophic pathways in carnivorous reef fishes across a major human disturbance gradient

# Authors: Matthew D. Ramirez [1,2,3], Kelton W. McMahon [2], Neil Rooney [4], Rana W. El-Sabaawi [1], and Julia K. Baum [1]
# Institution: [1] Department of Biology, University of Victoria, PO BOX 1700 Station CSC, Victoria, British Columbia, V8W 2Y2, Canada
# [2] Graduate School of Oceanography, University of Rhode Island, 215 South Ferry Road, Narragansett, Rhode Island, 02882, USA
# [3] Department of Biology and Marine Biology, University of North Carolina Wilmington, 601 South College Road, Wilmington, NC, 28403, USA
# [4] School of Environmental Sciences, University of Guelph, 50 Stone Road East, Guelph, Ontario, N1G 2W1, Canada

# Corresponding Author: Matthew D. Ramirez, Email: ramirezmd@uncw.edu

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
require(ggpubr)
require(nicheROVER)

## Set your working directory
# Make sure that this contains the "ki_bulk_sia.csv" file
#setwd("C:/Users/...") # If on a PC
#setwd("/Users/...") # If on a Mac
setwd("~/Desktop/KI_fish/data")


## Load the data
fish.bulk <- read.csv("ki_bulk_sia.csv")


## Acquire summary statistics for sample set
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

### Plot SIA versus size (Extended Data Fig. 3) -----------------------------------------------------------

new_labels <- c("AP.FURC" = "A. furca", 
                "CA.MELA" = "C. melampygus", 
                "CE.ARGU" = "C. argus", 
                "LU.BOHA" = "L. bohar", 
                "CE.UROD" = "C. urodeta", 
                "LU.FULV" = "L. fulvus")


ggplot(fish.bulk, aes(standard_length_.mm./10, d13C)) +
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


ggplot(fish.bulk, aes(standard_length_.mm./10, d15N)) +
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

# Assign numeric values to Disturbance Levels
fish.bulk$dist_cat[fish.bulk$dist_cat == "Very Low"] <- 1
fish.bulk$dist_cat[fish.bulk$dist_cat == "Medium"] <- 2
fish.bulk$dist_cat[fish.bulk$dist_cat == "Very High"] <- 3

# Assign numeric values to Species
fish.bulk$species[fish.bulk$species == "AP.FURC"] <- 1
fish.bulk$species[fish.bulk$species == "CA.MELA"] <- 2
fish.bulk$species[fish.bulk$species == "CE.ARGU"] <- 3
fish.bulk$species[fish.bulk$species == "LU.BOHA"] <- 4
fish.bulk$species[fish.bulk$species == "CE.UROD"] <- 5
fish.bulk$species[fish.bulk$species == "LU.FULV"] <- 6

# Subset and reorganize data for Species, Disturbance Level, d13C values, and d15N values
fish.bulk <- fish.bulk[,c(2,7,11,12)]
names(fish.bulk)[1:4] <- c("community","group","iso1","iso2") #community is species #group is disturbance level
fish.bulk <- fish.bulk[,c(3,4,2,1)]
fish.bulk <- fish.bulk[with(fish.bulk, order(community, group)),]

# Create separate datasets for each species
fish.bulk.AF <- subset(fish.bulk, community == "1")
fish.bulk.CM <- subset(fish.bulk, community == "2")
fish.bulk.CA <- subset(fish.bulk, community == "3")
fish.bulk.LB <- subset(fish.bulk, community == "4")
fish.bulk.CU <- subset(fish.bulk, community == "5")
fish.bulk.LF <- subset(fish.bulk, community == "6")



### Univariate isotope comparisons among disturbance levels -----------------------------------------------------------

# Perform Kruskal-wallis tests to evaluate differences in univariate d13C and d15N values among disturbance levels 
# 'group' is disturbance level

kruskal.test(iso1 ~ group, data = fish.bulk.AF) #p-value = 0.1138
kruskal.test(iso2 ~ group, data = fish.bulk.AF) #p-value = 0.0049  #DIFFERENCE

kruskal.test(iso1 ~ group, data = fish.bulk.CM) #p-value = 0.0202  #DIFFERENCE
kruskal.test(iso2 ~ group, data = fish.bulk.CM) #p-value = 0.1163

kruskal.test(iso1 ~ group, data = fish.bulk.CA) #p-value = 0.4654
kruskal.test(iso2 ~ group, data = fish.bulk.CA) #p-value = 0.6959

kruskal.test(iso1 ~ group, data = fish.bulk.LB) #p-value = 0.6324
kruskal.test(iso2 ~ group, data = fish.bulk.LB) #p-value = 0.5851

kruskal.test(iso1 ~ group, data = fish.bulk.CU) #p-value = 0.0560
kruskal.test(iso2 ~ group, data = fish.bulk.CU) #p-value = 0.0013  #DIFFERENCE

kruskal.test(iso1 ~ group, data = fish.bulk.LF) #p-value = 0.0165  #DIFFERENCE
kruskal.test(iso2 ~ group, data = fish.bulk.LF) #p-value = 0.1500  


# Perform pairwise comparisons using Wilcoxon rank sum exact test with bonferroni correction for multiple comparisons

pairwise.wilcox.test(fish.bulk.AF$iso2, fish.bulk.AF$group, p.adjust.method = "bonferroni", paired = FALSE) #1 (Very Low) v 2 (Medium) key difference

pairwise.wilcox.test(fish.bulk.CM$iso2, fish.bulk.CM$group, p.adjust.method = "bonferroni", paired = FALSE) # No differences

pairwise.wilcox.test(fish.bulk.CU$iso2, fish.bulk.CU$group, p.adjust.method = "bonferroni", paired = FALSE) #1 (Very Low) v 2 (Medium) AND 1 (Very Low) vs 3 (Very High) key differences

pairwise.wilcox.test(fish.bulk.LF$iso1, fish.bulk.LF$group, p.adjust.method = "bonferroni", paired = FALSE) #1 (Very Low) v 2 (Medium) key difference


### Quantify isotopic niche sizes using SIBER -----------------------------------------------------------

# create SIBER objects for each species (used for plotting + niche overlap analysis)
siber.AF <- createSiberObject(fish.bulk.AF)
siber.CM <- createSiberObject(fish.bulk.CM)
siber.CA <- createSiberObject(fish.bulk.CA)
siber.LB <- createSiberObject(fish.bulk.LB)
siber.CU <- createSiberObject(fish.bulk.CU)
siber.LF <- createSiberObject(fish.bulk.LF)

# Full SIBER model (used for all other analyses)
siber.full <- createSiberObject(fish.bulk)


# > Biplots w/ SEAc (Fig. 4) -----------------------------------------------------------

# Create d13C-d15N biplots for each species. Add sample size-correct maximum likelihood standard ellipse areas (SEAc).

# Set color palette
palette(c("#2A0BD9","#ABF8FF","#A60021"))

par(mfrow=c(2,3),mar=c(5,5,1,1))

# Create lists of plotting arguments to be passed onwards to each of the three plotting functions.
#community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
#group.ellipses.args  <- list(n = 100, p.interval = NULL, lty = 1, lwd = 2)
#group.hull.args      <- list(lty = 3, col = "black")

# Species-specific biplots with maximum likelihood standard ellipses (SEAc), which account for c. 40% of the data. 
# Standard ellipses were manually removed for disturbance levels with N < 3.
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
                  y.limits = c(8,17),
                  x.limits = c(-19,-7))
  points(fish.bulk.AF$iso1, fish.bulk.AF$iso2, bg=fish.bulk.AF$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.AF, m= c(9,4,3), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-19, 17, expression(paste(bold("(A)"), italic(" Aphareus furca"))), adj = c(0,1))
  legend("bottomright",c("Very Low", "Medium", "Very High"), pch=c(21,21,21), 
         col = c("black","black","black"), pt.cex=1,
         pt.bg = alpha(c("#2A0BD9","#ABF8FF","#A60021"),0.75))
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
                  y.limits = c(8,17),
                  x.limits = c(-19,-7))
  points(fish.bulk.CM$iso1, fish.bulk.CM$iso2, bg=fish.bulk.CM$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.CM, m= c(6,3,2), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-19, 17, expression(paste(bold("(B)"), italic(" Caranx melampygus"))), adj = c(0,1))
  #legend("bottomright",c("Very Low", "Medium", "Very High"), pch=c(21,21,21), 
  #       col = c("black","black","black"), pt.cex=1,
  #       pt.bg = alpha(c("#2A0BD9","#ABF8FF","#A60021"),0.75))
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
                  y.limits = c(8,17),
                  x.limits = c(-19,-7))
  points(fish.bulk.CA$iso1, fish.bulk.CA$iso2, bg=fish.bulk.CA$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.CA, m= c(29,28,11), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-19, 17, expression(paste(bold("(C)"), italic(" Cephalopholis argus"))), adj = c(0,1))
  #legend("bottomright",c("Very Low", "Medium", "Very High"), pch=c(21,21,21), 
  #     col = c("black","black","black"), pt.cex=1,
  #     pt.bg = alpha(c("#2A0BD9","#ABF8FF","#A60021"),0.75))
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
                  y.limits = c(8,17),
                  x.limits = c(-19,-7))
  points(fish.bulk.LB$iso1, fish.bulk.LB$iso2, bg=fish.bulk.LB$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.LB, m= c(20,21,23), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-19, 17, expression(paste(bold("(D)"), italic(" Lutjanus bohar"))), adj = c(0,1))
  #legend("bottomright",c("Very High", "Medium", "Very Low"), pch=c(21,21,21), 
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
                  y.limits = c(8,17),
                  x.limits = c(-19,-7))
  points(fish.bulk.CU$iso1, fish.bulk.CU$iso2, bg=fish.bulk.CU$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  text(-19, 17, expression(paste(bold("(E)"), italic(" Cephalopholis urodeta"))), adj = c(0,1))
  plotGroupEllipses(siber.CU, m= c(14,4,6), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  #legend("bottomright",c("Very High", "Medium", "Very Low"), pch=c(21,21,21), 
  #       col = c("black","black","black"), pt.cex=1,
  #       pt.bg = alpha(c("#A60021","#ABF8FF","#2A0BD9"),0.75))
}

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
                  y.limits = c(8,17),
                  x.limits = c(-19,-7))
  points(fish.bulk.LF$iso1, fish.bulk.LF$iso2, bg=fish.bulk.LF$group, 
         pch = c(21,21,21), cex = 1.2, lwd=1, col = c("black","black","black"))
  plotGroupEllipses(siber.LF, m= c(24,14,11), n = 100, p.interval = NULL,
                    lty = 1, lwd = 2)
  text(-19, 17, expression(paste(bold("(F)"), italic(" Lutjanus fulvus"))), adj = c(0,1))
  #legend("bottomleft",c("Very High", "Medium", "Very Low"), pch=c(21,21,21), 
  #       col = c("black","black","black"), pt.cex=1,
  #       pt.bg = alpha(c("#A60021","#ABF8FF","#2A0BD9"),0.75))
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


# calculate the SEA on the posterior distribution of covariance matrix for each group (i.e, Bayesian SEA or SEA-B) (Supplementary Table 3).
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


# > Plot SEA-B data (Extended Data Fi.g 1) -------------------------------------------------------
palette(c("#2A0BD9","#ABF8FF","#A60021"))

par(mfrow=c(1,1),mar=c(5,5,1,1))

siberDensityPlot(SEA.B, xticklabels = c(rep(c("V Lo","Med", "V Hi"),6)), 
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

text(2, 20, expression(italic("A. furca")))
text(5, 20, expression(italic("C. melampygus")))
text(8, 20, expression(italic("C. argus")))
text(11, 20, expression(italic("L. bohar")))
text(14, 20, expression(italic("C. urodeta")))
text(17, 20, expression(italic("L. fulvus")))



### Compare isotopic niche size (Extended Data Table 2) -------------------------------------------------------

# In order to test whether one group's ellipse is smaller or larger than another, we can simply calculate the probability 
# that its posterior distribution is smaller (or larger). This is achieved by comparing each pair of posterior draws for 
# both groups, and determining which is smaller in magnitude. We then find the proportion of draws that are smaller, and 
# this is a direct proxy for the probability that one group's posterior distribution (of ellipse size in this case) is 
# smaller than the other.

# calculate the proportion, and hence probability, of the SEA.B for group 1 being smaller than the SEA.B for group 2 (e.g., AF.1.2)
# calculate the proportion, and hence probability, of the SEA.B for group 1 being smaller than the SEA.B for group 3 (e.g., AF.1.3)
# calculate the proportion, and hence probability, of the SEA.B for group 2 being smaller than the SEA.B for group 2 (e.g., AF.2.3)

## A. furca
(AF.1.2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B))
(AF.1.3 <- sum( SEA.B[,1] < SEA.B[,3] ) / nrow(SEA.B))
(AF.2.3 <- sum( SEA.B[,2] < SEA.B[,3] ) / nrow(SEA.B))

## C. melampygus
(CM.1.2 <- sum( SEA.B[,4] < SEA.B[,5] ) / nrow(SEA.B))

## C. argus
(CA.1.2 <- sum( SEA.B[,7] < SEA.B[,8] ) / nrow(SEA.B))
(CA.1.3 <- sum( SEA.B[,7] < SEA.B[,9] ) / nrow(SEA.B))
(CA.2.3 <- sum( SEA.B[,8] < SEA.B[,9] ) / nrow(SEA.B))

## L. bohar
(LB.1.2 <- sum( SEA.B[,10] < SEA.B[,11] ) / nrow(SEA.B))
(LB.1.3 <- sum( SEA.B[,10] < SEA.B[,12] ) / nrow(SEA.B))
(LB.2.3 <- sum( SEA.B[,11] < SEA.B[,12] ) / nrow(SEA.B))

## C. urodeta
(CU.1.2 <- sum( SEA.B[,13] < SEA.B[,14] ) / nrow(SEA.B))
(CU.1.3 <- sum( SEA.B[,13] < SEA.B[,15] ) / nrow(SEA.B))
(CU.2.3 <- sum( SEA.B[,14] < SEA.B[,15] ) / nrow(SEA.B))

## L. fulvus
(LF.1.2 <- sum( SEA.B[,16] < SEA.B[,17] ) / nrow(SEA.B))
(LF.1.3 <- sum( SEA.B[,16] < SEA.B[,18] ) / nrow(SEA.B))
(LF.2.3 <- sum( SEA.B[,17] < SEA.B[,18] ) / nrow(SEA.B))


### Quantify isotopic niche overlap (Extended Data Table 3, Extended Data Fig. 2) -------------------------------------------------------

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



