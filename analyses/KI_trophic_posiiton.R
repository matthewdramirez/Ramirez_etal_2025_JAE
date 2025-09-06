
# R code for Ramirez et al. 2025 "Carbon pathways and trophic attributes are conserved in carnivorous reef fishes across a major human disturbance gradient" Journal of Animal Ecology

# Script to calculate and compare reef fish trophic positions using amino acid d15N values

# All figures were modified in Adobe Illustrator prior to publication.


##############################

## Load necessary packages
require(ggplot2)
require(tidyr)
require(RColorBrewer)
require(plyr)
require(cowplot)
require(propagate)
require(dplyr)
library(ggpubr)
require(arm)
require(car)

## Set your working directory
# Make sure that this contains the "ki_csia_n.csv" file
#setwd("C:/Users/...") # If on a PC
#setwd("/Users/...") # If on a Mac
setwd("~/Desktop/KI_fish/data")


# Load the AA d15N dataset
fish.N <- read.csv("ki_csia_n.csv")

fish.N$fish_code <- as.factor(fish.N$fish_code)
fish.N$dist_cat <- as.factor(fish.N$dist_cat)
fish.N$sci_name <- as.factor(fish.N$sci_name)

fish.N$dist_cat <- factor(fish.N$dist_cat, level = c("Very Low", "Medium", "Very High"))

# Rescale human disturbance to account for large gaps in values
fish.N$cp.z <- rescale(fish.N$continous.pressure.2km)

unique(fish.N$site)
unique(fish.N$pub.name)

fish.N$site <- factor(fish.N$site, level = c("10","11","16","15","19","14","8","35","34","27","30"))
# Note site 11 only in AA-C (Zp) data and site 14 NOT in AA-N data but in AA-C (He) and bulk SIA data
# So, AA-N data has VL = 4 sites, M = 3 sites, and VH = 2 sites

fish.N$pub.name <- factor(fish.N$pub.name, level = c("VL6","VL9","VL11","VL1","VL2","M4","M1","M2","M3","VH1","VH3"))
# Note site VL9 only in AA-C (Zp) data and site M4 NOT in AA-N data but in AA-C (He) and bulk SIA data
# So, AA-N data has VL = 4 sites, M = 3 sites, and VH = 2 sites


### Carnivorous fish source amino acid d15N patterns (Figure S7) -----------------------------------------------------------

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# By disturbance level
Phe <- ggplot(fish.N, aes(sci_name, Phe, fill = dist_cat)) + 
  #geom_boxplot(outlier.shape = 1, alpha = 0.75) +
  geom_point(position = position_dodge(0.75), pch=21, alpha = 0.75, size=2.5)+
  theme_test() +
  scale_fill_manual(values = c("#2A0BD9","#ABF8FF","#A60021")) +
  stat_summary(fun.data = function(x){return(data.frame(y = 10, label = length(x)))}, 
               geom = "text",
               aes(group=dist_cat),
               position = position_dodge(0.75),
               size = 3.5) +
  stat_summary(fun.data=data_summary, 
               color="black", pch = 16, size = 0.25,
               aes(group=dist_cat),
               position = position_dodge(0.75)) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  #theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  scale_y_continuous(limits=c(1,10),breaks = seq(1,9,2)) +
  labs(y=expression(Phenylalanine~{delta}^15*N~values~'(\u2030)')) +
  theme(plot.title.position = 'plot') +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(legend.position = 'top') + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Species") +
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  coord_cartesian(clip = 'off')


Lys <- ggplot(fish.N, aes(sci_name, Lys, fill = dist_cat)) + 
  #geom_boxplot(outlier.shape = 1,  alpha = 0.75) +
  geom_point(position = position_dodge(0.75), pch=21, alpha = 0.75, size=2.5)+
  theme_test() +
  scale_fill_manual(values = c("#2A0BD9","#ABF8FF","#A60021")) +
  stat_summary(fun.data = function(x){return(data.frame(y = 10, label = length(x)))}, 
               geom = "text",
               aes(group=dist_cat),
               position = position_dodge(0.75),
               size = 3.5) +
  stat_summary(fun.data=data_summary, 
               color="black", pch = 16, size = 0.25,
               aes(group=dist_cat),
               position = position_dodge(0.75)) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  #theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  scale_y_continuous(limits=c(1,10),breaks = seq(1,9,2)) +
  labs(y=expression(Lysine~{delta}^15*N~values~'(\u2030)')) +
  theme(plot.title.position = 'plot') +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(legend.position = 'top') + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Species") +
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  coord_cartesian(clip = 'off')

plot_grid(Phe, Lys, align = "h", nrow = 1, axis = "tb")



# by site

Phe2 <- ggplot(fish.N, aes(sci_name, Phe, fill = pub.name, colour=pub.name, shape = pub.name)) + 
  #geom_boxplot(outlier.shape = 1, alpha = 0.75, colour= "black") +
  geom_point(position = position_dodge(0.75), size=2.5)+
  theme_test() +
  scale_fill_manual(values = c(rep("#2A0BD9",4), rep("#ABF8FF",3), rep("#A60021",2)))+
  scale_colour_manual(values = c(rep("#2A0BD9",4), rep("#ABF8FF",3), rep("#A60021",2)))+
  scale_shape_manual(values = c(0,1,2,5,21:23,7,9))+
  stat_summary(fun.data = function(x){return(data.frame(y = 10, label = length(x)))}, 
               geom = "text",
               aes(group=pub.name),
               position = position_dodge(0.75),
               size = 3.5) +
  stat_summary(fun.data=data_summary, 
               color="black", pch = 16, size = 0.25,
               aes(group=pub.name),
               position = position_dodge(0.75)) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  #theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  scale_y_continuous(limits=c(1,10),breaks = seq(1,9,2)) +
  labs(y=expression(Phenylalanine~{delta}^15*N~values~'(\u2030)')) +
  theme(plot.title.position = 'plot') +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(legend.position = 'top') + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Species") +
  guides(fill = guide_legend(title = "Site", title.position = "top", title.hjust = 0.5)) +
  coord_cartesian(clip = 'off')

Lys2 <- ggplot(fish.N, aes(sci_name, Lys, fill = pub.name, colour=pub.name, shape = pub.name)) + 
  #geom_boxplot(outlier.shape = 1, alpha = 0.75, colour = "black") +
  geom_point(position = position_dodge(0.75), size=2.5)+
  theme_test() +
  scale_fill_manual(values = c(rep("#2A0BD9",4), rep("#ABF8FF",3), rep("#A60021",2)))+
  scale_colour_manual(values = c(rep("#2A0BD9",4), rep("#ABF8FF",3), rep("#A60021",2)))+
  scale_shape_manual(values = c(0,1,2,5,21:23,7,9))+
  stat_summary(fun.data = function(x){return(data.frame(y = 10, label = length(x)))}, 
               geom = "text",
               aes(group=pub.name),
               position = position_dodge(0.75),
               size = 3.5) +
  stat_summary(fun.data=data_summary, 
               color="black", pch = 16, size = 0.25,
               aes(group=pub.name),
               position = position_dodge(0.75)) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  #theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  scale_y_continuous(limits=c(1,10),breaks = seq(1,9,2)) +
  labs(y=expression(Lysine~{delta}^15*N~values~'(\u2030)')) +
  theme(plot.title.position = 'plot') +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(legend.position = 'top') + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Species") +
  guides(fill = guide_legend(title = "Site", title.position = "top", title.hjust = 0.5)) +
  coord_cartesian(clip = 'off')

plot_grid(Phe2, Lys2, align = "h", nrow = 1, axis = "tb")


### Estimate carnivorous fish trophic positions and error -----------------------------------------------------------

# Install devtools, if you haven't already.
#install.packages("devtools")
#library(devtools)
#install_github("anspiess/propagate")
#source("https://install-github.me/anspiess/propagate")


# Define TDF and Beta values:

# TDFs from Bradley et al. 2015 (marine teleosts), Betas from Ramirez et al. 2021
# Parms ordered as: (Beta, Beta_SD, TDF, TDF_SD, TDF2, TDF2_SD) 

#parms_glu_phe <- c(3.3, 1.8, 5.7, 0.3) # TP estimated via Glu-Phe only
parms_tr_sr <- c(3.0, 2.4, 5.5, 0.5) # TP estimate from mean of Trophic (Ala,Val,Leu,Pro,Glu) and Source (Phe,Lys) AAs

# Expression for which to propagate error:

# Single TDF Equations
#EXPR1 <- expression((Glu - Phe - Beta)/TDF + 1)
EXPR2 <- expression((Tr - Sr - Beta)/TDF + 1)


propResults <- lapply(1:nrow(fish.N), function(i){
  # For each row in the dataframe, i.e. for each individual, make a data.frame with the individual's values 
  # in the first row and the error in the second.
  #   
  # (This is a slightly non-intuitive format, but it is what `propagate` expects).
  # See `?propagate` for more info
  DAT <- data.frame(
    # Comment on/off for Glu-Phe only
    #Glu  = c(fish.N[i,"Glu"], fish.N[i, "Glu_SD"]),
    #Phe  = c(fish.N[i,"Phe"], fish.N[i, "Phe_SD"]),
    #Beta = c(parms_glu_phe[1],parms_glu_phe[2]),    
    #TDF  = c(parms_glu_phe[3],parms_glu_phe[4]),

    Tr  = c(fish.N[i,"Tr_avg"], fish.N[i, "Tr_SD"]),
    Sr  = c(fish.N[i,"Sr_avg"], fish.N[i, "Sr_SD"]),
    Beta = c(parms_tr_sr[1],parms_tr_sr[2]),    
    TDF  = c(parms_tr_sr[3],parms_tr_sr[4])
  )
  
  # Conduct uncertainty propagation using default settings (see `?propagate` for more options)
  # Remember to REPLACE first entry in formula with correct expression/equaiton for your analysis
  # nsim is 10000 here for demonstration, but was set to 1000000 for analyses
  res <- propagate(EXPR2, as.matrix(DAT), second.order=FALSE, do.sim=TRUE, cov=TRUE, df=NULL, 
                   nsim=10000, alpha=0.05)
  
  # Output the results from the uncertainty propagation for that individual, bound to the original data.frame columns:
  propResults_stitched <- cbind( fish.N[i, ], rbind(res$prop) )
  
  # Change the names of min and max CI's (`2.5%` and `97.5%` are bad column names!)
  names(propResults_stitched)[
    {ncol(propResults_stitched)-1}:ncol(propResults_stitched)
  ] <- c("CImin", "CImax")
  return(propResults_stitched)
})


#setwd("~/Desktop/KI_fish/data")

# Finally, put back as a data.frame.
TPcsia <- dplyr::bind_rows(propResults)
#write.csv(TPcsia, "TPcsia.csv")


### Trophic position comparisons -----------------------------------------------------------

#TPcsia <- read.csv("TPcsia.csv")

summary(TPcsia)

TPcsia$sci_name <- as.factor(TPcsia$sci_name)
TPcsia$dist_cat <- as.factor(TPcsia$dist_cat)

TPcsia$dist_cat <- factor(TPcsia$dist_cat, level = c("Very Low", "Medium", "Very High"))
TPcsia$sci_name <- factor(TPcsia$sci_name, level = c("Aphareus furca", "Caranx melampygus","Cephalopholis argus",  "Lutjanus bohar", "Cephalopholis urodeta", "Lutjanus fulvus"))

levels(TPcsia$sci_name)
new_labels <- c("Aphareus furca" = "A. furca", 
                "Caranx melampygus" = "C. melampygus", 
                "Cephalopholis argus" = "C. argus", 
                "Lutjanus bohar" = "L. bohar", 
                "Cephalopholis urodeta" = "C. urodeta", 
                "Lutjanus fulvus" = "L. fulvus")

#summary(TPcsia$Mean.1)


# > Plot standard length vs. TPcsia -----------------------------------------------------------
ggplot(TPcsia, aes(std_len_mm, Mean.1)) +
  geom_smooth(method="gam", colour="black", alpha=0.25) +
  geom_point(aes(fill=dist_cat), pch = 21, col = "black", size = 4, alpha = 0.75) + 
  scale_fill_manual(values = c("#2A0BD9", "#ABF8FF", "#A60021"))+
  theme_test() +
  #theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(legend.position = 'top') + 
  xlab("Standard Length (mm)") +
  ylab("Trophic Position") +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  theme(plot.title.position = 'plot') +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black")) +
  theme(legend.position = 'right') + 
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  guides(colour = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5))


# > Plot standard lenggth vs. TPcsia by species (Figure S1) -----------------------------------------------------------
ggplot(TPcsia, aes(std_len_mm/10, Mean.1)) +
  geom_smooth(method="glm", colour="black", alpha=0.25)+
  geom_point(aes(fill=dist_cat), pch = 21, col = "black", size = 4, alpha = 0.75) + 
  scale_fill_manual(values = c("#2A0BD9", "#ABF8FF", "#A60021"))+
  #stat_regline_equation(label.y = 4.6, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 4.5, aes(label = ..rr.label..), size = 3.25) + 
  theme_test() +
  #theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(legend.position = 'top') + 
  xlab("Standard Length (cm)") +
  ylab(expression(paste(TP[CSIA-AVG]))) +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black")) +
  theme(strip.text = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) + 
  facet_wrap(~sci_name, scales = "free_x", nrow=1, labeller=labeller(sci_name = new_labels)) 


# > Plot UNCORRECTED trophic position by species and dist_cat -----------------------------------------------------------
ggplot(TPcsia, aes(sci_name, Mean.1, fill = dist_cat, shape =dist_cat)) + 
  geom_errorbar(stat="summary", 
                fun.min=function(x) {mean(x)-1.96* sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+1.96*sd(x)/sqrt(length(x))},
                position = position_dodge(0.75), width = 0, lwd=1) + 
  stat_summary(fun.data = function(x){return(data.frame(y = 4.75, label = length(x)))}, 
               geom = "text", aes(group=dist_cat),
               position = position_dodge(0.75), 
               col = "black") +
  geom_point(stat="summary", fun="mean", position = position_dodge(0.75), pch = 21, col = "black", size=5, alpha = 0.75) + 
  scale_fill_manual(values = c("#2A0BD9", "#ABF8FF", "#A60021"))+
  scale_y_continuous(limits=c(3,4.8),breaks = seq(3,4.8,0.4)) +
  theme_classic()+theme(text = element_text(size=14))+
  #theme(axis.line = element_line(colour = 'black', size = .5))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  theme_test() +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  scale_x_discrete(labels = c("A. furca", "C. melampygus", "C. argus", "L. bohar", "C. urodeta", "L. fulvus")) +
  xlab("Species") + 
  ylab(expression(paste(TP[CSIA-AVG]))) +
  theme(legend.position = 'top') + 
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black")) +
  theme(strip.text = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5), lty="dotted")


# > Correct TPcsia to mean standard length for each species -----------------------------------------------------------

## GLMs to characeterize linear relationships between TP and standard length. 
## Coefficients used to correct TP data to mean standard length for each species
data.AF <- subset(TPcsia, species_code == "AF")
data.CM <- subset(TPcsia, species_code == "CM")
data.CA <- subset(TPcsia, species_code == "CA")
data.LB <- subset(TPcsia, species_code == "LB")
data.CU <- subset(TPcsia, species_code == "CU")
data.LF <- subset(TPcsia, species_code == "LF")

data.AF$std_len_cm <- data.AF$std_len_mm/10
data.CM$std_len_cm <- data.CM$std_len_mm/10
data.CA$std_len_cm <- data.CA$std_len_mm/10
data.LB$std_len_cm <- data.LB$std_len_mm/10
data.CU$std_len_cm <- data.CU$std_len_mm/10
data.LF$std_len_cm <- data.LF$std_len_mm/10

mod.AF <- lm(Mean.1 ~ std_len_cm, data = data.AF)
summary(mod.AF)

mod.CM <- lm(Mean.1 ~ std_len_cm, data = data.CM)
summary(mod.CM)

mod.CA <- lm(Mean.1 ~ std_len_cm, data = data.CA)
summary(mod.CA)

mod.LB <- lm(Mean.1 ~ std_len_cm, data = data.LB)
summary(mod.LB)

mod.CU <- lm(Mean.1 ~ std_len_cm, data = data.CU)
summary(mod.CU)

mod.LF <- lm(Mean.1 ~ std_len_cm, data = data.LF)
summary(mod.LF)


# Use coefficients from linear models to correct TPcsia to species mean standard length
data.AF$TP_cor <- data.AF$Mean.1+(mean(data.AF$std_len_cm)-data.AF$std_len_cm)*mod.AF$coefficients[2] 
data.CM$TP_cor <- data.CM$Mean.1+(mean(data.CM$std_len_cm)-data.CM$std_len_cm)*mod.CM$coefficients[2] 
data.CA$TP_cor <- data.CA$Mean.1+(mean(data.CA$std_len_cm)-data.CA$std_len_cm)*mod.CA$coefficients[2] 
data.LB$TP_cor <- data.LB$Mean.1+(mean(data.LB$std_len_cm)-data.LB$std_len_cm)*mod.LB$coefficients[2] 
data.CU$TP_cor <- data.CU$Mean.1+(mean(data.CU$std_len_cm)-data.CU$std_len_cm)*mod.CU$coefficients[2] 
data.LF$TP_cor <- data.LF$Mean.1+(mean(data.LF$std_len_cm)-data.LF$std_len_cm)*mod.LF$coefficients[2] 


TPcsia.cor <- rbind(data.AF,
                    data.CM,
                    data.CA,
                    data.LB,
                    data.CU,
                    data.LF)

summary(TPcsia.cor$TP_cor)

summary(data.AF$TP_cor)
summary(data.CM$TP_cor)
summary(data.CA$TP_cor)
summary(data.LB$TP_cor)
summary(data.CU$TP_cor)
summary(data.LF$TP_cor)

# Calculate range of TP_cor per species
diff1 <- round(max(data.AF$TP_cor)-min(data.AF$TP_cor),2)
diff2 <- round(max(data.CM$TP_cor)-min(data.CM$TP_cor),2)
diff3 <- round(max(data.CA$TP_cor)-min(data.CA$TP_cor),2)
diff4 <- round(max(data.LB$TP_cor)-min(data.LB$TP_cor),2)
diff5 <- round(max(data.CU$TP_cor)-min(data.CU$TP_cor),2)
diff6 <- round(max(data.LF$TP_cor)-min(data.LF$TP_cor),2)

# Calculate mean range of TP_cor per species
mean(c(diff1,diff2,diff3,diff4,diff5,diff6))

# Calculate diff in mean TP_cor between "Very High" and "Very Low" per species
diff1.2 <- round(mean(data.AF$TP_cor[data.AF$dist_cat == "Very High"])-mean(data.AF$TP_cor[data.AF$dist_cat == "Very Low"]),2)
diff2.2 <- round(mean(data.CM$TP_cor[data.CM$dist_cat == "Very High"])-mean(data.CM$TP_cor[data.CM$dist_cat == "Very Low"]),2)
diff3.2 <- round(mean(data.CA$TP_cor[data.CA$dist_cat == "Very High"])-mean(data.CA$TP_cor[data.CA$dist_cat == "Very Low"]),2)
diff4.2 <- round(mean(data.LB$TP_cor[data.LB$dist_cat == "Very High"])-mean(data.LB$TP_cor[data.LB$dist_cat == "Very Low"]),2)
diff5.2 <- round(mean(data.CU$TP_cor[data.CU$dist_cat == "Very High"])-mean(data.CU$TP_cor[data.CU$dist_cat == "Very Low"]),2)
diff6.2 <- round(mean(data.LF$TP_cor[data.LF$dist_cat == "Very High"])-mean(data.LF$TP_cor[data.LF$dist_cat == "Very Low"]),2)

# Calculate mean diff in mean TP_cor between "Very High" and "Very Low" per species
mean(c(diff1.2,abs(diff2.2),abs(diff3.2),diff4.2,diff5.2,diff6.2))

diff1.3 <- round(mean(data.AF$TP_cor[data.AF$dist_cat == "Very High"])-mean(data.AF$TP_cor[data.AF$dist_cat == "Medium"]),2)
diff2.3 <- round(mean(data.CM$TP_cor[data.CM$dist_cat == "Very High"])-mean(data.CM$TP_cor[data.CM$dist_cat == "Medium"]),2)
diff3.3 <- round(mean(data.CA$TP_cor[data.CA$dist_cat == "Very High"])-mean(data.CA$TP_cor[data.CA$dist_cat == "Medium"]),2)
diff4.3 <- round(mean(data.LB$TP_cor[data.LB$dist_cat == "Very High"])-mean(data.LB$TP_cor[data.LB$dist_cat == "Medium"]),2)
diff5.3 <- round(mean(data.CU$TP_cor[data.CU$dist_cat == "Very High"])-mean(data.CU$TP_cor[data.CU$dist_cat == "Medium"]),2)
diff6.3 <- round(mean(data.LF$TP_cor[data.LF$dist_cat == "Very High"])-mean(data.LF$TP_cor[data.LF$dist_cat == "Medium"]),2)

diff1.4 <- round(mean(data.AF$TP_cor[data.AF$dist_cat == "Very Low"])-mean(data.AF$TP_cor[data.AF$dist_cat == "Medium"]),2)
diff2.4 <- round(mean(data.CM$TP_cor[data.CM$dist_cat == "Very Low"])-mean(data.CM$TP_cor[data.CM$dist_cat == "Medium"]),2)
diff3.4 <- round(mean(data.CA$TP_cor[data.CA$dist_cat == "Very Low"])-mean(data.CA$TP_cor[data.CA$dist_cat == "Medium"]),2)
diff4.4 <- round(mean(data.LB$TP_cor[data.LB$dist_cat == "Very Low"])-mean(data.LB$TP_cor[data.LB$dist_cat == "Medium"]),2)
diff5.4 <- round(mean(data.CU$TP_cor[data.CU$dist_cat == "Very Low"])-mean(data.CU$TP_cor[data.CU$dist_cat == "Medium"]),2)
diff6.4 <- round(mean(data.LF$TP_cor[data.LF$dist_cat == "Very Low"])-mean(data.LF$TP_cor[data.LF$dist_cat == "Medium"]),2)

mean(c(abs(diff1.2),
     abs(diff2.2),
     abs(diff3.2),
     abs(diff4.2),
     abs(diff5.2),
     abs(diff6.2),
     abs(diff1.3),
     abs(diff2.3),
     abs(diff3.3),
     abs(diff4.3),
     abs(diff5.3),
     abs(diff6.3),
     abs(diff1.4),
     abs(diff2.4),
     abs(diff3.4),
     abs(diff4.4),
     abs(diff5.4),
     abs(diff6.4)))

min(c(abs(diff1.2),
       abs(diff2.2),
       abs(diff3.2),
       abs(diff4.2),
       abs(diff5.2),
       abs(diff6.2),
       abs(diff1.3),
       abs(diff2.3),
       abs(diff3.3),
       abs(diff4.3),
       abs(diff5.3),
       abs(diff6.3),
       abs(diff1.4),
       abs(diff2.4),
       abs(diff3.4),
       abs(diff4.4),
       abs(diff5.4),
       abs(diff6.4)))

max(c(abs(diff1.2),
       abs(diff2.2),
       abs(diff3.2),
       abs(diff4.2),
       abs(diff5.2),
       abs(diff6.2),
       abs(diff1.3),
       abs(diff2.3),
       abs(diff3.3),
       abs(diff4.3),
       abs(diff5.3),
       abs(diff6.3),
       abs(diff1.4),
       abs(diff2.4),
       abs(diff3.4),
       abs(diff4.4),
       abs(diff5.4),
       abs(diff6.4)))

# > Plot CORRECTED  trophic position by species and dist_cat (Figure 5) -----------------------------------------------------------
ggplot(TPcsia.cor, aes(sci_name, TP_cor, fill = dist_cat, shape =dist_cat)) + 
  geom_errorbar(stat="summary", 
                fun.min=function(x) {mean(x)-1.96* sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+1.96*sd(x)/sqrt(length(x))},
                position = position_dodge(0.75), width = 0, lwd=1) + 
  stat_summary(fun.data = function(x){return(data.frame(y = 4.45, label = length(x)))}, 
               geom = "text", aes(group=dist_cat),
               position = position_dodge(0.75), 
               col = "black") +
  geom_point(stat="summary", fun="mean", position = position_dodge(0.75), pch = 21, col = "black", size=5, alpha = 0.75) + 
  scale_fill_manual(values = c("#2A0BD9", "#ABF8FF", "#A60021"))+
  scale_y_continuous(limits=c(3.2,4.5),breaks = seq(3.2,4.4,0.4)) +
  theme_classic()+theme(text = element_text(size=14))+
  #theme(axis.line = element_line(colour = 'black', size = .5))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  theme_test() +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  scale_x_discrete(labels = c("A. furca", "C. melampygus", "C. argus", "L. bohar", "C. urodeta", "L. fulvus")) +
  xlab("Species") + 
  ylab(expression(paste(TP[CSIA]))) +
  theme(legend.position = 'top') + 
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black")) +
  guides(fill = guide_legend(title = "Disturbance Level", title.position = "top", title.hjust = 0.5)) +
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5), lty="dotted")


# > Plot CORRECTED  trophic position by species and SITE (Figure S6) -----------------------------------------------------------
ggplot(TPcsia.cor, aes(sci_name, TP_cor, fill = pub.name, colour=pub.name, shape = pub.name)) + 
  geom_errorbar(stat="summary", 
                fun.min=function(x) {mean(x)-1.96* sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+1.96*sd(x)/sqrt(length(x))},
                position = position_dodge(0.75), colour = "black", width = 0, lwd=1) + 
  stat_summary(fun.data = function(x){return(data.frame(y = 4.45, label = length(x)))}, 
               geom = "text", aes(group=pub.name),
               position = position_dodge(0.75), 
               col = "black") +
  geom_point(stat="summary", fun="mean", position = position_dodge(0.75),  size=5, alpha = 0.75) + 
  scale_fill_manual(values = c(rep("#2A0BD9",4), rep("#ABF8FF",3), rep("#A60021",2)))+
  scale_colour_manual(values = c(rep("#2A0BD9",4), rep("#ABF8FF",3), rep("#A60021",2)))+
  scale_shape_manual(values = c(0,1,2,5,21:23,7,9))+
  scale_y_continuous(limits=c(3.2,4.5),breaks = seq(3.2,4.4,0.4)) +
  theme_classic()+theme(text = element_text(size=14))+
  #theme(axis.line = element_line(colour = 'black', size = .5))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  theme_test() +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=12)) +
  scale_x_discrete(labels = c("A. furca", "C. melampygus", "C. argus", "L. bohar", "C. urodeta", "L. fulvus")) +
  xlab("Species") + 
  ylab(expression(paste(TP[CSIA]))) +
  theme(legend.position = 'top') + 
  theme(text=element_text(size=12)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(panel.border = element_rect(color = "black")) +
  guides(fill = guide_legend(title = "Site", title.position = "top", title.hjust = 0.5)) +
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5), lty="dotted")


# > Kruskal-Wallis test comparing TPcsia by species and dist_cat -----------------------------------------------------------

# Uncorrected TPcsia
# If the p-value is less than the significance level 0.05, we can conclude that there are significant differences between the treatment groups.
kruskal.test(Mean.1 ~ dist_cat, data = data.AF) #p-value = 0.005559
kruskal.test(Mean.1 ~ dist_cat, data = data.CM) #p-value = 0.09408
kruskal.test(Mean.1 ~ dist_cat, data = data.CA) #p-value = 0.1959
kruskal.test(Mean.1 ~ dist_cat, data = data.LB) #p-value = 0.8415
kruskal.test(Mean.1 ~ dist_cat, data = data.CU) #p-value = 0.5655
kruskal.test(Mean.1 ~ dist_cat, data = data.LF) #p-value = 0.5434

# Body size-corrected TPcsia
# If the p-value is less than the significance level 0.05, we can conclude that there are significant differences between the treatment groups.
kruskal.test(TP_cor ~ dist_cat, data = data.AF) #p-value = 0.009577
kruskal.test(TP_cor ~ dist_cat, data = data.CM) #p-value = 0.5159
kruskal.test(TP_cor ~ dist_cat, data = data.CA) #p-value = 0.6126
kruskal.test(TP_cor ~ dist_cat, data = data.LB) #p-value = 0.6938
kruskal.test(TP_cor ~ dist_cat, data = data.CU) #p-value = 0.1212
kruskal.test(TP_cor ~ dist_cat, data = data.LF) #p-value = 0.2491


# > General linear models comparing CORRECTED TP by species and dist_cat (Table 2) -----------------------------------------------------------

# Compare models with continuous vs. discrete dist_cat data
# Results of models with continuous dist_cat data reported in Table 2

# Aphareus furca
mod1.2 <- lm(TP_cor ~ cp.z, data = data.AF)
mod1.3 <- lm(TP_cor ~ dist_cat, data = data.AF)

summary(mod1.2)
summary(mod1.3)

AIC(mod1.2)
AIC(mod1.3)

anova(mod1.2,mod1.3,test = "F") # Results not significantly different

myAIC <- round(AIC(mod1.2,mod1.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]

# Compute the analysis of variance
AF.aov <- aov(TP_cor ~ dist_cat, data = data.AF)
# Summary of the analysis
summary(AF.aov)
# multiple pairwise-comparison between the means of groups
TukeyHSD(AF.aov)

# Test for Homogeneity of Variance
leveneTest(TP_cor ~ dist_cat, data = data.AF) #p > 0.05...no evidence to suggest that the variance across groups is statistically significantly differen

plot(AF.aov, 1)

# Test for Normality
# Extract the residuals
AF.aov_residuals <- residuals(object = AF.aov )
# Run Shapiro-Wilk test
shapiro.test(x = AF.aov_residuals ) # p > 0.05...no indication that normality is violated.

plot(AF.aov, 2)



# Caranx melampygus
mod2.2 <- lm(TP_cor ~ cp.z, data = data.CM)
mod2.3 <- lm(TP_cor ~ dist_cat, data = data.CM)

summary(mod2.2)
summary(mod2.3)

AIC(mod2.2)
AIC(mod2.3)

anova(mod2.2,mod2.3,test = "F") # Results not significantly different

myAIC <- round(AIC(mod2.2,mod2.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]

# Compute the analysis of variance
CM.aov <- aov(TP_cor ~ dist_cat, data = data.CM)
# Summary of the analysis
summary(CM.aov)
# multiple pairwise-comparison between the means of groups
TukeyHSD(CM.aov)

# Test for Homogeneity of Variance
leveneTest(TP_cor ~ dist_cat, data = data.CM) #p > 0.05...no evidence to suggest that the variance across groups is statistically significantly differen

plot(CM.aov, 1)

# Test for Normality
# Extract the residuals
CM.aov_residuals <- residuals(object = CM.aov )
# Run Shapiro-Wilk test
shapiro.test(x = CM.aov_residuals ) # p > 0.05...no indication that normality is violated.

plot(CM.aov, 2)




# Cephalopholis argus
mod3.2 <- lm(TP_cor ~ cp.z, data = data.CA)
mod3.3 <- lm(TP_cor ~ dist_cat, data = data.CA)

summary(mod3.2)
summary(mod3.3)

AIC(mod3.2)
AIC(mod3.3)

anova(mod3.2,mod3.3,test = "F") # Results not significantly different

myAIC <- round(AIC(mod3.2,mod3.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]

# Compute the analysis of variance
CA.aov <- aov(TP_cor ~ dist_cat, data = data.CA)
# Summary of the analysis
summary(CA.aov)
# multiple pairwise-comparison between the means of groups
TukeyHSD(CA.aov)

# Test for Homogeneity of Variance
leveneTest(TP_cor ~ dist_cat, data = data.CA) #p > 0.05...no evidence to suggest that the variance across groups is statistically significantly differen

plot(CA.aov, 1)

# Test for Normality
# Extract the residuals
CA.aov_residuals <- residuals(object = CA.aov )
# Run Shapiro-Wilk test
shapiro.test(x = CA.aov_residuals ) # p > 0.05...no indication that normality is violated.

plot(CA.aov, 2)



# Lutjanus bohar
mod4.2 <- lm(TP_cor ~ cp.z, data = data.LB)
mod4.3 <- lm(TP_cor ~ dist_cat, data = data.LB)

summary(mod4.2)
summary(mod4.3)

AIC(mod4.2)
AIC(mod4.3)

anova(mod4.2,mod4.3,test = "F") # Results not significantly different

myAIC <- round(AIC(mod4.2,mod4.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]

# Compute the analysis of variance
#data.LB <- subset(data.LB, fish_code != "LB8") # must drop fish LB8 of normality assumption to be upheld
LB.aov <- aov(TP_cor ~ dist_cat, data = data.LB)
# Summary of the analysis
summary(LB.aov)
# multiple pairwise-comparison between the means of groups
TukeyHSD(LB.aov)

# Test for Homogeneity of Variance
leveneTest(TP_cor ~ dist_cat, data = data.LB) #p > 0.05...no evidence to suggest that the variance across groups is statistically significantly differen

plot(LB.aov, 1)

# Test for Normality
# Extract the residuals
LB.aov_residuals <- residuals(object = LB.aov )
# Run Shapiro-Wilk test
shapiro.test(x = LB.aov_residuals ) # p > 0.05...no indication that normality is violated.

plot(LB.aov, 2)



# Cephalopholis urodeta
mod5.2 <- lm(TP_cor ~ cp.z, data = data.CU)
mod5.3 <- lm(TP_cor ~ dist_cat, data = data.CU)

summary(mod5.2)
summary(mod5.3)

AIC(mod5.2)
AIC(mod5.3)

anova(mod5.2,mod5.3,test = "F") # Model with categorical variable has lower AIC than model with continuous variable

myAIC <- round(AIC(mod5.2,mod5.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]

# Compute the analysis of variance
CU.aov <- aov(TP_cor ~ dist_cat, data = data.CU)
# Summary of the analysis
summary(CU.aov)
# multiple pairwise-comparison between the means of groups
TukeyHSD(CU.aov)

# Test for Homogeneity of Variance
leveneTest(TP_cor ~ dist_cat, data = data.CU) #p > 0.05...no evidence to suggest that the variance across groups is statistically significantly differen

plot(CU.aov, 1)

# Test for Normality
# Extract the residuals
CU.aov_residuals <- residuals(object = CU.aov )
# Run Shapiro-Wilk test
shapiro.test(x = CU.aov_residuals ) # p > 0.05...no indication that normality is violated.

plot(CU.aov, 2)



# Lutjanus fulvus
mod6.2 <- lm(TP_cor ~ cp.z, data = data.LF)
mod6.3 <- lm(TP_cor ~ dist_cat, data = data.LF)

summary(mod6.2)
summary(mod6.3)

AIC(mod6.2)
AIC(mod6.3)

anova(mod6.2,mod6.3,test = "F") # Results not significantly different

myAIC <- round(AIC(mod6.2,mod6.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]

# Compute the analysis of variance
LF.aov <- aov(TP_cor ~ dist_cat, data = data.LF)
# Summary of the analysis
summary(LF.aov)
# multiple pairwise-comparison between the means of groups
TukeyHSD(LF.aov)

# Test for Homogeneity of Variance
leveneTest(TP_cor ~ dist_cat, data = data.LF) #p > 0.05...no evidence to suggest that the variance across groups is statistically significantly differen

plot(LF.aov, 1)

# Test for Normality
# Extract the residuals
LF.aov_residuals <- residuals(object = LF.aov )
# Run Shapiro-Wilk test
shapiro.test(x = LF.aov_residuals ) # p > 0.05...no indication that normality is violated.

plot(LF.aov, 2)


# > General linear models comparing UNCORRECTED TP by species and dist_cat -----------------------------------------------------------

# >> Continuous data -----------------------------------------------------------

# Check for potential interactive effects
# Interaction term never significant
mod1.1 <- lm(Mean.1 ~ cp.z + std_len_cm + cp.z:std_len_cm, data = data.AF)
mod2.1 <- lm(Mean.1 ~ cp.z + std_len_cm + cp.z:std_len_cm, data = data.CM)
mod3.1 <- lm(Mean.1 ~ cp.z + std_len_cm + cp.z:std_len_cm, data = data.CA)
mod4.1 <- lm(Mean.1 ~ cp.z + std_len_cm + cp.z:std_len_cm, data = data.LB)
mod5.1 <- lm(Mean.1 ~ cp.z + std_len_cm + cp.z:std_len_cm, data = data.CU)
mod6.1 <- lm(Mean.1 ~ cp.z + std_len_cm + cp.z:std_len_cm, data = data.LF)

summary(mod1.1)
summary(mod2.1)
summary(mod3.1)
summary(mod4.1)
summary(mod5.1)
summary(mod6.1)

# Compare models with and without length term

# Aphareus furca
mod1.2 <- lm(Mean.1 ~ cp.z + std_len_cm, data = data.AF)
mod1.3 <- lm(Mean.1 ~ cp.z, data = data.AF)

summary(mod1.2)
summary(mod1.3)

AIC(mod1.2)
AIC(mod1.3)

anova(mod1.2,mod1.3,test = "F") # Adding body size DOES significantly improved fit

myAIC <- round(AIC(mod1.1, mod1.2,mod1.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Caranx melampygus
mod2.2 <- lm(Mean.1 ~ cp.z + std_len_cm, data = data.CM)
mod2.3 <- lm(Mean.1 ~ cp.z, data = data.CM)

summary(mod2.2)
summary(mod2.3)

AIC(mod2.2)
AIC(mod2.3)

anova(mod2.2,mod2.3,test = "F") # Adding body size DOES significantly improve fit (but close)

myAIC <- round(AIC(mod2.1, mod2.2,mod2.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Cephalopholis argus
mod3.2 <- lm(Mean.1 ~ cp.z + std_len_cm, data = data.CA)
mod3.3 <- lm(Mean.1 ~ cp.z, data = data.CA)

summary(mod3.2)
summary(mod3.3)

AIC(mod3.2)
AIC(mod3.3)

anova(mod3.2,mod3.3,test = "F") # Adding body size did NOT significantly improve fit; stick with simpler model

myAIC <- round(AIC(mod3.1, mod3.2,mod3.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Lutjanus bohar
mod4.2 <- lm(Mean.1 ~ cp.z + std_len_cm, data = data.LB)
mod4.3 <- lm(Mean.1 ~ cp.z, data = data.LB)

summary(mod4.2)
summary(mod4.3)

AIC(mod4.2)
AIC(mod4.3)

anova(mod4.2,mod4.3,test = "F") # Adding body size did NOT significantly improve fit; stick with simpler model

myAIC <- round(AIC(mod4.1, mod4.2,mod4.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Cephalopholis urodeta
mod5.2 <- lm(Mean.1 ~ cp.z + std_len_cm, data = data.CU)
mod5.3 <- lm(Mean.1 ~ cp.z, data = data.CU)

summary(mod5.2)
summary(mod5.3)

AIC(mod5.2)
AIC(mod5.3)

anova(mod5.2,mod5.3,test = "F") # Adding body size DID significantly improve fit

myAIC <- round(AIC(mod5.1, mod5.2,mod5.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Lutjanus fulvus
mod6.2 <- lm(Mean.1 ~ cp.z + std_len_cm, data = data.LF)
mod6.3 <- lm(Mean.1 ~ cp.z, data = data.LF)

summary(mod6.2)
summary(mod6.3)

AIC(mod6.2)
AIC(mod6.3)

anova(mod6.2,mod6.3,test = "F") # Adding body size DID significantly improve fit

myAIC <- round(AIC(mod6.1, mod6.2,mod6.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]


# >> Categorical data -----------------------------------------------------------

# Check for potential interactive effects
# Interaction term never significant
mod1.1 <- lm(Mean.1 ~ dist_cat + std_len_cm + dist_cat:std_len_cm, data = data.AF)
mod2.1 <- lm(Mean.1 ~ dist_cat + std_len_cm + dist_cat:std_len_cm, data = data.CM)
mod3.1 <- lm(Mean.1 ~ dist_cat + std_len_cm + dist_cat:std_len_cm, data = data.CA)
mod4.1 <- lm(Mean.1 ~ dist_cat + std_len_cm + dist_cat:std_len_cm, data = data.LB)
mod5.1 <- lm(Mean.1 ~ dist_cat + std_len_cm + dist_cat:std_len_cm, data = data.CU)
mod6.1 <- lm(Mean.1 ~ dist_cat + std_len_cm + dist_cat:std_len_cm, data = data.LF)

summary(mod1.1)
summary(mod2.1)
summary(mod3.1)
summary(mod4.1)
summary(mod5.1)
summary(mod6.1)

# Compare models with and without length term

# Aphareus furca
mod1.2 <- lm(Mean.1 ~ dist_cat + std_len_cm, data = data.AF)
mod1.3 <- lm(Mean.1 ~ dist_cat, data = data.AF)

summary(mod1.2)
summary(mod1.3)

AIC(mod1.2)
AIC(mod1.3)

anova(mod1.2,mod1.3,test = "F") # Adding body size did NOT significantly improve fit; stick with simpler model

myAIC <- round(AIC(mod1.1, mod1.2,mod1.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Caranx melampygus
mod2.2 <- lm(Mean.1 ~ dist_cat + std_len_cm, data = data.CM)
mod2.3 <- lm(Mean.1 ~ dist_cat, data = data.CM)

summary(mod2.2)
summary(mod2.3)

AIC(mod2.2)
AIC(mod2.3)

anova(mod2.2,mod2.3,test = "F") # Adding body size did NOT significantly improve fit (but close); stick with simpler model

myAIC <- round(AIC(mod2.1, mod2.2,mod2.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Cephalopholis argus
mod3.2 <- lm(Mean.1 ~ dist_cat + std_len_cm, data = data.CA)
mod3.3 <- lm(Mean.1 ~ dist_cat, data = data.CA)

summary(mod3.2)
summary(mod3.3)

AIC(mod3.2)
AIC(mod3.3)

anova(mod3.2,mod3.3,test = "F") # Adding body size did NOT significantly improve fit; stick with simpler model

myAIC <- round(AIC(mod3.1, mod3.2,mod3.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Lutjanus bohar
mod4.2 <- lm(Mean.1 ~ dist_cat + std_len_cm, data = data.LB)
mod4.3 <- lm(Mean.1 ~ dist_cat, data = data.LB)

summary(mod4.2)
summary(mod4.3)

AIC(mod4.2)
AIC(mod4.3)

anova(mod4.2,mod4.3,test = "F") # Adding body size did NOT significantly improve fit; stick with simpler model

myAIC <- round(AIC(mod4.1, mod4.2,mod4.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Cephalopholis urodeta
mod5.2 <- lm(Mean.1 ~ dist_cat + std_len_cm, data = data.CU)
mod5.3 <- lm(Mean.1 ~ dist_cat, data = data.CU)

summary(mod5.2)
summary(mod5.3)

AIC(mod5.2)
AIC(mod5.3)

anova(mod5.2,mod5.3,test = "F") # Adding body size DID significantly improve fit

myAIC <- round(AIC(mod5.1, mod5.2,mod5.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]



# Lutjanus fulvus
mod6.2 <- lm(Mean.1 ~ dist_cat + std_len_cm, data = data.LF)
mod6.3 <- lm(Mean.1 ~ dist_cat, data = data.LF)

summary(mod6.2)
summary(mod6.3)

AIC(mod6.2)
AIC(mod6.3)

anova(mod6.2,mod6.3,test = "F") # Adding body size DID significantly improve fit

myAIC <- round(AIC(mod6.1, mod6.2,mod6.3),2)
myAIC$dAIC <- myAIC$AIC - min(myAIC$AIC)
#myAIC$weights <- round(Weights(myAIC$AIC),2)
#myAIC[order(-myAIC$weights),]
myAIC[order(myAIC$dAIC),]

