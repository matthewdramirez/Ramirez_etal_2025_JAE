
# R code for Ramirez et al. 2025 "Carbon pathways and trophic attributes are conserved in carnivorous reef fishes across a major human disturbance gradient" Journal of Animal Ecology

# Script to classify individual carnivorous reef fish to carbon source groups via essential amino acid (EAA) d13C fingerprinting and simple and bootstrapped linear discriminant analysis (LDA)

# Code was modified from Fox et al. 2019. Trophic plasticity in a common reef‐building coral: Insights from δ<sup>13</sup>C analysis of essential amino acids. Functional Ecology, 33, 2203-2214.

# All figures were modified in Adobe Illustrator prior to publication.

##############################

### Linear Discriminant Analysis -----------------------------------------------------------

## Load necessary packages
require(reshape2)
require(ggplot2)
require(MASS)
require(ellipse)
require(dplyr)
require(plyr)


## Set your working directory
# Make sure that this contains the "ki_csia_c.csv" file
#setwd("C:/Users/...") # If on a PC
#setwd("/Users/...") # If on a Mac
setwd("~/Desktop/KI_fish/data")


# Load the AA d13C dataset
fish<-read.csv("ki_csia_c.csv") 

fish$dist_cat <- factor(fish$dist_cat, level = c("Very Low", "Medium", "Very High"))
fish$fish_code <- as.factor(fish$fish_code)
fish$species_code <- as.factor(fish$species_code)
fish$sci_name <- as.factor(fish$sci_name)

fish$site <- factor(fish$site, level = c("10","11","16","15","19","14","8","35","34","27","30"))
# Note site 11 only in AA-C data and site 14 NOT in AA-N data but in AA-C and bulk SIA data
# So, AA-N data has VL = 4 sites, M = 3 sites, and VH = 2 sites

fish$pub.name <- factor(fish$pub.name, level = c("VL6","VL9","VL11","VL1","VL2","M4","M1","M2","M3","VH1","VH3"))
# Note site VL9 only in AA-C data and site M4 NOT in AA-N data but in AA-C and bulk SIA data
# So, AA-N data has VL = 4 sites, M = 3 sites, and VH = 2 sites

# Setup dataframes

# Pull out just the essential amino acid values (Ile, Leu, Lys, Phe, Thr, Val); rename columns
# Exclude Met since missing from some fish
fish.ess <- data.frame(fish$fish_code, fish$group, fish$sci_name, fish$dist_cat, fish$Ile, fish$Leu, 
                       fish$Lys, fish$Phe, fish$Thr, fish$Val) 
colnames(fish.ess) <- c("fish_code", 'group', 'sci_name', 'dist_cat',"Ile", "Leu", "Lys","Phe","Thr", "Val") 
fish.ess$fish_code <- as.character(fish.ess$fish_code)

# Pull out the Non-essential amino acids (Ala, Asp, Gly, Glu, Pro, Ser); rename columns 
# Exclude Arg since missing from some fish
fish.Ness <- data.frame(fish$fish_code, fish$group, fish$sci_name, fish$dist_cat, fish$Ala, fish$Asp, 
                        fish$Gly, fish$Glu, fish$Pro, fish$Ser)
colnames(fish.Ness) <- c("fish_code", 'group', 'sci_name','dist_cat', "Ala", "Asx", "Gly", "Glx", "Pro", "Ser")
fish.Ness$fish_code <- as.character(fish.Ness$fish_code)


### > Plot d13C values for  essential amino acids of carbon sources -----------------------------------------------------------

# Melt AA groups into long-format and combine into a single dataframe
sourceE<-melt(fish.ess, id=c("fish_code","group", "sci_name","dist_cat"))
colnames(sourceE)<-c("fish_code","group", "sci_name","dist_cat","AA","d13C")

sourceN<-melt(fish.Ness, id=c("fish_code","group", "sci_name","dist_cat"))
colnames(sourceN)<-c("fish_code","group", "sci_name","dist_cat","AA","d13C")

aa_plot<-rbind(sourceN,sourceE)

# Adjust labels and re-order amino acids
aa_plot$AA<-ordered(aa_plot$AA,levels=c("Ala","Asx","Gly","Glx","Pro","Ser",
                                        "Ile","Leu","Lys","Phe","Thr","Val"))
aa_plot$group <- as.factor(aa_plot$group)

aa_plot_source <- subset(aa_plot, group == "Co" | group == "De" | group == "He" | group == "Zp" )
aa_plot_source <- subset(aa_plot_source, AA == "Thr" | AA == "Val" | AA == "Leu" | AA == "Ile" | AA == "Phe" | AA == "Lys")

aa_plot_fish <- subset(aa_plot, group == "Pi" | group == "GC")
aa_plot_fish <- subset(aa_plot_fish, AA == "Thr" | AA == "Val" | AA == "Leu" | AA == "Ile" | AA == "Phe" | AA == "Lys")

#min(aa_plot_source$d13C)
#max(aa_plot_source$d13C)

ggplot(aa_plot_source, aes(x=AA,y=d13C,fill=sci_name)) +
  geom_boxplot(outlier.shape = 1) +
  scale_fill_manual(labels = c(expression(paste("Coral (", italic("C. ornatissimus"),")")), 
                               expression(paste("Detritus (", italic("C. marginatus"),")")), 
                               expression(paste("Algae (", italic("A. nigricans"),")")), 
                               expression(paste("Zooplankton (", italic("P. olivaceus"),")"))),
                    values = c("#d01c8b","#a6611a", "#008837", "#2c7bb6")) +
  theme_test() + 
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA)) +
  theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  theme(axis.text = element_text(colour="black")) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  scale_y_continuous(limits=c(-28,-4),breaks = seq(-28,-4,4)) +
  labs(y=expression({delta}^13*C~'(\u2030)'),x=expression("Amino Acid")) +
  #theme(legend.position = 'right') +
  theme(legend.text.align = 0) +
  guides(fill = guide_legend(title = "Carbon Source", title.position = "top", title.hjust = 0.5))


### > Simple LDA -----------------------------------------------------------

# LDA with Co, De, He, Zp as sources and classify the carnivorous fish fractions

# Separate out carbon sources (Co, De, He, Zp) and carnivorous fish
fish.source <- fish.ess
fish.source <- fish.source[!(fish.source$group == "Pi" | fish.source$group == "GC"),]

fish.animal <- fish.ess[(fish.ess$group == "Pi" | fish.ess$group == "GC"),]
fish.animal$group <- NULL

## Calculate an error rate for our LDA model - using carbon source dataset
fish.lda <- lda(group ~ Ile + Leu + Lys + Phe + Thr + Val, data = fish.source, CV = TRUE)

# Create a table which compares the classification of the LDA model to the actual spp
ct.prod <- table(fish.source$group, fish.lda$class)

# Percent of samples correctly classified
sum(diag(prop.table(ct.prod)))  # total
diag(prop.table(ct.prod, 1))  # by source group

# Create a training lda function from the carbon source data - we will use this to classify the carnivorous fish.
# Examine coefficents of linear discriminants to determine AAs contributing to groups seperation.
# (Supplementary Table 2)
(fish.train <- lda(group ~ Ile + Leu + Lys + Phe + Thr + Val, data = fish.source))

# Create a dataframe with these LDA coordinates
datPred <- data.frame(fish_code = fish.source$fish_code, group=fish.source$group, sci_name = fish.source$sci_name, dist_cat = fish.source$dist_cat, predict(fish.train)$x) #create data.frame

# Predict the carnivorous fish fractions based on ESS
animal.res <- predict(fish.train, fish.animal)
class <- cbind(animal.res$class, animal.res$posterior, animal.res$x)
#write.csv(class,"classifications_allEAA.csv")

# Save the predicted coordinates 
datPred2 <- data.frame(fish_code = fish.animal$fish_code, dist_cat = fish.animal$dist_cat, sci_name = fish.animal$sci_name, group='Pi',
                       animal.res$x)

#write.csv(datPred2, "LDA_pred.csv")

# Add to the original animal dataframe which individual got classified as what
# (Supplementary Table 1)
fish.animal$class <- animal.res$class 

# Merge the source and animal dataframes for plotting
datPred3 <- rbind(datPred , datPred2)
colnames(datPred3)[1] <- "fish_code"



### > > Simple LDA plots (Figure 2) -----------------------------------------------------------

# LD1 vs. LD2
# LD2 vs. LD3 and LD1 vs. LD3 not reported given low contribution of LD3 -- plots also look like nonsense

dat_ell <- data.frame() 

datPred$group <- as.factor(datPred$group)
#datPred$group <- factor(datPred$group, level = c("He", "Co", "De", "Zp"))
datPred$group <- factor(datPred$group, level = c("Co", "De", "He", "Zp"))
datPred$group <- factor(datPred$group, level = c("Co", "De", "He", "Zp"))

for(g in levels(datPred$group)){ 
  dat_ell <- rbind(dat_ell, cbind(as.data.frame(with(datPred[datPred$group==g,], 
                                                     ellipse(cor(LD1, LD2), 
                                                             scale=c(sd(LD1), sd(LD2)), 
                                                             centre=c(mean(LD1), mean(LD2))))), group=g)) } 

#levels(datPred3$sci_name)
datPred3$sci_name <- factor(datPred3$sci_name, level = c("Chaetodon ornatissimus","Ctenochaetus marginatus","Acanthurus nigricans","Pseudanthias olivaceus",
                                                 "Aphareus furca", "Caranx melampygus","Cephalopholis argus",  "Lutjanus bohar", "Cephalopholis urodeta", "Lutjanus fulvus"))



ggplot(datPred3, aes(x=LD1, y=LD2, shape = group) ) + 
  geom_point(data=subset(datPred3, group == "Pi" | group == "GC"),  aes(col = sci_name), size=5, alpha = 0.85) + #aes(color = sci_name, shape = sci_name)
  geom_path(data=subset(dat_ell, group != "Pi" | group == "GC"), aes(x=x, y=y), size =0.75, linetype = 2) + 
  geom_point(data=subset(datPred3, group != "Pi" | group == "GC"),  size=4) + #aes(color = sci_name, shape = sci_name)
  scale_shape_manual(values = c(0,2,5,19,6))+ # not sure why, buy regardless of left, GC and Pi have same symbol within plot
  scale_color_manual(values = c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677","#AA3377"))+ # not sure why, buy regardless of left, GC and Pi have same symbol within plot
  theme_test() +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=14)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(plot.title.position = 'plot') +
  theme(legend.position = 'top') + 
  guides(colour = guide_legend(title = "Species", title.position = "top"), 
         shape = guide_legend(title = "Carbon Source", title.position = "top")) +
  xlab(expression(paste(LD[1]~"(83.47%)"))) +
  ylab(expression(paste(LD[2]~"(13.92%)"))) +
  coord_cartesian(clip = 'off') 


ggplot(datPred3, aes(x=LD1, y=LD2, shape = group) ) + 
  geom_path(data=subset(dat_ell, group != "Pi" | group == "GC"), aes(x=x, y=y), size =0.75, linetype = 2) + 
  geom_point(data=subset(datPred3, group != "Pi" | group == "GC"),  size=4) + #aes(color = sci_name, shape = sci_name)
  geom_point(data=subset(datPred3, group == "Pi" | group == "GC"),  aes(col = sci_name), size=5, alpha = 0.85) + #aes(color = sci_name, shape = sci_name)
  scale_shape_manual(values = c(0,2,5,19,6))+ # not sure why, buy regardless of left, GC and Pi have same symbol within plot
  scale_color_manual(values = c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677","#AA3377"))+ # not sure why, buy regardless of left, GC and Pi have same symbol within plot
  theme_test() +
  theme(axis.text = element_text(colour="black")) +
  theme(text=element_text(size=14)) +
  theme(axis.ticks.length=unit(2.5,"mm")) +
  theme(plot.title.position = 'plot') +
  theme(legend.position = 'top') + 
  guides(colour = guide_legend(title = "Species", title.position = "top"), 
         shape = guide_legend(title = "Carbon Source", title.position = "top")) +
  xlab(expression(paste(LD[1]~"(83.47%)"))) +
  ylab(expression(paste(LD[2]~"(13.92%)"))) +
  coord_cartesian(clip = 'off') + 
  facet_wrap(~dist_cat)


### > Bootstrapped LDA (Table S2) -----------------------------------------------------------

# This part of the code calls on outputs from the original (raw data) LDA above.
# Specifically, the projected LDA coordinates for the carnivorous fish (datPred2) are required for the permutational distance analysis. 
# This variable is re-written as "poc" below.

# Isolate carbon sources
He <- subset(fish.source,group=="He")
Co <- subset(fish.source,group=="Co")
De <- subset(fish.source,group=="De")
Zp <- subset(fish.source,group=="Zp")

# Create blank objects to fill with data
He.boot <- NULL
Co.boot <- NULL
De.boot <- NULL
Zp.boot <- NULL
source.boot <- NULL
animal.class <- NULL
a <- NULL
b <- NULL
animal.prop <- NULL
prop <- NULL
source.prop <- NULL

# Set number of permutation
# 10 is shown here as an example to make the code run faster, but 10,000 was used in the manuscript
n <- 1000  

# Permutation analysis
for(i in 1:n){  
  
  # Randomly sample each of the 4 source populations with replacement for the number of samples for each group
  He.boot[[i]]<-sample_n(He,15,replace=TRUE)
  He.boot<-as.data.frame(He.boot[[i]])
  
  Co.boot[[i]]<-sample_n(Co,15,replace=TRUE)
  Co.boot<-as.data.frame(Co.boot[[i]])
  
  De.boot[[i]]<-sample_n(De,15,replace=TRUE)
  De.boot<-as.data.frame(De.boot[[i]])
  
  Zp.boot[[i]]<-sample_n(Zp,15,replace=TRUE)
  Zp.boot<-as.data.frame(Zp.boot[[i]])
  
  # Combine the newly resampled source data
  source.boot<-rbind(He.boot,Co.boot,De.boot,Zp.boot)
  
  ## Calculating an error rate for our LDA model - using source dataset
  # first... lets run an LDA with a jacknifing model fit, to look at error rate
  # in this function the line 'CV = TRUE' makes the LDA do jacknifed (leave one out) model fit
  fish.lda <- lda(group ~ Ile + Leu + Lys + Phe + Thr + Val, data = source.boot, CV = TRUE)
  
  # Create a table which compares the classification of the LDA model to the actual spp
  ct.prod <- table(source.boot$group, fish.lda$class)
  
  # Report what % of each species is being correctly classified; compile over runs
  prop<-diag(prop.table(ct.prod, 1))
  source.prop<-rbind(source.prop,prop)
  
  # Create a training lda function from the source data - we will use this to classify the carnivorous fish
  fish.train <- lda(group ~ Ile + Leu + Lys + Phe + Thr + Val, data = source.boot)
  
  # Create a dataframe with these LDA coordinates
  datPred <- data.frame(fish_code = source.boot$fish_code, group=source.boot$group, sci_name = source.boot$sci_name, dist_cat = source.boot$dist_cat, predict(fish.train)$x) #create data.frame
  
  # Predict the carnivorous fish fractions based on ESS... and save this into a dataframe
  animal.res <- predict(fish.train, fish.animal)
  
  # Extract classifications into a table 
  a <- table(animal.res$class) #this calculates percentages of classification
  
  # pair each individual with their classification
  # in b 4=Zp, 3=De, 2=Co, 1=He
  b <- as.data.frame(cbind(fish.animal$fish_code,animal.res$class))
  
  # Build a df that contains all the successive classifications 
  animal.prop<-as.data.frame(rbind(animal.prop,a))
  
  # Compile b so we can count how many times an individual was classified as a particular 
  # group
  animal.class<-as.data.frame(rbind(animal.class,b))
}


# Calculate percentages of classification for each individual
# (Supplementary Table 1)
id_prop <- ddply(animal.class, c("V1"), summarise,
               Co = (sum(V2=="1")/n),
               De = (sum(V2=="2")/n),
               He = (sum(V2=="3")/n),
               Zp = (sum(V2=="4")/n))

#write.csv(id_prop, "LDA_boot.csv")

# Calculate percentages of classification based on total sample size
animal.prop$Co.prop <- animal.prop$Co / nrow(fish.animal)
animal.prop$De.prop <- animal.prop$De / nrow(fish.animal)
animal.prop$He.prop <- animal.prop$He / nrow(fish.animal)
animal.prop$Zp.prop <- animal.prop$Zp / nrow(fish.animal)

# Calculate mean and CI for each carbon source
(LDA_boot_sum <- sapply(animal.prop, function(animal.prop) 
  c("Mean"= mean(animal.prop,na.rm=TRUE),
    "Median" = median(animal.prop,na.rm=TRUE),
    "CI bound" = (sd(animal.prop)*1.96),
    "Upper CI" = quantile(animal.prop,0.975),
    "Lower CI" = quantile(animal.prop,0.025))))

#write.csv(LDA_boot_sum, "LDA_boot_sum.csv")