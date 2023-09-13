
# Robust trophic pathways in carnivorous reef fishes across a major human disturbance gradient

# Authors: Matthew D. Ramirez [1,2,3], Kelton W. McMahon [2], Neil Rooney [4], Rana W. El-Sabaawi [1], and Julia K. Baum [1]
# Institution: [1] Department of Biology, University of Victoria, PO BOX 1700 Station CSC, Victoria, British Columbia, V8W 2Y2, Canada
# [2] Graduate School of Oceanography, University of Rhode Island, 215 South Ferry Road, Narragansett, Rhode Island, 02882, USA
# [3] Department of Biology and Marine Biology, University of North Carolina Wilmington, 601 South College Road, Wilmington, NC, 28403, USA
# [4] School of Environmental Sciences, University of Guelph, 50 Stone Road East, Guelph, Ontario, N1G 2W1, Canada

# Corresponding Author: Matthew D. Ramirez, Email: ramirezmd@uncw.edu

# Script to quantify proportional contributions of carbon source to individual carnivorous reef fish via essential amino acid (EAA) d13C analysis and Bayesian stable isotope mixing modeling

# Code was modified from Stock et al. 2018. Analyzing mixing systems using a new generation of Bayesian tracer mixing models. PeerJ. e5096. doi: 10.7717/peerj.5096
# -- including, vignettes provided at: https://github.com/brianstock/MixSIAR


##############################

### Load necessary packages
require(R2jags)
require(MixSIAR)
require(tidyr)
require(dplyr)

### Set your working directory to folder containing this script.
#setwd("C:/Users/.../Cephalopholis urodeta") # If on a PC
#setwd("/Users/.../Cephalopholis urodeta") # If on a Mac
setwd("~/Desktop/KI_fish/analyses/SIMMs/Cephalopholis urodeta")


## Run the 'output_diagnostics' and 'output_stats' scripts to be able to extract data from a fit MixSIAR model (be sure to set file path accordingly)
source("~/Desktop/KI_fish/analyses/SIMMs/output_diagnostics.R")
source("~/Desktop/KI_fish/analyses/SIMMs/output_stats.R")


### Load mixture, source, and discrimination data (be sure to set file path accordingly)
mix.filename  <- "~/Desktop/KI_fish/analyses/SIMMs/Cephalopholis urodeta/KI_fish_CSIA_C_CU.csv"
source.filename <- "~/Desktop/KI_fish/analyses/SIMMs/Cephalopholis urodeta/KI_sources.csv"
discr.filename <- "~/Desktop/KI_fish/analyses/SIMMs/Cephalopholis urodeta/KI_discrimination.csv"


### Choose output options 
#options(max.print=1000000)
output_options <- list(summary_save = TRUE,                    # Save the summary statistics as a txt file?
                       summary_name = "summary_statistics",    # If yes, specify the base file name (.txt will be appended later)
                       sup_post = TRUE,                       # Suppress posterior density plot output in R?
                       plot_post_save_pdf = TRUE,              # Save posterior density plots as pdfs?
                       plot_post_name = "posterior_density",   # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                       sup_pairs = FALSE,                      # Suppress pairs plot output in R?
                       plot_pairs_save_pdf = TRUE,             # Save pairs plot as pdf?
                       plot_pairs_name = "pairs_plot",         # If yes, specify the base file name (.pdf/.png will be appended later)
                       sup_xy = TRUE,                         # Suppress xy/trace plot output in R?
                       plot_xy_save_pdf = TRUE,                # Save xy/trace plot as pdf?
                       plot_xy_name = "xy_plot",               # If yes, specify the base file name (.pdf/.png will be appended later)
                       gelman = TRUE,                          # Calculate Gelman-Rubin diagnostic test?
                       heidel = FALSE,                         # Calculate Heidelberg-Welch diagnostic test?
                       geweke = TRUE,                          # Calculate Geweke diagnostic test?
                       diag_save = TRUE,                       # Save the diagnostics as a txt file?
                       diag_name = "diagnostics",              # If yes, specify the base file name (.txt will be appended later)
                       indiv_effect = FALSE,                   # Is Individual a random effect in the model? (already specified)
                       plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
                       plot_pairs_save_png = FALSE,            # Save pairs plot as png?
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = FALSE,
                       return_obj = TRUE)


### Define prior and error structure
# "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
alpha.gen <- rep(1, 4)
resid_err <- TRUE
process_err <- TRUE


### Run model

# Define subdirectory to store model results 
mainDir <- getwd()
subDir <- paste0("model")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

mix_CU <- load_mix_data(filename=mix.filename, 
                        iso_names=c("Thr","Val","Ile","Leu","Phe","Lys"), 
                        factors=c("fish_code"), 
                        fac_random=c(TRUE), 
                        fac_nested=c(FALSE), 
                        cont_effects=NULL)

source <- load_source_data(filename=source.filename, source_factors=NULL, 
                           conc_dep=FALSE, data_type="means", mix_CU)

discr <- load_discr_data(filename=discr.filename, mix_CU)


### Write JAGS file
model_filename <- "MixSIAR_mix_CU.txt"
write_JAGS_model(model_filename, resid_err, process_err, mix_CU, source)

# MCMC run set to "test" for demonstration but was run as "very long" in final analyses
jags_CU <- run_model(run="test", mix_CU, source, discr, model_filename, 
                     alpha.prior = alpha.gen, resid_err, process_err)


### Analyze diagnostics and output
output_JAGS(jags_CU, mix_CU, source, output_options)

# Diagnostics output
diag <- output_diagnostics(jags_CU, mix_CU, source, output_options)
names(diag)
head(diag$gelman)
head(diag$geweke)

# Summary statistics 
df.stats <- output_stats(jags_CU, mix_CU, source, output_options)
df.stats <- as.data.frame(df.stats)

# extract data for individual fish
df.stats2 <- df.stats
df.stats2$row_names <- row.names(df.stats2)
df.stats2 <- df.stats2[-c(1:11),]
df.stats2 <- separate(data=df.stats2, col=row_names, into=c("p","fish_code","source"), sep="\\.")
df.stats2 <- df.stats2[,-10]
rownames(df.stats2) <- NULL
colnames(df.stats2)[3:9] <- c("p2.5", "p5", "p25", "p50", "p75", "p95", "p97.5")

# extract data for species
df.stats3 <- df.stats
df.stats3$row_names <- row.names(df.stats3)
df.stats3 <- df.stats3[c(8:11),]
df.stats3 <- separate(data=df.stats3, col=row_names, into=c("p","fish_code","source"), sep="\\.")
df.stats3 <- df.stats3[,-10]
rownames(df.stats3) <- NULL
colnames(df.stats3)[3:9] <- c("p2.5", "p5", "p25", "p50", "p75", "p95", "p97.5")
df.stats3$species_code <- "CU"
df.stats3$sci_name <- "Cephalopholis urodeta"


### Save results
saveRDS(jags_CU, "mix_CU.rds")
#mix_CU <- readRDS("mix_CU.rds")


### Move back up to root directory
setwd(mainDir)


### Append sample metadata to results for individuals
meta <- read.csv("KI_fish_CSIA_C_CU.csv") 
meta <- meta[,c(1:10)]
df.stats2 <- merge(df.stats2, meta, by ="fish_code")


### Save summary statistics to main SIMMs folder
setwd("~/Desktop/KI_fish/analyses/SIMMs")

write.csv(df.stats2, "sum_stat_CU.csv")
write.csv(df.stats3, "sum_stat_CU_global.csv")
