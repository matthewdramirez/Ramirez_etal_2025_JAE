# Ramirez_etal_KIFishTrophic

This repository contains data and R code associated with the article "Robust trophic pathways in carnivorous reef fishes across a major human disturbance gradient," currently aunder review.

_Reference:_ Ramirez, MD, KW McMahon KW, N Rooney, R W. E-S, and JK Baum  (in review) Robust trophic pathways in carnivorous reef fishes across a major human disturbance gradient.

### File Descriptions (data folder) ###

* ***ki_bulk_sia_final.csv*** is a comma-delimited file containing sample meta-data and bulk stable isotope data for sampled carnivorous reef fish.

* ***ki_csia_c.csv*** is a comma-delimited file containing sampled meta-data and amino acid-specific stable carbon isotope data for sampled carnivorous reef fish.

* ***ki_csia_n.csv*** is a comma-delimited file containing sampled meta-data and amino acid-specific stable nitrogen isotope data for sampled carnivorous reef fish.

* ***ki_fish_data_sum.csv*** is a comma-delimited file containing summary data for fish observed in underwater visual census surveys (from Magel et al. 2020).

### File Descriptions (analysis folder) ###
* ***KI_bulk_SIA.R*** includes code to characterize variation in reef fish isotopic niche sizes and positions (i.e., Standard Ellipse Areas) via bulk muscle stable carbon (d13C) and nitrogen (d15N) isotope data.

* ***KI_LDA.R*** includes code to classify individual carnivorous reef fish to carbon source groups via essential amino acid (EAA) d13C fingerprinting and simple and bootstrapped linear discriminant analysis (LDA).

* ***KI_trophic_position.R*** includes code to calculate and compare reef fish trophic positions using amino acid d15N values

* ***KI_fish_community_metrics.R*** includes code to create Fig. 1 B-D (multi-panel figure plotting relative biomass and abundance of fish trophic groups at each local human disturbance level)

* ***summarySE.R*** is a function to calculate summary statistics

### File Descriptions (SIMMs sub-folder) ###
* ***KI_SIMM_plots.R*** includes code to collate and plot SIMM results across species.

* ***output_diagnostics.R*** includes code to returns diagnostics for a fit MixSIAR model.

* ***output_stats.R*** includes code to return summary statistics from a fit MixSIAR model.

* Species-specific subfolders:
  * ***KI_mixsiar_XX.R*** includes code to quantify proportional contributions of carbon sources to individual carnivorous reef fish via essential amino acid (EAA) d13C analysis and Bayesian stable isotope mixing modeling.

  * ***KI_discrimination.csv*** is a comma-delimited file containing mean and SDs of discrimination factors used in the SIMMs.

  * ***KI_sources.csv*** is a comma-delimited file containing mean and SDs of essential amino acid (EAA) d13C values from the carbon source proxies.

  * ***KI_fish_CSIA_C_XX.csv*** is a comma-delimited file containing sample meta-data and amino acid-specific stable carbon isotope data for sampled carnivorous reef fish.


Note: Most figures were further modified using Adobe Illustrator before final publication.
