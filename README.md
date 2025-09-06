# Ramirez_etal_2025_JournalOfAnimalEcology

This repository contains data and R code associated with the article "Robust trophic pathways in carnivorous reef fishes across a major human disturbance gradient," published in teh Journal of Animal Ecology.

_Reference:_ Ramirez MD, Besser AC, Newsome SD, McMahon KW (in review) Meta-analysis of primary producer amino acid δ<sup>15</sup>N values and their influence on trophic position estimation. Methods in Ecology and Evolution XX, XXX-XXX. https://doi.org/

### File Descriptions ###

* ***Beta_data.xlsx*** includes β value data, the difference between trophic and source amino acid (AA) δ<sup>15</sup>N values within primary producers, resulting from a meta-analysis of the published primary producer literature. See _Metadata_ tab within this file for full description of the dataset. 

  Note: As of June 1, 2021, "Beta_data.xlsx" excludes unpublished data from Chen et al. 2020 (ice algae, n = 6) and A. C. Besser [arid habitat primary producer; freshwater eukaryotic microalgae (n = 4), grass (n = 7), forb (n = 6), cactus (n = 5), shrub (n = 6)]. Please direct data inquiries to Shaomin Chen (ice algae; shaomin.chen[at]dal.ca) or Alexi C. Besser (arid primary producer; acbesser[at]unm.edu). 

* ***Fig1 - Conceptual Simulation.R*** includes code to reproduce Figure 1B, which is a simulation demonstrating how variation in β values and AA-specific TDFs propagate through a simple hypothetical food chain to influence consumer trophic position estimates.

* ***Fig5 - Sensitivity Analysis.R*** includes code to reproduce Figure 5, which illustrates how consumer trophic position estimates change as a function of variation in mean β values and trophic discrimination factors. 

  Note: This code produces consumer-specific results panels that were compiled by the authors using Adobe Illustrator.


* ### File Descriptions (data folder) ###

* ***ki_bulk_sia_final.csv*** is a comma-delimited file containing sample meta-data and bulk stable isotope data for sampled carnivorous reef fish.

* ***ki_csia_c.csv*** is a comma-delimited file containing sampled meta-data and amino acid-specific stable carbon isotope data for sampled carnivorous reef fish.

* ***ki_csia_n.csv*** is a comma-delimited file containing sampled meta-data and amino acid-specific stable nitrogen isotope data for sampled carnivorous reef fish.

* ***ki_fish_data_sum.csv*** is a comma-delimited file containing summary data for fish observed in underwater visual census surveys (from Magel et al. 2020).

### File Descriptions (analysis folder) ###

* ***KI_bulk_SIA.R*** includes code to characterize variation in reef fish isotopic niche sizes and positions (i.e., Standard Ellipse Areas) via bulk muscle stable carbon (d13C) and nitrogen (d15N) isotope data.

* ***KI_LDA.R*** includes code to classify individual carnivorous reef fish to carbon source groups via essential amino acid (EAA) d13C fingerprinting and simple and bootstrapped linear discriminant analysis (LDA).

* ***KI_trophic_position.R*** includes code to calculate and compare reef fish trophic positions using amino acid d15N values

* ***KI_fish_community_metrics.R*** includes code to create Fig. 1 B-D (multi-panel figure plotting relative biomass and abundance of fish trophic groups at each local human disturbance level)

* ***summarySE.R*** is a function to calculate summary statistics

### File Descriptions (SIMMs sub-folder) ###

* ***KI_SIMM_plots.R includes code to collate and plot SIMM results across species.

* ***output_diagnostics.R includes code to returns diagnostics for a fit MixSIAR model.

* ***output_stats.R includes code to return summary statistics from a fit MixSIAR model.

### Species-specific subfolders: ###

* ***KI_mixsiar_XX.R*** includes code to quantify proportional contributions of carbon sources to individual carnivorous reef fish via essential amino acid (EAA) d13C analysis and Bayesian stable isotope mixing modeling.

* ***KI_discrimination.csv*** is a comma-delimited file containing mean and SDs of discrimination factors used in the SIMMs.

* ***KI_sources.csv*** is a comma-delimited file containing mean and SDs of essential amino acid (EAA) d13C values from the carbon source proxies.

* ***KI_fish_CSIA_C_XX.csv*** is a comma-delimited file containing sample meta-data and amino acid-specific stable carbon isotope data for sampled carnivorous reef fish.

Note: Most figures were further modified using Adobe Illustrator before final publication.

