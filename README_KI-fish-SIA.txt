<!---
This README uses Markdown syntax, though it is saved in .txt format.
To view the file with its intended text formatting, save a copy as a .md file.
Then you can open it with the browser or text viewer of your choice.
For details on Markdown syntax, see <https://daringfireball.net/projects/markdown/>.
--->

Reference Information
=====================

Provenance for this README
--------------------------

* File name: README_KI-fish-SIA.txt
* Authors: Matthew D. Ramirez, Kelton W. McMahon, Neil Rooney, Rana W. El-Sabaawi, Julia K. Baum
* Other contributors: REDACTED
* Date created: 2023-09-07
* Date modified: 2025-08-29

Dataset Version and Release History
-----------------------------------

* Current Version:
  * Number: 4.0.0
  * Date: 2024-07-31
  * Persistent identifier: n/a
  * Summary of changes: n/a

* Embargo Provenance: n/a
  * Scope of embargo: n/a
  * Embargo period: n/a

Dataset Attribution and Usage
-----------------------------

* Dataset Title: Data for the article "Carbon pathways and trophic attributes are conserved in carnivorous reef fishes across a major human disturbance gradient."

* Persistent Identifier: n/a

* Dataset Contributors:

  * Creators: Matthew D. Ramirez, Kelton W. McMahon, Neil Rooney, Rana W. El-Sabaawi, Julia K. Baum

* Date of Issue: 2025-08-29

* Publisher: Journal of Animal Ecology

* License: Use of these data is covered by the following license:
  * Title: CC0 1.0 Universal (CC0 1.0)
  * Specification: n/a

* Dataset citation:
    > Ramirez, M.D., K.W. McMahon, N. Rooney, R.W. El-Sabaawi, and J.K. Baum. 2025. Data for the article "Carbon pathways and trophic attributes are conserved in carnivorous reef fishes across a major human disturbance gradient", Zenodo, Dataset, https://doi.org/10.5281/zenodo.13948016

  * Corresponding publication:
    > Ramirez, M.D., K.W. McMahon, N. Rooney, R.W. El-Sabaawi, and J.K. Baum. 2025. Robust trophic pathways in carnivorous reef fishes across a major human disturbance gradient. Journal of Animal Ecology. 


Contact Information
-------------------

  * Name: Matthew D. Ramirez
  * Affiliations: Department of Biology and Marine Biology, University of North Carolina Wilmington; Department of Biology, University of Victoria; Graduate School of Oceanography, University of Rhode Island 
  * ORCID ID: https://orcid.org/0000-0002-9628-8517
  * Email: ramirezmd@uncw.edu
  * Address: 5600 Marvin K. Moss Lane, Wilmington, NC, 28409

* Alternative Contact: 
  * Name: Julia K. Baum
  * Affiliations: Department of Biology, University of Victoria
  * ORCID ID: https://orcid.org/0000-0002-9827-1612
  * Email: baum@uvic.ca
  * Address: PO BOX 1700 Station CSC, Victoria, British Columbia, V8W 2Y2, Canada

* Contributor ORCID IDs:
  * Kelton W. McMahon: https://orcid.org/0000-0002-9648-4614
  * Neil Rooney: https://orcid.org/0000-0003-4360-4034
  * Rana W. El-Sabaawi: https://orcid.org/0000-0002-0561-1068
- - -

Additional Dataset Metadata
===========================

Dates and Locations
-------------------

* Dates of data collection: Field data and samples were July–August of 2010–2012 on Kiritimati Atoll. Bulk SIA analyses were conducted between 203 and 2020. CSIA-AA analyses were conducted in 2022.

* Geographic locations of data collection: Kiritimati Atoll. See Figure 1 in manuscript.

* Other locations pertaining to dataset contents: Bulk SIA analyses were performed at the University of Victoria and University of Windsor. CSIA-AA analyses were performed at the University of Rhode Island.

- - -


Methodological Information
==========================

* Methods of data collection/generation: see manuscript for details

- - -

Data and File Overview
======================

Summary Metrics
---------------

* File count: 33
* Total file size: 61 KB
* Range of individual file sizes: 16 - 25 KB
* File formats: .csv, .R

Naming Conventions
------------------

* File naming scheme: prefix denotes data types. 

Table of Contents
-----------------

* data (folder)
  * ki_bulk_sia_final.csv
  * ki_csia_c.csv
  * ki_csia_n.csv
  * ki_fish_data_sum.csv

* analyses (folder)
  * KI_bulk_SIA.R
  * KI_LDA.R
  * KI_trophic_position.R
  * KI_fish_community_metrics.R
  * summarySE.R
  * SIMMs (folder)
    * KI_SIMM_plots.R
    * output_diagnostics.R
    * output_stats.R
    * Cephalopholis argus (folder)
      * KI_mixsiar_CA.R
      * KI_discrimination.csv
      * KI_sources.csv
      * KI_fish_CSIA_C_CA.csv
    * Cephalopholis urodeta (folder)
      * KI_mixsiar_CU.R
      * KI_discrimination.csv
      * KI_sources.csv
      * KI_fish_CSIA_C_CU.csv
    * Aphareus furca (folder)
      * KI_mixsiar_AF.R
      * KI_discrimination.csv
      * KI_sources.csv
      * KI_fish_CSIA_C_AF.csv
    * Caranx melampygus (folder)
      * KI_mixsiar_CM.R
      * KI_discrimination.csv
      * KI_sources.csv
      * KI_fish_CSIA_C_CM.csv
    * Lutjanus bohar (folder)
      * KI_mixsiar_LB.R
      * KI_discrimination.csv
      * KI_sources.csv
      * KI_fish_CSIA_C_LB.csv
    * Lutjanus fulvus (folder) 
      * KI_mixsiar_LF.R
      * KI_discrimination.csv
      * KI_sources.csv
      * KI_fish_CSIA_C_LF.csv   

Setup
-----

* Unpacking instructions: n/a

* Relationships between files/folders: n/a

* Recommended software/tools: R

- - -

File/Folder Details
===================

Details for: ki_bulk_sia.csv
---------------------------------------

* Description: a comma-delimited file containing sample meta-data and bulk stable isotope data for sampled carnivorous reef fish.

* Format(s): .csv

* Size(s): 20 KB

* Dimensions: 233 rows x 13 columns

* Variables:
  * id: alphanumeric unique identification code for an individual
  * species: identification code for species 
  * sci_name: scientific name
  * group: fish trophic group (Pi = piscivore, GC = generalist carnivore)
  * year: calendar year of sample collection
  * site: identification code for sampling location
  * dist_cat: local human disturbance level category (see manuscript for details)
  * weight_.g.: fish weight in grams
  * total_length_.mm.: total fish length in millimeters
  * standard_length_.mm.: standard fish length in millimeters
  * d13C (‰): bulk stable carbon isotope values; permil
  * d15N (‰): bulk stable nitrogen isotope value; permil
  * C.N: ratio of percent carbon to percent nitrogen in bulk SIA sample

* Missing data codes: NA


Details for: ki_csia_c.csv
---------------------------------------

* Description: a comma-delimited file containing sampled meta-data and amino acid-specific stable carbon isotope data for sampled carnivorous reef fish.

* Format(s): .csv

* Size(s): 37 KB

* Dimensions: 145 rows x 37 columns

* Variables:
  * fish_code: condensed alphanumeric unique identification code for an individual
  * species_code: condensed identification code for species 
  * id: alphanumeric unique identification code for an individual
  * sci_name: scientific name
  * group: fish trophic group (Pi = piscivore, GC = generalist carnivore, He = herbivore, De = detritivore, Co = corallivore, Zp = zooplanktivore)
  * site: identification code for sampling location
  * dist_cat: local human disturbance level category (see manuscript for details)
  * isotope: stable isotope measured
  * n_inj: number of sample injections on GC-C-IRMS per sample
  * Columns 10 through 37: amino acid-specific stable carbon isotope values for sampled fish. Three letter codes denote each amino acid: Ala = alanine, Gly = glycine, Thr = threonine, Ser = serine, Val = valine, Leu = leucine, Ile = isoleucine, Pro = proline, Asp = aspartic acid, Glu = glutamic acid, Phe = phenylalanine, Lys = lysine. Values are means of triplicate injections. Columns with the suffix "SD" denote the standard deviation of the triplicate injections for each isotope and amino acid.

* Missing data codes: NA


Details for: ki_csia_n.csv
---------------------------------------

* Description: a comma-delimited file containing sampled meta-data and amino acid-specific stable nitrogen isotope data for sampled carnivorous reef fish.

* Format(s): .csv

* Size(s): 41 KB

* Dimensions: 85 rows x 46 columns

* Variables:
  * fish_code: condensed alphanumeric unique identification code for an individual
  * species_code: condensed identification code for species 
  * id: alphanumeric unique identification code for an individual
  * sci_name: scientific name
  * group: fish trophic group (Pi = piscivore, GC = generalist carnivore, He = herbivore, De = detritivore, Co = corallivore, Zp = zooplanktivore)
  * site: identification code for sampling location
  * dist_cat: local human disturbance level category (see manuscript for details)
  * pub.name: alternative site name
  * continous.pressure.2km: quantitative human pressure value (see manuscript for details)
  * wgt_g: fish weight in grams
  * tot_len_mm: total fish length in millimeters
  * std_len_mm: standard fish length in millimeters
  * isotope: stable isotope measured
  * n_inj: number of sample injections on GC-C-IRMS per sample
  * sumV = microbial resynthesis 
  * Columns 15 through 46: amino acid-specific stable nitrogen isotope values for sampled fish. Three letter codes denote each amino acid: Ala = alanine, Gly = glycine, Thr = threonine, Ser = serine, Val = valine, Leu = leucine, Ile = isoleucine, Pro = proline, Asp = aspartic acid, Glu = glutamic acid, Phe = phenylalanine, Lys = lysine. Values are means of triplicate injections. Columns with the suffix "SD" denote the standard deviation of the triplicate injections for each isotope and amino acid.

* Missing data codes: NA



Details for: ki_fish_data_sum.csv
---------------------------------------

* Description: a comma-delimited file containing summary data for fish observed in underwater visual census surveys (from Magel et al. 2020)

* Format(s): .csv

* Size(s): 12 KB

* Dimensions: 66 rows x 23 columns

* Variables:
  * year: calendar year of sample collection
  * date: date in day_month_year format
  * site: identification code for sampling location
  * observer: initials of the observer who conducted the survey
  * dist_cat: local human disturbance level category (see manuscript for details)
  * Columns 6 through 23: biomass (BM) and abundance (AB) estimates for all fish (total) and by trophic group (corallivores = 'coral', detritivores = 'det', generalist carnivores = 'gen', herbivores = 'herb', invertivores = 'inv', omnivores = 'omn', piscivores = 'pisc', planktivores = 'plank')

* Missing data codes: NA



Details for: R scripts in 'analyses' folder
---------------------------------------

* Description: R scripts for completing statistical analyses.

* Format(s): .R

* R scripts:
  * KI_bulk_SIA.R: Code to characterize variation in reef fish isotopic niche sizes and positions (i.e., Standard Ellipse Areas) via bulk muscle stable carbon (d13C) and nitrogen (d15N) isotope data. 
  * KI_LDA.R: Code to classify individual carnivorous reef fish to carbon source groups via essential amino acid (EAA) d13C fingerprinting and simple and bootstrapped linear discriminant analysis (LDA)
  * KI_trophic_position.R: Code to calculate and compare reef fish trophic positions using amino acid d15N values
  * KI_fish_community_metrics.R: Code to create Fig. 1 B-D (multi-panel figure plotting relative biomass and abundance of fish trophic groups at each local human disturbance level)
  * summarySE.R: Function to calculate summary statistics



Details for: 'SIMMs' sub-folder within the 'analyses' folder
---------------------------------------

* Description: R scripts and data inputs for implementing and collating results from species-specific stable isotope mixing models (SIMMs).

* Format(s): .R, .csv

* R scripts:
  * KI_SIMM_plots.R: Code to collate and plot SIMM results across species.
  * output_diagnostics.R: Code to returns diagnostics for a fit MixSIAR model
  * output_stats.R: Code to return summary statistics from a fit MixSIAR model

* Species-specific subfolders:
  * KI_mixsiar_XX.R: Code to quantify proportional contributions of carbon sources to individual carnivorous reef fish via essential amino acid (EAA) d13C analysis and Bayesian stable isotope mixing modeling
  * KI_fish_CSIA_C_XX.csv: a comma-delimited file containing sample meta-data and amino acid-specific stable carbon isotope data for sampled carnivorous reef fish.
    * Variables:
      * fish_code: condensed alphanumeric unique identification code for an individual
      * species_code: condensed identification code for species 
      * id: alphanumeric unique identification code for an individual
      * sci_name: scientific name
      * group: fish trophic group (Pi = piscivore, GC = generalist carnivore, He = herbivore, De = detritivore, Co = corallivore, Zp = zooplanktivore)
      * site: identification code for sampling location
      * dist_cat: local human disturbance level category (see manuscript for details)
      * wgt_g: fish weight in grams
      * tot_len_mm: total fish length in millimeters
      * std_len_mm: standard fish length in millimeters
      * isotope: stable isotope measured
      * n_inj: number of sample injections on GC-C-IRMS per sample
      * Columns 13 through 39: amino acid-specific stable carbon isotope values for sampled growth layers. Three letter codes denote each amino acid: Ala = alanine, Gly = glycine, Thr = threonine, Ser = serine, Val = valine, Leu = leucine, Ile = isoleucine, Pro = proline, Asp = aspartic acid, Glu = glutamic acid, Phe = phenylalanine, Lys = lysine. Values are means of triplicate injections. Columns with the suffix "SD" denote the standard deviation of the triplicate injections for each isotope and amino acid.
  * KI_discrimination.csv: a comma-delimited file containing mean and SDs ofdiscrimination factors used in the SIMMs.
    * Variables:
      * Source: code for fish trophic group used as a carbon source proxy (He = herbivore, De = detritivore, Co = corallivore, Zp = zooplanktivore)
      * Columns 2 through 13: Mean and SDs for amino acid-specific discrimination factors used in SIMMs. Three letter codes denote each amino acid: Thr = threonine, Val = valine, Leu = leucine, Ile = isoleucine, Phe = phenylalanine, Lys = lysine. 
  * KI_sources.csv: a comma-delimited file containing mean and SDs of essential amino acid (EAA) d13C values from the carbon source proxies.
    * Variables:
      * Source: code for fish trophic group used as a carbon source proxy (He = herbivore, De = detritivore, Co = corallivore, Zp = zooplanktivore)
      * Columns 2 through 13: Mean and SDs for amino acid-specific discrimination factors used in SIMMs. Three letter codes denote each amino acid: Thr = threonine, Val = valine, Leu = leucine, Ile = isoleucine, Phe = phenylalanine, Lys = lysine. 
      * n: sample size


- - -
END OF README
