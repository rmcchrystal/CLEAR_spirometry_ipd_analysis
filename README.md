# Unsupervised home spirometry versus supervised clinic spirometry: Analysis of longitudinal participant-level data from the CLEAR trial in patients with bronchiectasis
Public repository containing analysis materials for "Unsupervised home spirometry versus supervised clinic spirometry: Analysis of longitudinal participant-level data from the CLEAR trial in patients with bronchiectasis"

## Overview
We have provided the analysis pipeline for the analysis, comprising R scripts written to prepare, analyse and summarise spirometry data. We are unable to provide the individual participant data used for the analysis, and as such the analysis with current materials cannot be replicated externally. We share these materials solely for the transparency of our analysis setup and for the availability of outputs that formed the main and supplementary results.

## Materials
**Scripts** - R scripts written for data preparation and analysis
- 00_config.R: Configuration for workflow (i.e. clear memory, turn off scientific number notations)
- 00_functions_regex.R: Functions written for analysis pipeline and regex patterns for string detection
- 01_prepare_patients.R: Prepares participant-level metadata (IPD unavailable publicly, cannot be replicated)
- 02_prepare_home.R: Prepares home spirometry data
- 03_prepare_clinic.R: Prepares clinic spirometry data
- 04_join_home.R: Harmonise home and clinic spirometry into a single dataset
- 05_create_denom.R: Updated participant-level metadata
- 06_pex_tidy.R: Prepares exacerbation data
- 07_join.R: Harmonises data sources (metadata, spirometry, exacerbation)
- 08_add_lf_proxies.R: Creates proxy exacerbation dates (i.e. 7d/3d rules)
- 09_prep_data_obj1_2.R: Prepares data for objectives 1 and 2 (agreement and variation over time)
- 10a_obj1_sum.R: Summary statistics for agreement analysis
- 10b_obj1_spiro_cor.R: Correlation summaries for home and clinic spirometry
- 10c_obj1_ba.R: Bland-Altman analysis
- 10e_obj1_ba_sum_tidy.R: Tidy summaries of Bland-Altman analysis
- 11a_obj2.R: Variation over time analysis
- 11b_obj2_sum_tidy.R: Tidy summaries for variation analysis
- 12_obj3.R: Weekly home spirometry analysis
- 13_obj4.R: Spirometry quality analysis
- 14_obj5_3d.R: Exacerbation responsiveness analysis for 3d proxy rule
- 14_obj5_7d.R: Exacerbation responsiveness analysis for 7d proxy rule

**Outputs** - Analysis outputs, including plots and summary statistics
- 01_agreement: Agreement metrics (plots and summaries of lung function, correlation and Bland-Altman analysis)
- 02_variation: Variation over time (plots and summaries of lung function)
- 03_adherence: Weekly home spirometry adherence summaries
- 04_quality: Spirometry quality summaries
- 05_responsiveness: Exacerbation responsiveness summaries
