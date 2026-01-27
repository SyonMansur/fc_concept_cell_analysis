# Frontal Cortex Concept Cell Analysis
by: Syon Mansur

**Date:** UPDATED 2026-01-26

---

## Overview
This repository contains a MATLAB pipeline designed to identify and analyze "concept cells" within the Frontal Cortex (FC). While these units are well-documented in the Medial Temporal Lobe (MTL), this project focuses on a novel identification and characterization of similar selective responses in frontal regions using the ZETA-test protocol.

---

## Pipeline Architecture
The analysis follows a five-step execution flow:

* **1. Extraction (`extractSelectiveResponses.m`)**: Crawls subject directories to load `clusterCharacteristics.mat` files (pre-processed session data) and executes ZETA-tests to identify units with significant stimulus-evoked responses.
* **2. Regional Mapping (`mapBrainRegions.m`)**: Standardizes raw cluster region codes into cleaned anatomical labels and broad area categories (FC vs. MTL).
* **3. Channel Filtering (`filterBestClusters.m`)**: Selects the "best" cluster per electrode channel based on the highest average ZETA score.
* **4. Anatomical Verification (`verifyLocalizations.m`)**: Cross-references unit locations against BIDS-derivative clinical localization Excel sheets for verified anatomical ground truth.
* **5. Visualization (`generatePipelineFigures.m`)**: Produces the final suite of histograms, bar plots, and population statistics comparing selective response dynamics across brain regions.

---

## Methods

### Selectivity Detection (ZETA-test)
Unit selectivity is determined using the ZETA-test, a parameter-free method that detects stimulus-evoked activity by comparing the observed spike-time distribution against a null-hypothesis linear distribution. This allows for robust detection of selective cells without the need for manual binning or arbitrary thresholds.
For more information about the ZETA test, read the original paper here by Montjin et al. (https://elifesciences.org/articles/71969)
### Anatomical Targeting
The project specifically contrasts recorded units in the **Frontal Cortex** (including Anterior Cingulate and Orbitofrontal areas) against the **Medial Temporal Lobe** (Hippocampus, Amygdala, Entorhinal) to identify potential concept cell populations outside of traditional memory centers.

---

## Dependencies
This project requires the following third-party toolbox:
* **ZETA Toolbox**: [Download from GitHub](https://github.com/JorritMontijn/ZETA).

---

## Usage
1. Configure `raw_data_path` and `localization_dir` in `main_pipeline.m`.
2. Add the ZETA toolbox to your MATLAB path.
3. Run `main_pipeline.m`.
