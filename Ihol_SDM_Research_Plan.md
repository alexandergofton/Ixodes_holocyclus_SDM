# Research, Development, and Implementation Plan
# Species Distribution Model for *Ixodes holocyclus* in Eastern Australia

**Author:** Alexander W. Gofton
**Date:** 2026-02-17
**Version:** 2.0

---

## 1. Background and Rationale

*Ixodes holocyclus* (the eastern paralysis tick) is the most medically significant tick in Australia, responsible for tick paralysis in companion animals and humans, mammalian meat allergy (alpha-gal syndrome), and transmission of several pathogens. Its distribution spans the eastern coastal fringe from north Queensland to eastern Victoria, constrained by climatic factors—particularly temperature, moisture, and rainfall.

### 1.1 Key Findings from the Literature

**Teo et al. (2021) — Climatic requirements and distribution modelling:**

- Compared CLIMEX process-based modelling with a simpler "climatic-range method" (based on Shelford's Law of Tolerance) for mapping the fundamental niche of *I. holocyclus*.
- The climatic-range method explained 96% of the geographic distribution using monthly min/max temperature, rainfall, and relative humidity—outperforming CLIMEX (93%).
- Key physiological thresholds from Heath (1974): lower developmental threshold 8–13°C, optimum 25–28°C, upper limit ~32°C.
- CLIMEX parameters (Table 1): moisture index limits (SM0=0.2, SM1=0.5, SM2=1.5, SM3=2.5), temperature limits (DV0=8, DV1=15, DV2=28, DV3=32), cold stress threshold (-10°C), heat stress threshold (33°C), dry stress threshold (0.1), wet stress threshold (2.5).
- The climatic-range method uses the observed min/max range of each climate variable at known occurrence points to define a climatic envelope—cells where ALL variables fall within the envelope are classified as suitable.
- Under climate change scenarios, the species' range may shift southward toward Melbourne or contract depending on moisture changes.

**Teo et al. (2023) — Phenology and seasonality of tick paralysis:**

- Analysed 22,840 tick paralysis cases in companion animals across 4 Australian regions over 22 years.
- Season timing varies with latitude: earlier onset in the north (QLD) and later in the south (VIC).
- Victoria has a bimodal pattern (autumn: Feb–May and spring: Sep–Dec); Queensland and NSW have a single primary season.
- An 18°C mean daily temperature threshold appears to trigger the tick paralysis season.
- Peak incidence occurs approximately 1 month after season onset.
- The life cycle is approximately 1 year in QLD/NSW but potentially 1.5 years in Victoria (explaining the bimodal pattern).
- Negative correlation between successive seasons in Victoria may reflect herd immunity dynamics in host populations.

**Teo et al. (2024) — Weather-based prediction of tick paralysis:**

- Developed mixed-effects models to predict tick paralysis prevalence and season timing.
- Key predictors for annual prevalence: minimum temperature in summer, saturation deficit in summer, rainfall in summer, and ectoparasiticide uptake.
- Key predictors for season start timing: mean temperature in winter, saturation deficit in winter.
- Used SILO climate data (Australian Bureau of Meteorology's gridded dataset).
- Demonstrated that weather data from the preceding season can predict the following year's tick season.

**Clark et al. (2022) — Ensemble forecasting of tick paralysis:**

- Developed an iterative ensemble model combining ARIMA, GAM, GARCH, and Prophet forecasts for near-term (biweekly) tick paralysis incidence at a single veterinary clinic.
- Key environmental predictors: moisture-vegetation index anomalies, max temperature anomalies, NDVI, SOI (Southern Oscillation Index).
- Used seasonally-decomposed time series with MSIS-optimised ensemble weights.
- Ensemble outperformed benchmark ETS model, especially for peak season forecasting.
- Validated using rolling-window horizon error calculations and CRPS scoring.

### 1.2 What This Means for Our Model

For the current objective (Layer 1: modelling the current distribution of *I. holocyclus*), the literature strongly supports:

1. **A correlative SDM approach** using climate variables as the primary drivers, since this species is fundamentally climate-limited.
2. **An ensemble of multiple SDM algorithms** rather than a single model, to capture uncertainty and improve prediction robustness.
3. **Careful variable selection** informed by the species' known physiology, avoiding purely data-driven variable selection that may overfit.
4. **Inclusion of moisture/humidity variables** alongside temperature and rainfall—saturation deficit and vapour pressure deficit are repeatedly identified as important.
5. **Spatial thinning** of occurrence records to reduce sampling bias (the existing thinned dataset at 10 km is appropriate).
6. **Background point generation** using a target-group approach or geographic constraints to reduce sampling bias in pseudo-absences.

---

## 2. Objectives

### 2.1 Current Objective (This Analysis)

Model the current potential distribution (habitat suitability) of *Ixodes holocyclus* across eastern Australia using an ensemble of species distribution models.

### 2.2 Future Objectives (Not Addressed Here)

- Layer 2: Calibrate spatial risk using alpha-gal syndrome and tick-bite data.
- Layer 3: Dynamic temporal forecasting of tick activity.
- Layer 4: Integrated tick-bite risk index.

---

## 3. Data

### 3.1 Occurrence Data

**Source:** GBIF (already downloaded to `data/ihol_occurrences_gbif.csv`)

- 719 records with coordinates, dates, and state information.
- After spatial thinning at 10 km: ~236 records (already processed in `processed_data/thinned_occurrences/`).

**Cleaning steps required:**

1. Remove duplicate coordinates.
2. Remove records falling in the ocean or outside Australia.
3. Remove records with coordinate uncertainty >10 km (if available).
4. Assess temporal bias and optionally restrict to post-1970 records (to match climate normals).
5. Spatial thinning to reduce sampling bias (10 km thinning distance, already done).

### 3.2 Environmental Predictors

**Available data (already downloaded):**

| Source | Variables | Resolution | Path |
|--------|-----------|------------|------|
| WorldClim v2.1 | 19 bioclimatic variables (bio1–bio19) | 2.5 arcmin (~5 km) | `data/climate_data/climate/wc2.1_2.5m/` |
| WorldClim v2.1 | Monthly tmin, tmax, precipitation | 2.5 arcmin | `data/climate_data/climate/wc2.1_2.5m/` |
| ENVIREM | PET, aridity index, climatic moisture index, continentality, thermicity, etc. | 2.5 arcmin | `data/climate_data/climate/envirem/` |
| CHELSA v2.1 | Monthly vapour pressure deficit (VPD) | 30 arcsec (aggregate to 2.5 arcmin) | `data/climate_data/climate/chelsa/` |

**Variable selection strategy (informed by literature):**

Rather than using all 19 bioclimatic variables (which introduces multicollinearity and overfitting), we will use a hypothesis-driven approach:

**Core variables (physiologically justified):**

- **bio1** — Annual mean temperature (overall thermal suitability)
- **bio5** — Max temperature of warmest month (heat stress limit; ~32–33°C upper limit)
- **bio6** — Min temperature of coldest month (cold stress/developmental limit; ~8°C lower limit)
- **bio12** — Annual precipitation (moisture requirement for off-host survival)
- **bio15** — Precipitation seasonality (moisture reliability)
- **vpd_annual** — Mean annual vapour pressure deficit from CHELSA (desiccation risk; Teo et al. 2021, 2024 identified saturation deficit as key)
- **aridity_idx** — Thornthwaite aridity index from ENVIREM (integrated moisture availability)
- **cmi_idx** — Climatic moisture index from ENVIREM (moisture balance)

**Variable reduction procedure:**

1. Extract values at occurrence points and background points.
2. Calculate pairwise Pearson correlations.
3. Remove variables with |r| > 0.7, retaining the one with stronger ecological justification.
4. Assess variable importance in preliminary models and remove any with negligible contribution.
5. Check VIF (variance inflation factor); remove variables with VIF > 10.

Target: 5–7 uncorrelated, ecologically meaningful predictors in the final model.

### 3.3 Study Area

Eastern Australia, defined as longitude 140°E–155°E, latitude 10°S–40°S. This encompasses Queensland, New South Wales, ACT, Victoria, and the eastern portions of South Australia. The extent is slightly wider than the known range to capture the full environmental gradient and allow the models to identify range boundaries.

---

## 4. Modelling Approach

### 4.1 Overview

We will fit an **ensemble of five SDM algorithms**, each with different assumptions and strengths, and combine their predictions into a weighted ensemble. This is current best practice in SDM (Araujo & New 2007; Hao et al. 2019).

### 4.2 Background / Pseudo-absence Generation

Since we have presence-only data, we need to generate background points (pseudo-absences) for the discriminative models.

**Strategy:**

1. Generate 10,000 random background points within the study extent but constrained to land areas.
2. Optionally apply a **target-group bias correction**: use the geographic distribution of all Australian tick GBIF records as a sampling density surface, so that background points are drawn proportionally to survey effort. This corrects for the strong coastal and urban sampling bias in GBIF data.
3. Ensure background points are at least 5 km from any occurrence point (to avoid placing absences at presence locations).

### 4.3 Model Algorithms

| # | Algorithm | Type | R Package | Rationale |
|---|-----------|------|-----------|-----------|
| 1 | **MaxEnt** | Machine learning (maximum entropy) | `maxnet` | Gold standard for presence-background SDMs; handles complex response curves; regularisation prevents overfitting |
| 2 | **Random Forest** | Machine learning (ensemble of trees) | `randomForest` | Captures non-linear relationships and interactions; robust to outliers |
| 3 | **Boosted Regression Trees (BRT/GBM)** | Machine learning (gradient boosting) | `gbm` via `dismo` | Strong predictive performance; handles interactions; built-in variable importance |
| 4 | **Generalised Additive Model (GAM)** | Statistical (semi-parametric) | `mgcv` | Smooth, interpretable response curves; well-suited for ecological data |
| 5 | **Generalised Linear Model (GLM)** | Statistical (parametric) | base R | Simple, interpretable baseline; linear and quadratic terms allow testing of physiological thresholds |

### 4.4 Model Fitting

**For each algorithm:**

1. Fit model using thinned occurrence points (presence=1) and background points (absence=0).
2. Use the selected environmental predictors as covariates.
3. Tune hyperparameters where applicable:
   - MaxEnt: regularisation multiplier (test 0.5, 1, 2, 3, 5) and feature classes (L, LQ, LQH, LQHP).
   - Random Forest: number of trees (500–1000), mtry (sqrt of number of predictors).
   - BRT: learning rate (0.001–0.01), tree complexity (3–5), bag fraction (0.75).
   - GAM: basis dimension (k) and smoothing penalty.
   - GLM: linear + quadratic terms for temperature variables.

### 4.5 Model Evaluation

**Spatial cross-validation** using a 5-fold spatially blocked design (`blockCV` package):

- Blocks are defined geographically (e.g., ~200 km block size) so that training and test folds are spatially separated. This avoids inflated performance metrics from spatial autocorrelation.
- Calculate for each fold:
  - **AUC** (Area Under the ROC Curve): discrimination ability (>0.7 acceptable, >0.8 good).
  - **TSS** (True Skill Statistic): threshold-dependent accuracy (>0.4 acceptable, >0.6 good).
  - **Boyce Index** (continuous): evaluates predicted-to-expected ratio across suitability bins (positive = good).

**Additional checks:**

- Response curves: verify that each variable's modelled effect is ecologically plausible (e.g., suitability declines below 8°C and above 32°C).
- Variable importance: quantify each predictor's contribution.
- Compare model predictions against the known range map from Teo et al. (2021) as a qualitative sanity check.

### 4.6 Ensemble Construction

1. Weight individual models by their mean spatial cross-validation AUC.
2. Calculate the weighted mean predicted suitability across all five models for each grid cell.
3. Also calculate prediction uncertainty (standard deviation across models) for each cell.
4. Generate a binary presence/absence map using a threshold that maximises TSS.

---

## 5. Implementation Steps (Quarto Document Structure)

The analysis will be implemented as a single Quarto (.qmd) document with the following sections, each corresponding to a code chunk or series of chunks:

### Step 1: Setup and Configuration
- Load packages.
- Define paths and study area extent.
- Set random seed for reproducibility.

### Step 2: Occurrence Data Preparation
- Load GBIF occurrence data.
- Clean: remove duplicates, ocean points, restrict to study area.
- Optionally filter by date (post-1970).
- Spatial thinning at 10 km (use pre-thinned data or re-thin).
- Map occurrence points.

### Step 3: Environmental Data Preparation
- Load all WorldClim bioclimatic, ENVIREM, and CHELSA VPD rasters.
- Crop and mask to study area (eastern Australia land only).
- Resample ENVIREM and CHELSA to match WorldClim grid.
- Stack into a single multi-layer raster.
- Select candidate variables based on ecological justification (Section 3.2).

### Step 4: Variable Selection and Multicollinearity
- Extract environmental values at occurrence and background points.
- Pairwise correlation matrix; remove |r| > 0.7 keeping ecologically justified variable.
- VIF analysis; remove VIF > 10.
- Finalise predictor set.

### Step 5: Background Point Generation
- Generate random background points within study extent on land.
- Optionally: target-group bias file approach.
- Minimum 5 km buffer from occurrence points.
- Extract environmental values at background points.
- Combine occurrence + background into modelling dataset.

### Step 6: Spatial Cross-Validation Setup
- Create spatially-blocked folds using `blockCV`.
- Block size ~200 km (based on spatial autocorrelation range).
- 5 folds.

### Step 7: Individual Model Fitting and Evaluation
- For each of the 5 algorithms:
  - Fit model on training folds.
  - Predict on test fold.
  - Calculate AUC, TSS, Boyce Index.
  - Generate response curves.
  - Calculate variable importance.
- Summarise cross-validation results in a table.

### Step 8: Model Tuning (MaxEnt and BRT)
- Tune MaxEnt regularisation and feature classes using `ENMeval`.
- Tune BRT learning rate and tree complexity using `dismo::gbm.step`.
- Refit best configurations.

### Step 9: Ensemble Prediction
- Refit all 5 models on the full dataset.
- Predict habitat suitability across study area for each model.
- Compute AUC-weighted ensemble mean.
- Compute prediction uncertainty (SD across models).
- Apply TSS-optimised threshold for binary map.

### Step 10: Visualisation and Output
- Map: Continuous ensemble suitability surface.
- Map: Binary suitable/unsuitable with occurrence points overlay.
- Map: Prediction uncertainty.
- Map: Individual model predictions (panel).
- Response curves (panel for each variable).
- Variable importance bar chart.
- Cross-validation performance summary table.
- Save all GeoTIFF outputs.

### Step 11: Ecological Validation
- Compare predicted distribution against known range from Teo et al. (2021).
- Check if predicted suitability aligns with the 8°C cold limit, 32°C heat limit, and moisture constraints.
- Discuss any areas of overprediction or underprediction.

---

## 6. R Packages Required

```r
# Core spatial
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Data manipulation
library(tidyverse)
library(here)

# SDM fitting
library(maxnet)          # MaxEnt
library(randomForest)    # Random Forest
library(gbm)             # Boosted Regression Trees
library(dismo)           # BRT helper functions (gbm.step)
library(mgcv)            # GAMs

# SDM evaluation and tuning
library(blockCV)         # Spatial cross-validation
library(ENMeval)         # MaxEnt tuning
library(predicts)        # SDM evaluation metrics
library(ecospat)         # Boyce Index

# Variable selection
library(usdm)            # VIF analysis
library(corrplot)        # Correlation visualisation

# Spatial thinning
library(spThin)

# Visualisation
library(ggplot2)
library(patchwork)
library(viridis)
library(tmap)
```

---

## 7. Expected Outputs

1. **Ensemble habitat suitability map** — continuous 0–1 suitability across eastern Australia at 2.5 arcmin (~5 km) resolution.
2. **Binary distribution map** — threshold-based predicted presence/absence.
3. **Prediction uncertainty map** — standard deviation across model predictions.
4. **Individual model predictions** — for comparison and transparency.
5. **Cross-validation performance table** — AUC, TSS, Boyce Index for each algorithm.
6. **Variable importance rankings** — contribution of each predictor.
7. **Response curves** — fitted relationships between suitability and each predictor.
8. **GeoTIFF files** — all spatial outputs saved for future integration with Layers 2–4.

---

## 8. Quality Assurance

- All code will be reproducible from the Quarto document.
- Random seed set at the top of the document.
- Spatial cross-validation used (not random CV) to provide realistic performance estimates.
- Response curves checked against known physiological limits.
- Multiple algorithms with ensemble approach to reduce model-specific bias.
- Variables selected based on ecological knowledge, not purely statistical criteria.

---

## 9. Limitations and Assumptions

1. **Presence-only data:** GBIF records represent observed occurrences, not true absences. Background points are a proxy for available habitat.
2. **Sampling bias:** Records are biased toward urban and coastal areas. Spatial thinning and target-group background selection partially mitigate this.
3. **Climate normals:** WorldClim data represent 1970–2000 averages; the actual current distribution may reflect more recent climate conditions.
4. **Fundamental vs. realised niche:** SDMs model the realised niche (where the species is found given current biotic and dispersal constraints), not necessarily the full fundamental niche.
5. **Static model:** This Layer 1 model is a time-averaged snapshot. Temporal dynamics (seasonality, interannual variation) will be addressed in Layers 3 and 4.
6. **Resolution:** 2.5 arcmin (~5 km) may not capture fine-scale habitat variation (microclimate refugia, urban heat islands).

---

## 10. Timeline

| Step | Description | Estimated Effort |
|------|-------------|-----------------|
| 1–3 | Data preparation | Small |
| 4–5 | Variable selection and background generation | Small |
| 6 | Cross-validation setup | Small |
| 7–8 | Model fitting and tuning | Moderate |
| 9 | Ensemble construction | Small |
| 10–11 | Visualisation and validation | Moderate |

---

## 11. References

- Clark, N.J., et al. (2022). Near-term forecasting of companion animal tick paralysis incidence: An iterative ensemble model. *PLoS Computational Biology*.
- Teo, E.J.M., et al. (2021). Climatic requirements of the eastern paralysis tick, *Ixodes holocyclus*, with a consideration of its possible range shift. *Parasitology Research*.
- Teo, E.J.M., et al. (2023). Phenology of tick paralysis in dogs and cats in eastern Australia. *Parasites & Vectors*.
- Teo, E.J.M., et al. (2024). Weather predicts tick paralysis. *Parasites & Vectors*.
- Heath, A.C.G. (1974). The temperature and humidity preferences of *Haemaphysalis longicornis*, *Ixodes holocyclus* and *Rhipicephalus sanguineus*. *Journal of Medical Entomology*.
