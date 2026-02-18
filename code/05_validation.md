# 05 Validation
Alexander W. Gofton
2026-02-18

- [1. Setup and Configuration](#1-setup-and-configuration)
- [2. Load Outputs from Previous
  Scripts](#2-load-outputs-from-previous-scripts)
- [11. Ecological Validation](#11-ecological-validation)
  - [11.1 Cross-Validation Performance
    Summary](#111-cross-validation-performance-summary)
- [12. Session Info](#12-session-info)

# 1. Setup and Configuration

``` r
# --- Core spatial ---
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# --- Data manipulation ---
library(tidyverse)

# --- SDM fitting ---
#library(maxnet)
#library(randomForest)
#library(gbm)
#library(dismo)
#library(mgcv)

# --- SDM evaluation and tuning ---
#library(blockCV)
#library(ENMeval)
#library(ecospat)

# --- Variable selection ---
#library(corrplot)

# --- Spatial thinning ---
#library(spThin)

# --- Visualisation ---
library(patchwork)
library(viridis)

# --- Set seed for reproducibility ---
set.seed(7990)

# --- Paths ---
base_dir <- "/Users/gof005/Library/CloudStorage/OneDrive-CSIRO/OneDrive - Docs/01_Projects/alpha_gal/02_SAP_2025-6/Ihol_SDM"
data_dir <- file.path(base_dir, "data")
processed_dir <- file.path(base_dir, "processed_data")
output_dir <- file.path(base_dir, "outputs")
figures_dir <- file.path(base_dir, "figures")


# --- Study area extent (Eastern Australia) ---
study_extent <- ext(140, 155, -40, -10)
```

# 2. Load Outputs from Previous Scripts

``` r
# --- From 02_env_vars.qmd ---
env_all    <- rast(file.path(output_dir, "env_layers", "env_all.tif"))

# --- From 03_bg_points_model_fitting.qmd ---
bg_out_dir   <- file.path(output_dir, "bg_model_inputs")
model_data   <- readRDS(file.path(bg_out_dir, "model_data.rds"))
pred_cols    <- readRDS(file.path(bg_out_dir, "pred_cols.rds"))
algorithms   <- readRDS(file.path(bg_out_dir, "algorithms.rds"))
algo_names   <- readRDS(file.path(bg_out_dir, "algo_names.rds"))
cv_summary   <- read_csv(file.path(output_dir, "cv_performance_summary.csv"),
                         show_col_types = FALSE)

# --- From 04_full_model_fitting.qmd ---
sdm_out_dir    <- file.path(output_dir, "sdm_results")
ensemble_binary <- rast(file.path(sdm_out_dir, "ensemble_binary.tif"))
weights         <- readRDS(file.path(sdm_out_dir, "ensemble_weights.rds"))
tss_result      <- readRDS(file.path(sdm_out_dir, "tss_result.rds"))

cat("Loaded all inputs for validation.\n")
```

    Loaded all inputs for validation.

# 11. Ecological Validation

``` r
# --- Check against known physiological thresholds ---

# Extract environmental values at predicted suitable locations
suitable_cells <- as.data.frame(ensemble_binary, xy = TRUE, na.rm = TRUE)
suitable_cells <- suitable_cells[as.integer(suitable_cells[[3]]) == 1L, ]
suitable_vect <- vect(suitable_cells[, 1:2], geom = c("x", "y"), crs = "EPSG:4326")
suitable_env <- terra::extract(env_all, suitable_vect)

cat("=== Ecological Validation ===\n\n")

# Temperature checks (from Heath 1974, Teo et al. 2021)
if ("bio5" %in% names(suitable_env)) {
  bio5_range <- range(suitable_env$bio5, na.rm = TRUE)
  cat(sprintf("Max temp warmest month (bio5) in suitable area: %.1f - %.1f °C\n",
              bio5_range[1], bio5_range[2]))
  cat(sprintf("  Expected upper limit: ~32-33°C (Heat stress threshold)\n"))
  cat(sprintf("  Check: %s\n\n",
              ifelse(bio5_range[2] <= 35, "PASS", "WARNING - exceeds expected limit")))
}

if ("bio6" %in% names(suitable_env)) {
  bio6_range <- range(suitable_env$bio6, na.rm = TRUE)
  cat(sprintf("Min temp coldest month (bio6) in suitable area: %.1f - %.1f °C\n",
              bio6_range[1], bio6_range[2]))
  cat(sprintf("  Expected lower limit: ~8°C (Developmental threshold)\n"))
  cat(sprintf("  Check: %s\n\n",
              ifelse(bio6_range[1] >= -2, "PASS - within reasonable range",
                     "WARNING - well below expected limit")))
}

if ("bio12" %in% names(suitable_env)) {
  bio12_range <- range(suitable_env$bio12, na.rm = TRUE)
  cat(sprintf("Annual precipitation (bio12) in suitable area: %.0f - %.0f mm\n",
              bio12_range[1], bio12_range[2]))
  cat(sprintf("  Expected: >500 mm (moisture requirement for off-host survival)\n\n"))
}

# --- Geographic range check ---
cat("Geographic extent of predicted suitable habitat:\n")
cat(sprintf("  Latitude: %.1f to %.1f S\n",
            abs(max(suitable_cells$y)), abs(min(suitable_cells$y))))
cat(sprintf("  Longitude: %.1f to %.1f E\n",
            min(suitable_cells$x), max(suitable_cells$x)))
cat("\n")
cat("Expected range (Teo et al. 2021):\n")
cat("  Coastal fringe from north QLD (~15°S) to eastern VIC (~38°S)\n")
cat("  Predominantly east of the Great Dividing Range\n")
```

    === Ecological Validation ===

    Max temp warmest month (bio5) in suitable area: 21.8 - 32.4 °C
      Expected upper limit: ~32-33°C (Heat stress threshold)
      Check: PASS

    Min temp coldest month (bio6) in suitable area: 0.6 - 19.3 °C
      Expected lower limit: ~8°C (Developmental threshold)
      Check: PASS - within reasonable range

    Annual precipitation (bio12) in suitable area: 645 - 4496 mm
      Expected: >500 mm (moisture requirement for off-host survival)

    Geographic extent of predicted suitable habitat:
      Latitude: 12.6 to 39.0 S
      Longitude: 143.2 to 153.6 E

    Expected range (Teo et al. 2021):
      Coastal fringe from north QLD (~15°S) to eastern VIC (~38°S)
      Predominantly east of the Great Dividing Range

## 11.1 Cross-Validation Performance Summary

``` r
cat("\n")
cat("===================================================\n")
cat("  FINAL MODEL PERFORMANCE SUMMARY\n")
cat("===================================================\n\n")

print(cv_summary, row.names = FALSE)

cat("\n")
cat("Ensemble weights:\n")
for (i in seq_along(algo_names)) {
  cat(sprintf("  %15s: %.3f\n", algo_names[i], weights[i]))
}
cat(sprintf("\nTSS threshold for binary map: %.3f\n", tss_result$threshold))
cat(sprintf("Number of occurrence records used: %d\n", sum(model_data$presence == 1)))
cat(sprintf("Number of background points used: %d\n", sum(model_data$presence == 0)))
cat(sprintf("Number of environmental predictors: %d\n", length(pred_cols)))
cat(sprintf("Predictors: %s\n", paste(pred_cols, collapse = ", ")))
```


    ===================================================
      FINAL MODEL PERFORMANCE SUMMARY
    ===================================================

    # A tibble: 5 × 7
      Algorithm     AUC_mean AUC_sd TSS_mean TSS_sd Boyce_mean Boyce_sd
      <chr>            <dbl>  <dbl>    <dbl>  <dbl>      <dbl>    <dbl>
    1 MaxEnt           0.968  0.014    0.881  0.059      0.711    0.174
    2 Random Forest    0.964  0.014    0.857  0.054      0.559    0.168
    3 BRT              0.962  0.014    0.838  0.058      0.519    0.152
    4 GAM              0.967  0.012    0.873  0.062      0.558    0.173
    5 GLM              0.966  0.015    0.871  0.068      0.662    0.145

    Ensemble weights:
               MaxEnt: 0.201
        Random Forest: 0.200
                  BRT: 0.199
                  GAM: 0.200
                  GLM: 0.200

    TSS threshold for binary map: 0.334
    Number of occurrence records used: 224
    Number of background points used: 5489
    Number of environmental predictors: 6
    Predictors: bio5, bio6, bio12, aridity_idx, cmi_idx, vpd_annual

# 12. Session Info

``` r
sessionInfo()
```

    R version 4.4.1 (2024-06-14)
    Platform: aarch64-apple-darwin20
    Running under: macOS 26.2

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

    locale:
    [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

    time zone: Australia/Brisbane
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] viridis_0.6.5           viridisLite_0.4.2       patchwork_1.3.2        
     [4] lubridate_1.9.4         forcats_1.0.0           stringr_1.5.1          
     [7] dplyr_1.1.4             purrr_1.1.0             readr_2.1.5            
    [10] tidyr_1.3.1             tibble_3.3.0            ggplot2_4.0.0          
    [13] tidyverse_2.0.0         rnaturalearthdata_1.0.0 rnaturalearth_1.1.0    
    [16] sf_1.0-21               terra_1.8-93           

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6       xfun_0.53          tzdb_0.5.0         vctrs_0.6.5       
     [5] tools_4.4.1        generics_0.1.4     parallel_4.4.1     proxy_0.4-27      
     [9] pkgconfig_2.0.3    KernSmooth_2.23-26 RColorBrewer_1.1-3 S7_0.2.0          
    [13] lifecycle_1.0.4    compiler_4.4.1     farver_2.1.2       codetools_0.2-20  
    [17] htmltools_0.5.8.1  class_7.3-23       yaml_2.3.10        pillar_1.11.0     
    [21] crayon_1.5.3       classInt_0.4-11    tidyselect_1.2.1   digest_0.6.37     
    [25] stringi_1.8.7      fastmap_1.2.0      grid_4.4.1         archive_1.1.12.1  
    [29] cli_3.6.5          magrittr_2.0.3     utf8_1.2.6         e1071_1.7-16      
    [33] withr_3.0.2        scales_1.4.0       bit64_4.6.0-1      timechange_0.3.0  
    [37] rmarkdown_2.29     bit_4.6.0          gridExtra_2.3      hms_1.1.3         
    [41] evaluate_1.0.5     knitr_1.50         rlang_1.1.6        Rcpp_1.1.0        
    [45] glue_1.8.0         DBI_1.2.3          vroom_1.6.5        jsonlite_2.0.0    
    [49] R6_2.6.1           units_0.8-7       
