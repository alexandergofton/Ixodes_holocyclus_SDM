# Project Plan: Forecasting *Ixodes holocyclus* Tick-Bite Risk in Eastern Australia

## 1. Project Goal

The primary goal of this project is to develop a scientifically robust, spatially explicit tick-bite risk index for the eastern paralysis tick, *Ixodes holocyclus*, across Eastern Australia. This index will be resolved to a 5x5km grid and will represent the probability of a human or animal being bitten by a tick. The model will integrate climate data, environmental features, and novel data on tick-borne disease distribution to produce dynamic risk maps that are useful for the public, veterinarians, and public health officials.

## 2. Core Scientific Questions

This project will address the following key questions:
- Can we define and model a "tick-bite risk index" that is more granular and spatially explicit than simple case counts?
- How can we effectively integrate quantitative data on alpha-gal syndrome to represent and validate spatial variations in tick-human contact rates?
- What are the key environmental, geographic, and climatic drivers of this tick-bite risk index at a fine spatial scale (5x5km)?
- How does this risk vary seasonally and geographically across the diverse landscape of Eastern Australia?

## 3. Methodology: A Multi-layered Modeling Approach

The model will be constructed as a hierarchical framework in R. Each layer will model a different component of tick risk, building upon the outputs of the preceding layers.

### Layer 1: Foundational Habitat Suitability (The "Where")

*   **Objective:** To map the fundamental niche of *I. holocyclus*, defining the geographical areas where it can establish and sustain populations.
*   **Method:**
    1.  Replicate and refine the Species Distribution Models (SDMs) from Teo et al. (2021). The "climatic-range method" outlined in the paper will be prioritized due to its parsimony and high performance.
    2.  **Input Data:**
        *   Climate Data: Mean monthly temperature, rainfall, and relative humidity layers from a reputable source like WorldClim or the `climatedata` R package.
        *   Occurrence Data: *I. holocyclus* presence records from `data/ihol_occurrences_gbif.csv` will be used for model training and validation.
    3.  **R Packages:** `terra` for spatial data handling, `dismo` and `predicts` for model fitting and projection.
*   **Output:** A static 5x5km raster map of Eastern Australia with a habitat suitability score (0 to 1). This map will serve as the baseline constraint for all subsequent risk modeling.

### Layer 2: Spatial Risk Calibration (The "Human-Tick Interaction Hotspots")

*   **Objective:** To refine the baseline habitat suitability with real-world data on tick-human interaction, identifying where bites are most likely to occur.
*   **Method:**
    1.  Develop a spatial regression model where the density of alpha-gal syndrome cases is the response variable. This will allow us to "weight" the habitat suitability by observed public health outcomes.
    2.  **Input Data:**
        *   **Response Variable:** The provided quantitative data on the geographic distribution of alpha-gal syndrome cases will be processed into a case density surface on the 5x5km grid.
        *   **Predictor Variables:**
            *   Layer 1 Habitat Suitability Score.
            *   Geographic Features: Elevation, slope, and aspect from a Digital Elevation Model (e.g., from the `elevatr` R package).
            *   Environmental Features: Land cover data (e.g., from Copernicus or MODIS) to extract proportions of forest, shrubland, and urban interfaces.
            *   Socio-economic Data: Population density as a potential covariate to account for observation/reporting bias.
    3.  **Model:** A spatial generalized linear model (GLM) or a machine learning model like Random Forest will be fitted to explain the variation in alpha-gal cases.
    4.  **R Packages:** `sf` and `terra` for data processing, `mgcv` for GAMs, or `randomForest` for the RF model.
*   **Output:** A "Human-Tick Interaction" risk surface. This layer adjusts the baseline habitat suitability based on where people are most frequently bitten.

### Layer 3: Dynamic Tick Activity Forecast (The "When" and "How Much")

*   **Objective:** To forecast the seasonal and short-term (weekly) changes in tick activity across the landscape.
*   **Method:**
    1.  Build an ensemble forecasting model, inspired by the work of Clark et al. (2022), for each grid cell.
    2.  **Ensemble Components:** The ensemble will combine forecasts from several models to improve accuracy and robustness.
        *   Time Series Models: `ARIMA` and `GARCH` models (from the `forecast` and `rugarch` packages) to capture temporal autocorrelation.
        *   Generative Models: `GAM` (`mgcv`) and `Prophet` to model non-linear relationships with environmental drivers.
    3.  **Input Data (as time-series for each grid cell):**
        *   **Response Variable Proxy:** A "tick activity index" will be generated based on the phenology models from Teo et al. (2023), incorporating the described latitudinal shifts and the bimodal pattern for Victoria.
        *   **Predictor Variables:** Weather data (temperature, rainfall, saturation deficit), NDVI, and the Southern Oscillation Index (SOI).
*   **Output:** Weekly forecasts of the "tick activity index" for each 5x5km grid cell.

### Layer 4: Final Tick-Bite Risk Index Synthesis

*   **Objective:** To combine the outputs of the three layers into a single, final, and easily interpretable risk index.
*   **Method:**
    1.  The final risk index for each grid cell `i` at time `t` will be calculated as a product of the preceding layers:
        `RiskIndex(i, t) = HabitatSuitability(i) * HumanTickInteraction(i) * TickActivityForecast(i, t)`
    2.  This multiplicative approach ensures that risk is only non-zero in areas where ticks can live, are known to interact with humans, and are seasonally active.
    3.  The final index will be scaled (e.g., 0-100) and classified into intuitive risk categories (e.g., Low, Moderate, High, Very High) for public communication.

## 4. Model Validation

-   **Layer 1 (Habitat Suitability):** Standard spatial cross-validation techniques will be used, with evaluation metrics such as AUC (Area Under the Curve) and TSS (True Skill Statistic).
-   **Layer 2 (Spatial Risk):** A portion of the alpha-gal syndrome data will be held out as a test set to validate the model's ability to predict hotspots in new areas.
-   **Layer 3 (Activity Forecast):** The temporal patterns of the forecasted tick activity index will be validated against regional tick paralysis case data from Teo et al. (2023).

## 5. Project Implementation in R

-   **Project Structure:** A clear project directory will be established (`data/`, `R/`, `outputs/`).
-   **Data Management:** R scripts will be written to automate the downloading and preprocessing of all data, aligning it to the master 5x5km grid for Eastern Australia.
-   **Coding:** The entire workflow will be implemented in R. The code will be version-controlled using Git.

## 6. Deliverables

1.  **A commented set of R scripts** that allows for the full reproduction of the analysis.
2.  **The final tick-bite risk model** saved as an R object or set of model files.
3.  **A series of GeoTIFF raster files** representing the weekly tick-bite risk index.
4.  **This markdown document**, updated with the project's results, figures, and conclusions.
