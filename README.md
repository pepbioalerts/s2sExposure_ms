# s2sExposure

Scripts to run an exposure model from subseasonal-to-seasonal (S2S) forecast data for an example species, and to compute a meta-ranking for regional prioritization.

This repository provides:
- An exposure modeling workflow using GEOS-S2S-V2 forecasts for a focal species.
- A meta-ranking pipeline to prioritize regions based on multi-criteria scoring.
- Modular R functions and an example end-to-end workflow that you can adapt to other species or regions.

R is the only language used in this repository.

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Repository Structure](#repository-structure)
- [Data Access](#data-access)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
  - [Workflow: Exposure Forecasts](#workflow-exposure-forecasts)
  - [Workflow: Meta-Ranking](#workflow-meta-ranking)
- [Reproducibility](#reproducibility)
- [Performance Notes](#performance-notes)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Acknowledgements](#acknowledgements)

---

## Overview

Subseasonal-to-seasonal (S2S) forecasts bridge the gap between weather forecasts and climate projections. This repository demonstrates how to:
1. Ingest GEOS-S2S-V2 forecast data for target variables relevant to species exposure (e.g., temperature in our specific case).
2. Overlay forecasts with species expert range maps (e.g., IUCN) to estimate exposure metrics.
3. Aggregate and rank regions based on exposure and additional indicators (meta-ranking).

The workflows here can be adapted to different species, variables, spatial domains, and temporal horizons.

---

## Key Features

- End-to-end example workflow (`workflowExample.R`) for exposure computation.
- Modular functions (`funcs_*.R`) for:
  - Data ingestion and preprocessing
  - Spatial operations and masking with species ranges
  - Exposure metric calculation and summarization
  - Ranking and meta-ranking assembly
  - Visualization and export
- Clear separation of scripts and outputs to facilitate reproducibility.

---

## Repository Structure

- `code/`
  - `funcs_*.R`: Modular functions used throughout the workflows.
    - Common modules you’ll find (names may vary slightly by version):
      - `funcs_data_ingest.R`: Helpers to read GEOS-S2S-V2 data and harmonize dimensions.
      - `funcs_spatial.R`: Spatial utilities for reprojection, masking, and zonal stats with species ranges.
      - `funcs_exposure.R`: Exposure metric definitions and aggregations.
      - `funcs_ranking.R`: Ranking and meta-ranking computation utilities.
      - `funcs_plotting.R`: Basic plotting and export functions.
  - `workflowExample.R`: Example script showing how to run the exposure workflow end-to-end on a focal species.
  - `ranking/*.R`: Scripts for building meta-ranking across regions (weights, normalization, compositing).
- `data/` (optional; not versioned here)
  - Paths you configure locally for input and intermediate files.
- `outputs/` (optional; created by workflow)
  - Saved results, figures, and tables from runs.

Note: Folder names may differ depending on your setup; the key components are the modular functions and example workflows inside `code/`.

---

## Data Access

You will need two primary data sources:

1. Species expert range maps (IUCN):
   - Available from IUCN: https://www.iucnredlist.org/resources/spatial-data-download
   - Ensure you have the required permissions and comply with IUCN data usage terms.

2. GEOS-S2S-V2 forecasts (NCCS Discover):
   - Available on the Discover server of NASA NCCS.
   - Consult your local documentation or NCCS guidance for access credentials and data paths.
   - Typical variables: temperature, precipitation, humidity, etc., depending on your exposure metric.

If you do not have access to Discover, you can adapt the ingestion functions to other S2S datasets (ECMWF, CFSv2, etc.), provided you harmonize variable names, grids, and time coordinates.


Species expert range maps from IUCN are available at available at https://www.iucnredlist.org/resources/spatial-data-download. The GEOS-S2S-V2 data are available on the Discover server of NCCS at https://www.nccs.nasa.gov/systems/ data-portal, and GEOS-S2S-V2 forecasts output data are presently available at https://gmao.gsfc.nasa.gov/gmaoftp/gmaofcst/. Historical climate reanalysis data from ERA5 are available at https://cds.climate.copernicus.eu/ ERA5. Administrative units available at gadm.org v. 4.1.  Climate exposure results are synthetized by large groups (e.g. mammals, birds, amphibians, reptiles) to respect species expert ranges licensing and protect vulnerable species. These and resulting data are available at a dedicated GitHub repository https://github.com/pepbioalerts/s2sExposure_ms

---

## Prerequisites

- R (≥ 4.1 recommended)
- Suggested R packages:
  - `terra` or `sf` (spatial operations)
  - `stars` or `raster` (gridded data handling)
  - `ncdf4` (NetCDF IO for S2S data)
  - `tidyverse` (data manipulation and plotting)
  - `lubridate` (time handling)
  - `exactextractr` (zonal stats; optional)
  - `data.table` (fast aggregation; optional)
  - `ggplot2` (visualization)
  - `yaml` or `config` (configuration; optional)

Install with:
```r
install.packages(c("terra", "sf", "stars", "ncdf4", "tidyverse", "lubridate", "ggplot2"))
# Optional:
install.packages(c("exactextractr", "data.table", "yaml", "config"))
```

You may need system libraries for `sf/terra` (GEOS, GDAL, PROJ). On Linux/macOS, use your package manager; on Windows, install via RTools or pre-compiled binaries.



## Usage

### Workflow: Exposure Forecasts

1. Prepare inputs:
   - Download IUCN range maps for your target species.
   - Ensure GEOS-S2S-V2 data is accessible and the path is correctly set.
   - Define your exposure variable and threshold.

2. Run the example workflow:
   - Open `code/workflowExample.R`.
   - Update paths and parameters to match your environment.
   - Source required function files:
     ```r
     source("code/funcs_data_ingest.R")
     source("code/funcs_spatial.R")
     source("code/funcs_exposure.R")
     source("code/funcs_plotting.R") # optional
     ```
   - Execute the workflow to:
     - Read S2S forecasts (NetCDF).
     - Reproject/mask with species range polygons.
     - Calculate exposure metrics (e.g., fraction of days exceeding threshold within range).
     - Aggregate by time windows and export tables/plots to `outputs/`.

3. Outputs:
   - Tables (CSV) with exposure metrics by lead time and region.
   - Figures (PNG/PDF) visualizing exposure across time and space.
   - Intermediate RDS files for reproducibility.

### Workflow: Meta-Ranking

1. Prepare regional shapefile and indicators:
   - Regional boundaries for prioritization (e.g., administrative regions).
   - Indicators (exposure, variability, trends) either computed from the exposure workflow or other sources.

2. Run ranking scripts:
   - Source ranking utilities:
     ```r
     source("code/funcs_ranking.R")
     ```
   - Configure weights and normalization (see `config.yaml` example).
   - Compute composite scores and rankings per region.
   - Export results to `outputs/ranking/`.

3. Outputs:
   - Ranked list of regions with composite scores.
   - Diagnostic plots and sensitivity tests (if implemented).




## Citation

If you use this repository or derived workflows in a publication, please cite:
- IUCN Red List spatial data: [IUCN Red List Spatial Data](https://www.iucnredlist.org/resources/spatial-data-download)
- GEOS-S2S-V2: Cite the appropriate NASA NCCS/GEOS references for S2S data access and model details (consult NCCS documentation).

Additionally, cite this repository as:
- pepbioalerts (Year). s2sExposure_ms: Scripts to run an exposure model from S2S data. GitHub repository. https://github.com/pepbioalerts/s2sExposure_ms

---

## License

Please refer to the repository’s license file (if present). If no license is specified, usage may be restricted; contact the maintainers for permissions.

---

## Acknowledgements

- IUCN for providing species expert range maps.
- NASA NCCS for Discover access and GEOS-S2S-V2 data resources.
- Contributors to the R spatial ecosystem (`sf`, `terra`, `stars`, etc.).
