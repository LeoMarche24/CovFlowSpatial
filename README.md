# Covariance modeling for flow-driven spatial domains via directed linear networks

## Overview

This repository contains the R code and data to reproduce the analyses in:

> *Covariance modeling for flow-driven spatial domains via directed linear networks*

We introduce a convolution-based covariance framework for geostatistical domains constrained by physical barriers and directed flows. The marine domain (northern Tyrrhenian Sea) is discretized into a directed linear network encoding ocean current directionality via a Markovian transition matrix. A penalized estimator is proposed for the resulting weighted covariance structure. The model is embedded in a Monte Carlo simulation framework applied to RCP 4.5 SST projections to identify thermal hot spots of ecological risk.

The repository implements the analysis on two domains:

- **Sardinia domain** (`*Sardinia*`): the full-scale domain used for all results reported in the paper. The main estimation and simulation steps are computationally intensive and are intended to be run on an HPC cluster.
- **Tyrrhenian (Temp) domain** (`*Temp*`): a reduced toy domain suitable for local exploration and testing. Domain-specific settings are defined at the top of each analysis script.

---

## Repository Structure

```
CovFlowSpatial/
├── Data/
│   ├── Pre_Processing/                      # Raw data preprocessing scripts
│   ├── Data_Sardinia.RData                  # Processed SST observations – Sardinia domain
│   ├── Data_Temp.RData                      # Processed SST observations – Tyrrhenian domain
│   ├── Covariance_Data_45_*.RData           # Final fitted covariance parameters – RCP 4.5
│   ├── Covariance_Estimation_45_*.RData     # Yearly covariance-estimation results
│   ├── Cross_Validation_*.RData             # Rolling-block cross-validation results
│   ├── Projections_Sardinia.RData           # RCP 4.5 SST projections – Sardinia domain
│   ├── Projections_Temp.RData               # RCP 4.5 SST projections – Tyrrhenian domain
│   ├── data_projection_Sardinia_2050.RData  # Simulated SST fields – Sardinia, 2050
│   ├── data_projection_Temp_2050.RData      # Simulated SST fields – Tyrrhenian, 2050
│   ├── bbox_Sardinia.RData                  # Bounding box – Sardinia domain
│   └── bbox_Temp.RData                      # Bounding box – Tyrrhenian domain
├── jobs/                                    # PBS cluster submission scripts
├── functions_network.R                      # Core functions: directed network construction,
│                                            #   transition matrix, path enumeration (Sections 2.2–2.3)
├── functions_observations.R                 # Covariance estimation, fitting and prediction functions
├── Exploratory.R                            # Exploratory spatial analysis and network visualization
│                                            #   (Figures 1, 3)
├── covariance_estimation.R                  # Penalized covariance estimator (Section 4,
│                                            #   Proposition 5); reproduces Figure 5
├── covariance_function_plotting.R           # Covariance diagnostics and final parameter estimates
├── simulation_study.R                       # Simulation study validating the estimation procedure
│                                            #   (Section 5 and Supplementary Material)
├── sensibility_estimator.R                  # Sensitivity analysis for the covariance estimator
├── cross_validation_rolling_blocks.R        # Rolling spatial-block cross-validation
├── plot_results.R                           # Cross-validation plots and numerical summaries
├── hot_spot_identification.R                # Monte Carlo SST simulation + excursion sets
├── plot_hot_spot.R                          # Final hot-spot visualization (Figure 11)
├── ConvProcessLinearNetwork.Rproj           # RStudio project file
└── README.md
```

---

## Data Access

Two data sources are required. Both are freely accessible upon registration:

### 1. RCP Projections — Copernicus Climate Change Service (C3S)

- **Source:** Copernicus Climate Change Service (2020). *Marine biogeochemistry data for the Northwest European Shelf and Mediterranean Sea from 2006 up to 2100 derived from climate projections.*
- **DOI:** [10.24381/cds.dcc9295c](https://doi.org/10.24381/cds.dcc9295c)
- **Variables needed:** Monthly mean sea surface temperature (SST) and sea current velocity components (u, v) for August, RCP 4.5 scenario, Mediterranean domain, 2006–2099.
- **Registration:** Free account required at [cds.climate.copernicus.eu](https://cds.climate.copernicus.eu)

### 2. Satellite Observations — Copernicus Marine Service (CMEMS)

- **Source:** E.U. Copernicus Marine Service Information (2020). *Multi Observation Global Ocean 3D Temperature Salinity Height Geostrophic Current and MLD.*
- **DOI:** [10.48670/moi-00052](https://doi.org/10.48670/moi-00052)
- **Variables needed:** Monthly mean SST and current velocity for the Mediterranean Sea, August, 2006–2022.
- **Registration:** Free account required at [marine.copernicus.eu](https://marine.copernicus.eu)

Once downloaded, place the raw files in `Data/Pre_Processing/` and run the scripts therein to generate the processed `.RData` files expected by the main scripts.

> **Note:** The processed `.RData` files in `Data/` are included in this repository and are sufficient to reproduce all paper results without re-downloading the raw data.

---

## Requirements

- **R** ≥ 4.2.0
- Key packages:

```r
install.packages(c(
  "sp",
  "gstat",
  "ncdf4",
  "ggplot2",
  "spatstat",
  "igraph",
  "sf",
  "progress",
  "ggmap",
  "ggforce",
  "tidyverse",
  "plotly",
  "reshape2",
  "geosphere"
))
```

---

## How to Reproduce — Step-by-Step

Run scripts in the following order. All scripts assume the working directory is the repository root. The helper files `functions_network.R` and `functions_observations.R` are sourced automatically by the analysis scripts.

| Step | Script | What it does | Paper reference |
|------|--------|--------------|-----------------|
| 1 | `Exploratory.R` | Exploratory analysis; builds and visualizes the linear network | Figures 1, 2 and 3 |
| 2 | `covariance_estimation.R` | Fits the penalized covariance estimator over 2006–2022 residuals | Section 4 and Section 6.1 |
| 3 | `covariance_function_plotting.R` | Produces covariance diagnostics and saves the final fitted parameters | Section 6.1, Figure 10 |
| 4 | `simulation_study.R` | Runs the simulation study validating the estimator | Section 5 and Supplementary Material |
| 5 | `sensibility_estimator.R` | Evaluates sensitivity to the number of covariance bins | Supplementary Material |
| 6 | `cross_validation_rolling_blocks.R` | Performs rolling spatial-block cross-validation | Supplementary Material |
| 7 | `plot_results.R` | Produces cross-validation plots and numerical summaries | Supplementary Material |
| 8 | `hot_spot_identification.R` | Generates M = 500 Monte Carlo SST realizations; computes excursion sets | Section 6.2 |
| 9 | `plot_hot_spot.R` | Produces final hot-spot maps | Figure 11 |

> ⚠️ **The estimation and simulation scripts on the Sardinia domain are computationally intensive.** PBS submission files are provided in `jobs/`. The active domain can be switched from `"Sardinia"` to `"Temp"` at the top of each analysis script, with settings adjusted to the size of the selected domain when needed.

---

## HPC

The computationally intensive analyses on the Sardinia domain were executed on an HPC cluster using [Apptainer](https://apptainer.org/) (formerly Singularity) for containerized execution.

### Container

The runtime environment is packaged as an Apptainer image built from a Docker base image available at:  
[https://app.docker.com/accounts/leonardo2405](https://app.docker.com/accounts/leonardo2405)

To pull and convert the Docker image to an Apptainer `.sif` file:

```bash
apptainer pull r-geo-phd_v3.0.sif docker://leonardo2405/r-geo-phd:v3
```

### Running the Analysis

```bash
apptainer exec r-geo-phd_v3.0.sif Rscript simulation_study.R
```

Submission scripts for the PBS scheduler are provided in `jobs/`. They assume that jobs are submitted from that directory and move to the repository root before running the corresponding R script:

```bash
cd jobs
qsub simulation_study.job
```

Before submission, update `container_img` and the Apptainer executable path in the `.job` files to match the local cluster configuration. Requested CPUs, memory and wall time are specified in each submission file.

---

## License

This repository is licensed under the [MIT License](LICENSE.txt).

---
