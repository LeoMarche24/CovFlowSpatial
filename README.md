# Covariance modeling for flow-driven spatial domains via directed linear networks

## Overview

This repository contains the R code and data to reproduce the analyses in:

> *Covariance modeling for flow-driven spatial domains via directed linear networks*

We introduce a convolution-based covariance framework for geostatistical domains constrained by physical barriers and directed flows. The marine domain (northern Tyrrhenian Sea) is discretized into a directed linear network encoding ocean current directionality via a Markovian transition matrix. A penalized estimator is proposed for the resulting weighted covariance structure. The model is embedded in a Monte Carlo simulation framework applied to RCP 4.5 SST projections to identify thermal hot spots of ecological risk.

The repository implements the analysis on two domains:

- **Sardinia domain** (`*Sardinia*`): the full-scale domain used for all results reported in the paper. The hot-spot identification step (Step 6) is computationally intensive and requires an HPC cluster.
- **Tyrrhenian (Temp) domain** (`*Temp*`): a reduced toy domain designed to allow the complete pipeline to be run on a local machine, suitable for exploration and testing.

---

## Repository Structure

```
ConvProcHotSpot/
├── Data/
│   ├── Pre_Processing/                      # Raw data preprocessing scripts
│   ├── Data_Sardinia.RData                  # Processed SST residuals – Sardinia domain
│   ├── Data_Temp.RData                      # Processed SST residuals – Tyrrhenian domain
│   ├── Covariance_Data_45_Sardinia.RData    # Covariance estimation inputs – Sardinia, RCP 4.5
│   ├── Covariance_Data_45_Temp.RData        # Covariance estimation inputs – Tyrrhenian, RCP 4.5
│   ├── Projections_Sardinia.RData           # RCP 4.5 SST projections – Sardinia domain
│   ├── Projections_Temp.RData               # RCP 4.5 SST projections – Tyrrhenian domain
│   ├── data_projection_Sardinia_2050.RData  # Simulated SST fields – Sardinia, 2050
│   ├── data_projection_Temp_2050.RData      # Simulated SST fields – Tyrrhenian, 2050
│   ├── bbox_Sardinia.RData                  # Bounding box – Sardinia domain
│   └── bbox_Temp.RData                      # Bounding box – Tyrrhenian domain
├── HPC/                                     # Cluster submission scripts (see HPC section)
├── functions_network.R                      # Core functions: directed network construction,
│                                            #   transition matrix, path enumeration (Sections 2.2–2.3)
├── functions_observations.R                 # Core functions: data preprocessing, spatial alignment,
│                                            #   residual computation (Section 5.1)
├── exploratory.R                            # Exploratory spatial analysis and network visualization
│                                            #   (Figures 1, 3)
├── covariance_estimation.R                  # Penalized covariance estimator (Section 4,
│                                            #   Proposition 5); reproduces Figure 6
├── simulation_study.R                       # Simulation study validating the estimation procedure
│                                            #   (Supplementary Material)
├── hot_spot_identification.R                 # Monte Carlo SST simulation + excursion sets
│                                            #   (Section 5.2); reproduces Figures 7, 8
├── plot_hotSpot.R                           # Final hot-spot visualization (Figure 8)
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
  "igraph",      # directed network construction and path enumeration
  "sf",          # spatial data handling
  "terra",       # raster operations
  "ggplot2",     # visualization
  "tidyverse"    # data manipulation
))

---

## How to Reproduce — Step-by-Step

Run scripts in the following order. All scripts assume the working directory is the repository root.

| Step | Script | What it does | Paper reference |
|------|--------|--------------|-----------------|
| 1 | `functions_network.R` | **Source first.** Defines all network construction and covariance functions | Sections 2.2, 2.3, 3 |
| 2 | `functions_observations.R` | **Source first.** Defines data loading, alignment and residual utilities | Section 5.1 |
| 3 | `Exploratory.R` | Exploratory analysis; builds and visualizes the linear network | Figures 1, 3 |
| 4 | `Covariance Estimation.R` | Fits the penalized covariance estimator over 2006–2022 residuals | Section 4, Figure 6 |
| 5 | `Simulation Study.R` | Runs the simulation study validating the estimator | Supplementary Material |
| 6 | `HotSpot Identification.R` | Generates M = 500 Monte Carlo SST realizations; computes excursion sets | Section 5.2, Figures 7, 8 |
| 7 | `Plot HotSpot.R` | Produces final hot-spot maps | Figure 8 |

> ⚠️ **Steps 5 and 6 on the Sardinia domain are computationally intensive** and require an HPC cluster. The Tyrrhenian (Temp) domain can be used to run the full pipeline locally for testing purposes. The active domain can be switched from `"Sardinia"` to `"Temp"` at the top of each script.

---

## HPC

Steps 5–6 on the Sardinia domain were executed on an HPC cluster using [Apptainer](https://apptainer.org/) (formerly Singularity) for containerized and fully reproducible execution.

### Container

The runtime environment is packaged as an Apptainer image built from a Docker base image available at:  
[https://app.docker.com/accounts/leonardo2405](https://app.docker.com/accounts/leonardo2405)

To pull and convert the Docker image to an Apptainer `.sif` file:

```bash
apptainer pull myimage.sif docker://leonardo2405/r-geo-phd:v3
```

### Running the Analysis

```bash
apptainer exec myimage.sif Rscript "Name_Script"
```

Submission scripts for the PBS scheduler are provided in `HPC/`. The Sardinia hot-spot identification step requires approximately **30 GB of RAM** and completes in **under 30 minutes** on the cluster.

---

## License

This repository is licensed under the [MIT License](LICENSE).

---
