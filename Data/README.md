# Data provenance and attribution

The `.RData` files cover the `Sardinia` and `Temp` domains. The source citations and reuse terms below apply to the Copernicus material in these files.

## Climate projections

Kay, S., Clark, J., Hall, A., Marsh, J., and Fernandes, J. (2020). *Marine biogeochemistry data for the Northwest European Shelf and Mediterranean Sea from 2006 up to 2100 derived from climate projections*. Copernicus Climate Change Service (C3S), Climate Data Store. [DOI 10.24381/cds.dcc9295c](https://doi.org/10.24381/cds.dcc9295c) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

`Projections_*.RData` contains modified C3S information. Preprocessing selects August surface temperature and current components from the POLCOMS–ERSEM projections for 2006–2022, 2050, and 2099 under RCP 4.5 and RCP 8.5. It crops the data to the two study domains, adds coordinates in EPSG:3035, and saves the tables as R objects. [Extraction notebook](Pre_Processing/DataClean.ipynb) · [Data treatment script](Pre_Processing/Data_Treatment.R).

## Ocean observations

*Multi Observation Global Ocean 3D Temperature Salinity Height Geostrophic Current and MLD* (`MULTIOBS_GLO_PHY_TSUV_3D_MYNRT_015_012`). E.U. Copernicus Marine Service Information, Marine Data Store. [DOI 10.48670/moi-00052](https://doi.org/10.48670/moi-00052) · [Product page](https://data.marine.copernicus.eu/product/MULTIOBS_GLO_PHY_TSUV_3D_MYNRT_015_012/description).

`Data_*.RData` contains processed data from the European Union Copernicus Marine Service, © Mercator Ocean. Preprocessing selects August sea surface temperature and geostrophic current components for 2006–2022, crops the data to the study domains, and adds coordinates in EPSG:3035. [Observation extraction script](Pre_Processing/Treatment_RealData.R) · [Data treatment script](Pre_Processing/Data_Treatment.R).

Reuse of the Copernicus Marine material is subject to the service licence and its attribution requirements. [Copernicus Marine Service licence](https://marine.copernicus.eu/user-corner/service-commitments-and-licence) · [Attribution guidance](https://help.marine.copernicus.eu/en/articles/4444611-citing-copernicus-marine-products-and-services).

## Derived results and domain settings

| Files | Contents and provenance |
| --- | --- |
| `Covariance_Estimation_45_*.RData`, `Covariance_Data_45_*.RData` | Yearly covariance estimates and final fitted parameters from the Copernicus Marine observations and C3S RCP 4.5 projections. |
| `Cross_Validation_Rolling_Blocks_SST_*.RData` | Cross-validation results using the Copernicus Marine observations and C3S RCP 4.5 projections. |
| `data_projection_*_2050.RData` | The 2050 hot-spot analysis output, combining C3S RCP 4.5 projections with model parameters derived from both sources. |
| `range_heatmap_*.RData` | A range parameter computed from the network built using the 2022 Copernicus Marine data. |
| `bbox_*.RData` | Bounding boxes defined in the preprocessing script for the two study domains. |

The derived results that use marine observations were generated using E.U. Copernicus Marine Service Information, [DOI 10.48670/moi-00052](https://doi.org/10.48670/moi-00052). The results based on climate projections also use C3S information, [DOI 10.24381/cds.dcc9295c](https://doi.org/10.24381/cds.dcc9295c).
