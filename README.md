# s2sExposure_ms
This repository contains the code to run an exposure model from subseasonal-to-seasonal forecasts data for an example species, and to run the metaranking for regional prioritization. It also contains the data outptuts of the analysis aggregated at species groups for manuscript reproducibility.

Full data can be obtained here:
Species expert range maps from IUCN are available at available at https://www.iucnredlist.org/resources/spatial-data-download. The GEOS-S2S-V2 data are available on the Discover server of NCCS at https://www.nccs.nasa.gov/systems/ data-portal, and GEOS-S2S-V2 forecasts output data are presently available at https://gmao.gsfc.nasa.gov/gmaoftp/gmaofcst/. Historical climate reanalysis data from ERA5 are available at https://cds.climate.copernicus.eu/ ERA5. Administrative units available at gadm.org v. 4.1.  Climate exposure results are synthetized by large groups (e.g. mammals, birds, amphibians, reptiles) to respect species expert ranges licensing and protect vulnerable species. These and resulting data are available at a dedicated GitHub repository https://github.com/pepbioalerts/s2sExposure_ms

The /code folder contains all scripts needed to run the forecats: It contains the functions (funcs_*.R) as well as the exammple workflow for the exposure forecasts (workflowExample.R) and the ranking procedure (workflowRanking.R).

