# ETOPO1 topography

This directory should contain the ETOPO1 global topography file:

`etopo1_ice_g_i2.bin`

The file is not included in this repository because it is approximately 445 MB.

## Download

Download the dataset from:

- https://data.meereisportal.de/data/topography/etopo1_ice_g_i2/
- https://www.ncei.noaa.gov/products/etopo-global-relief-model

After downloading, place the file at:

`external/topo/etopo1_ice_g_i2.bin`

## Use in this repository

The topography file is required only to regenerate the Arctic sea-ice density map.

Most analyses can be rerun without it when the previously generated map is available:

`figures/Map_density_lab.png`

The large binary file is excluded from version control through `.gitignore`.
