# CryoSat-2 data

This directory contains the CryoSat-2 NetCDF files used to reproduce the satellite-based analyses in this repository.

The NetCDF files are not included in this repository because of their size.

The expected folder is:

`data/cs2/`

## Required files

Download and place the following files in:

`data/cs2/`

Required filenames:

`uit_sit_v2northpolarstereo_80km_feb_2011-2023_v3p0.nc`

`uit_sit_v2northpolarstereo_80km_jul_2011-2023_v3p0.nc`

## Download

The CryoSat-2 fields used in this study are archived on Zenodo:

https://doi.org/10.5281/zenodo.21839354

After downloading, place both NetCDF files directly in:

`data/cs2/`

## Use in this repository

The CryoSat-2 NetCDF files are required to reproduce:

- Figure 3
- Figure S5
- the sea-ice volume and thickness-difference checks in `scripts/checks/cs2_manuscript_checks.m`

The files should be excluded from version control through `.gitignore`.