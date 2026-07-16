\# SnowModel-LG model data



This directory contains SnowModel-LG snow products used in the analysis.



The large NetCDF files are not included in this repository because of their size.



The expected folder structure is:



`data/model/ERA5/`



and



`data/model/MERRA2/`



\## ERA5-forced SnowModel-LG



Place ERA5-forced SnowModel-LG files in:



`data/model/ERA5/`



The analysis uses snow-depth and snow-density NetCDF files.



Typical filenames include:



`SM\_snod\_ERA5\_\*.nc`



`SM\_sden\_ERA5\_\*.nc`



\## MERRA-2-forced SnowModel-LG



Place MERRA-2-forced SnowModel-LG files in:



`data/model/MERRA2/`



The analysis uses snow-depth and snow-density NetCDF files.



Typical filenames include:



`SM\_snod\_MERRA2\_\*.nc`



`SM\_sden\_MERRA2\_\*.nc`



\## Download



SnowModel-LG products are available from the NSIDC data archive:



https://doi.org/10.5067/27A0P5M6LZBI



After downloading, place the files in the corresponding ERA5 or MERRA2 folder.



\## Use in this repository



The SnowModel-LG files are required to recreate the model matchups from the original NetCDF products.



Most analysis and figure scripts can be rerun without the full NetCDF archive when the processed matchup file is already available:



`data/final/snow\_model\_matchups\_with\_models.mat`



The large NetCDF files should be excluded from version control through `.gitignore`.

