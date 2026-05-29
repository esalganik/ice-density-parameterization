# Arctic sea-ice density parameterization

This repository contains MATLAB scripts used to process observational sea-ice density datasets, generate a merged Arctic sea-ice density database, derive observation-based density parameterizations, and reproduce the figures presented in:

> Salganik, E., Divine, D., Landy, J. C., Kauker, F., Nicolaus, M., & Stroeve, J.  
> *An observation-based year-round parameterization of Arctic sea-ice density.*  
> Geophysical Research Letters (submitted).

---

## Repository structure

```text
repository/
├── run_all.m
├── scripts/
│   ├── process/
│   ├── analysis/
│   │   ├── manuscript_statistics.m
│   │   ├── analyze_density_model.m
│   │   └── analyze_snow_dependence.m
│   └── helpers/
├── data/
│   ├── raw/
│   ├── processed/
│   ├── final/
│   ├── model/
│   └── colormaps/
├── external/
│   ├── m_map/
│   └── topo/
├── figures/
├── results/
│   └── manuscript_statistics.txt
├── README.md
└── LICENSE
```

---

## Workflow

The complete processing and analysis workflow can be executed using:

```matlab
run_all
```

The workflow consists of:

1. Import and processing of observational datasets
2. Merging of processed datasets
3. Addition of SnowModel-LG snow products
4. Generation of manuscript statistics
5. Generation of manuscript figures
6. Derivation of sea-ice density parameterizations

Large auxiliary datasets (ETOPO1 topography and SnowModel-LG files) are optional for lightweight reruns if previously processed outputs already exist.

---

## Optional large datasets

Some processing steps depend on large external datasets:

- ETOPO1 topography (`~445 MB`)
- SnowModel-LG NetCDF files

If these files are unavailable, the workflow can still run using previously generated intermediate products:

- `figures/Map_density_lab.png`
- `data/final/snow_model_matchups_with_models.mat`

This allows lightweight reruns without downloading all large external files.

---

## Source datasets

This repository uses publicly available observational and model datasets from multiple Arctic field campaigns and reanalysis products.

Raw observational datasets should be placed in:

```text
data/raw/
```

Model snow products should be placed in:

```text
data/model/
```

---

### Observational datasets

| Repository dataset | Source | DOI / Link |
|---|---|---|
| CONTRASTS | Divine et al. (2026), first- and second-year sea-ice density during CONTRASTS expedition | https://doi.org/10.1594/PANGAEA.993687 |
| GoNorth | Salganik et al. (2023), first- and second-year sea-ice density during GoNorth leg 1 | https://doi.org/10.1594/PANGAEA.962567 |
| MOSAiC FYI | Oggier et al. (2023a), MOSAiC first-year ice dataset | https://doi.org/10.1594/PANGAEA.956732 |
| MOSAiC SYI | Oggier et al. (2023b), MOSAiC second-year ice dataset | https://doi.org/10.1594/PANGAEA.959830 |
| MOSAiC Leg 5 | Salganik et al. (2024), MOSAiC Leg 5 sea-ice density dataset | https://doi.org/10.1594/PANGAEA.971266 |
| Nansen Legacy | Divine et al. (2025), Barents Sea and Arctic Basin sea-ice properties | https://doi.org/10.1002/gdj3.70001 |
| Svalbard fjord ice observations | Kongsfjorden, Svalbard sea-ice density observations | See manuscript references |
| SUDARCO | Additional observational dataset used in this study | See manuscript references |

---

### Snow products

| Product | Source | DOI / Link |
|---|---|---|
| SnowModel-LG (ERA5 forcing) | Liston et al. (2021) | https://doi.org/10.5067/27A0P5M6LZBI |
| SnowModel-LG (MERRA-2 forcing) | Liston et al. (2021) | https://doi.org/10.5067/27A0P5M6LZBI |

---

### Additional datasets

| Dataset | Source | DOI / Link |
|---|---|---|
| CryoSat-2 sea-ice thickness | Landy and Dawson (2022) | https://doi.org/10.5285/D8C66670-57AD-44FC-8FEF-942A46734ECB |

---

## External topography data

`plot_density_map.m` requires the ETOPO1 global topography dataset:

```text
etopo1_ice_g_i2.bin
```

Download from:

- https://data.meereisportal.de/data/topography/etopo1_ice_g_i2/
- https://www.ncei.noaa.gov/products/etopo-global-relief-model

Place the downloaded file in:

```text
external/topo/
```

The topography file is not included in this repository because of its size (~445 MB).

---

## MATLAB requirements

Tested with MATLAB R2024b.

Required MATLAB products and toolboxes:

- MATLAB
- Mapping Toolbox
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox

External MATLAB packages:

- M_Map mapping package

The repository includes a lightweight redistributed version of `m_map` required for map generation.

---

## Colormaps

Scientific colormaps (`lipari`, `buda`, `deep`) follow recommendations from:

> Crameri, F., Shephard, G. E., & Heron, P. J. (2020).  
> *The misuse of colour in science communication.*  
> Nature Communications, 11, 5444.  
> https://doi.org/10.1038/s41467-020-19160-7

and

> Thyng, K.M., Greene, C.A., Hetland, R.D., Zimmerle, H.M., & DiMarco, S.F. (2016).  
> *True colors of oceanography: Guidelines for effective and accurate colormap selection.*  
> Oceanography, 29(3), 9–13.  
> https://doi.org/10.5670/oceanog.2016.66

---

## Figures

The repository reproduces:

- Figure 1: Arctic sea-ice density map and temperature-thickness parameterization
- Figure 2: Snow-thickness and freeboard-dependent density parameterizations

Generated figures are exported to:

```text
figures/
```

---

## Citation

If you use this repository, please cite:

> Salganik, E., Divine, D., Landy, J. C., Kauker, F., Nicolaus, M., & Stroeve, J.  
> *An observation-based year-round parameterization of Arctic sea-ice density.*  
> Geophysical Research Letters (submitted).

---

## License

Please see `LICENSE` for repository licensing information.

External datasets and packages retain their original licenses and citation requirements.
