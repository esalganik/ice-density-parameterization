# Arctic sea-ice density parameterization

This repository contains MATLAB scripts used to process observational sea-ice density datasets, generate a merged Arctic sea-ice density database, derive observation-based density parameterizations, and reproduce the figures presented in:

> Salganik, E., Divine, D., Landy, J. C., Kauker, F., Nicolaus, M., & Stroeve, J.
> *An observation-based year-round parameterization of Arctic sea-ice density.*
> Geophysical Research Letters (submitted).

The repository reproduces all manuscript figures, Supporting Information figures, summary statistics, empirical uncertainty estimates, and cross-validation tables from the submitted study.

---

## Repository structure

```text
repository/
├── run_all.m
├── scripts/
│   ├── process/
│   ├── analysis/
│   ├── checks/
│   ├── merge_processed_summaries.m
│   ├── add_SMLG_model_snow.m
│   └── plot_density_map.m
├── data/
│   ├── raw/
│   ├── processed/
│   ├── final/
│   ├── model/
│   ├── cs2/
│   └── colormaps/
├── external/
│   ├── m_map/
│   ├── ncpolarm/
│   ├── inpaint_nans/
│   ├── viridis/
│   └── topo/
├── figures/
├── results/
├── README.md
└── LICENSE
```

This is a selective overview rather than a complete file listing. Processing scripts are stored in `scripts/process/`, analysis and figure scripts in `scripts/analysis/`, and optional diagnostic tests in `scripts/checks/`. The main workflow is controlled by `run_all.m`.

---

## Workflow

Run the complete workflow from the repository root using:

```matlab
run_all
```

By default, `run_all.m` executes the complete processing and analysis workflow. Individual stages can be disabled using the workflow switches at the beginning of the script. The required raw and external datasets must be available in the directories described below.

The workflow consists of:

1. Import and processing of observational datasets
2. Merging of processed datasets
3. Addition of SnowModel-LG snow products
4. Generation of manuscript statistics
5. Derivation, cross-validation, and empirical uncertainty evaluation of the sea-ice density parameterizations
6. Generation of manuscript and supporting information figures

Large auxiliary datasets (ETOPO1 topography and SnowModel-LG files) are optional for lightweight reruns if previously processed outputs already exist.

---

## Additional checks

The scripts in `scripts/checks/` are standalone diagnostic analyses. They support specific statements or methodological choices in the manuscript but are not part of the default `run_all.m` workflow.

- `cs2_manuscript_checks.m` reproduces the CryoSat-2 sea-ice volume and thickness-difference values reported in the manuscript using the February and July 2011–2023 CryoSat-2 products.
- `check_temperature_thickness_independence.m` compares temperature-only and temperature-thickness models and tests whether thickness explains residual density variability.
- `snow_density_source_sensitivity_freeboard_only.m` tests the sensitivity of the freeboard-based parameterization to the source of snow density.
- `freezeup_skill_test.m` evaluates the influence of October and November FYI freeze-up observations on model skill.
- `analyze_density_model_script_no_svalbard.m` tests whether excluding Svalbard landfast-ice observations changes the fitted temperature-thickness parameterization.
- `temperature_threshold_by_ice_age.m` tests whether separate transition temperatures are supported for FYI and SYI/MYI.

---

## External and large datasets

Some processing steps depend on large external datasets:

- ETOPO1 topography (`~445 MB`)
- SnowModel-LG NetCDF files
- CryoSat-2 NetCDF files archived at https://doi.org/10.5281/zenodo.21839354

If the ETOPO1 or SnowModel-LG source files are unavailable, the workflow can still use previously generated intermediate products:

- `figures/Map_density_lab.png`
- `data/final/snow_model_matchups_with_models.mat`

This allows lightweight reruns without downloading the ETOPO1 and SnowModel-LG source files when the corresponding intermediate products already exist.

The CryoSat-2 NetCDF files must be downloaded from https://doi.org/10.5281/zenodo.21839354 and placed in `data/cs2/` to regenerate Figures 3 and S5 and to run `scripts/checks/cs2_manuscript_checks.m`.

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

CryoSat-2 products required for Figures 3 and S5 are archived at:

https://doi.org/10.5281/zenodo.21839354

Download the files and place them in:

```text
data/cs2/
```

The required files are:

```text
uit_sit_v2northpolarstereo_80km_feb_2011-2023_v3p0.nc
uit_sit_v2northpolarstereo_80km_jul_2011-2023_v3p0.nc
```

### Observational datasets

| Repository dataset | Source | DOI / Link |
|---|---|---|
| CONTRASTS | Divine et al. (2026), first- and second-year sea-ice density during CONTRASTS expedition | https://doi.org/10.1594/PANGAEA.993687 |
| GoNorth | Salganik et al. (2023), first- and second-year sea-ice density during GoNorth leg 1 | https://doi.org/10.1594/PANGAEA.962567 |
| MOSAiC FYI | Oggier et al. (2023a), MOSAiC first-year ice dataset | https://doi.org/10.1594/PANGAEA.956732 |
| MOSAiC SYI | Oggier et al. (2023b), MOSAiC second-year ice dataset | https://doi.org/10.1594/PANGAEA.959830 |
| MOSAiC Leg 5 | Salganik et al. (2024), MOSAiC Leg 5 sea-ice density dataset | https://doi.org/10.1594/PANGAEA.971266 |
| Nansen Legacy | Divine et al. (2025), Barents Sea and Arctic Basin sea-ice properties | https://doi.org/10.1002/gdj3.70001 |
| N-ICE2015–SUDARCO–Svalbard compilation | Combined sea-ice density observations from N-ICE2015, SUDARCO, and Kongsfjorden landfast ice | https://doi.org/10.21334/NPOLAR.2026.7B676D79 |

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
| CryoSat-2 fields used in this study    | Landy and Salganik (2026) | [https://doi.org/10.5281/zenodo.21839354](https://doi.org/10.5281/zenodo.21839354) |
| ETOPO1 topography | Amante and Eakins (2009) | https://doi.org/10.7289/V5C8276M |

---

## External topography data

Generation of the sea-ice density map requires the ETOPO1 global topography dataset:

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

External MATLAB packages included in this repository:

- `M_Map`
- `ncpolarm`
- `inpaint_nans`
- `viridis`

The repository includes redistributed copies of these packages under `external/`. Each package retains its original copyright notice and license in the corresponding subdirectory.

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

The repository reproduces the figures presented in the manuscript and Supporting Information:

Manuscript:

- Figure 1: Arctic sea-ice density map and temperature-thickness parameterization
- Figure 2: Snow-thickness and freeboard-dependent density parameterizations
- Figure 3: CryoSat-2 sea-ice density fields and corresponding sea-ice thickness distributions

Supporting Information:

- Figure S1: Hydrostatically derived effective sea-ice density during MOSAiC
- Figure S2: Directly measured MOSAiC ridge-core densities
- Figure S3: Snow-thickness and freeboard-dependent density parameterizations using MERRA-2 SM-LG
- Figure S4: Comparison of the sea-ice density formulations used in the CryoSat-2 sensitivity analysis
- Figure S5: CryoSat-2 sea-ice thickness differences between density parameterizations
- Figure S6: Evaluation of ERA5- and MERRA-2-forced SnowModel-LG snow thickness
- Table S2: Leave-one-campaign-out cross-validation statistics

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
