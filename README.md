# Transferring Sentinel-2 and Landsat 8/9 water-quality retrievals to Appalachian coalfield reservoirs

Data and code for:

> Akwemoh, A. F., Gyawali, B. R., Basnet, U., & Zourarakis, D. P. Transferring Sentinel-2 and Landsat 8/9 water-quality retrievals to data-scarce Appalachian coalfield reservoirs: A withheld-reservoir test. *Remote Sensing* (submitted).

Random forest models for turbidity, Secchi depth, total suspended solids (TSS) and chlorophyll-a were trained on satellite–in situ matchups from lakes in seven U.S. states (KY, VA, NC, OH, PA, TN, WV) and tested on 12 withheld reservoirs in and adjacent to the Eastern Kentucky Coal Field. This repository contains the screened matchup data, the Google Earth Engine scripts used to extract them, and Python scripts that reproduce the statistics and figures in the article.

## Repository structure

```
data/            screened matchup datasets, EKY test data, 2017–2026 retrievals
gis/             reservoir watersheds (USGS StreamStats, EPSG:5070)
earth_engine/    Google Earth Engine scripts (matchup extraction and time series)
scripts/         analysis scripts, run in numbered order
results/         output tables written by the scripts
figures/         output figures written by the scripts (600 dpi)
run_all.py       runs scripts 01–07 in order
```

## Quick start

```bash
pip install -r requirements.txt
python run_all.py
```

The full run takes about 15–30 minutes on a laptop. Each script can also be run on its own from the repository root, e.g. `python scripts/02_eky_transfer_test.py`.

The coalfield-holdout check in `05_robustness_checks.py` is optional. It needs the Kentucky Geological Survey physiographic regions shapefile (https://www.uky.edu/KGS/gis/regions.htm), which is not redistributed here:

```bash
python run_all.py --coalfield path/to/regions.shp
```

## Scripts and the results they reproduce

| Script | Reproduces | Output |
|---|---|---|
| `01_regional_cv.py` | Station-grouped five-fold cross-validation (Table 2, regional columns); lake-grouped sensitivity | `results/regional_cv.csv` |
| `02_eky_transfer_test.py` | Withheld EKY test (Table 2), bootstrap confidence intervals (Table 5), leave-one-reservoir-out and within-reservoir ρ (Table S2) | `results/eky_test_*.csv` |
| `03_paired_sensor_comparison.py` | Sentinel-2 vs Landsat 8/9 on identical samples (Tables 3 and 4) | `results/paired_comparison_*.csv` |
| `04_permutation_importance.py` | Permutation importance (Figure 6) | `results/permutation_importance.csv` |
| `05_robustness_checks.py` | Cloud limits (Table 7), state-balanced and Appalachian-only training, coalfield holdout (Table S4) | `results/robustness_checks.csv` |
| `06_crosssensor_agreement.py` | Same-day Sentinel-2/Landsat agreement (Table 9) | `results/crosssensor_agreement.csv` |
| `07_figures.py` | Figures 4–8 | `figures/*.png` |

All scripts use fixed random seeds, so repeated runs give identical results. Point estimates (R², ρ, MAPE, bias ratio) match the article. Bootstrap confidence bounds are Monte Carlo estimates and can differ from the published values by a few hundredths (more for TSS, which has only 8–10 test samples). Small differences can also occur across scikit-learn versions.

## Data files

All concentrations are in their original units; models use log10-transformed values.

| File | Description |
|---|---|
| `training_sentinel2.csv` | 7,109 Sentinel-2 Level-2A matchups at 992 stations (2017–2026) |
| `training_landsat89.csv` | 6,200 Landsat 8/9 Collection 2 Level-2 matchups at 941 stations (2013–2026) |
| `training_paired_s2_landsat.csv` | 2,695 field samples with a valid matchup from both sensors |
| `eky_test_sentinel2.csv`, `eky_test_landsat89.csv` | Withheld EKY test matchups (29 Sentinel-2, 28 Landsat 8/9) with screened field values |
| `eky_field_visits_screened.csv` | Screened EKY field visits, 2015–2023 |
| `eky_screening_log.csv` | Records removed during screening, by reason |
| `eky_stations.csv` | The 37 EKY monitoring stations |
| `eky_reservoirs.csv` | Reservoir characteristics (Table 1) |
| `eky_timeseries_sentinel2.csv`, `eky_timeseries_landsat89.csv` | Retrievals at EKY stations, 2017–2026: median and 10th–90th percentile of tree predictions, and an in-domain flag for each parameter |

**Key columns**

- `station`, `lake`, `state`, `latitude`, `longitude`: Water Quality Portal monitoring location and its coordinates (WGS84). `state` uses federal state codes (21 = KY, 37 = NC, 39 = OH, 42 = PA, 47 = TN, 51 = VA, 54 = WV).
- `date`: field sampling date. `img_date`: image acquisition date. `dt_days`: image minus sampling date, in days.
- `fold`: fixed station-grouped fold (0–4), used by `04_permutation_importance.py`.
- `B1`–`B12`: Sentinel-2 surface reflectance. `coastal`, `blue`, `green`, `red`, `nir`, `swir1`, `swir2`: Landsat 8/9 surface reflectance. Values are the median of a 3 × 3 pixel window over water pixels (MNDWI > 0).
- `cloud_pct`: scene cloud cover (%). `n_px`: number of valid water pixels.
- `turb_ntu` (NTU), `secchi_m` (m), `tss_mgl` (mg/L), `chla_ugl` (µg/L), `fe_mgl` (total iron, mg/L): field measurements.
- `L_*` (paired file): the Landsat 8/9 values for the same field sample; `img_gap_days`: days between the two images.
- `*_in_domain` (time series): TRUE when every model input lies within the 1st–99th percentile of the training data.

## Data sources

- Field data: Water Quality Portal (Read et al., 2017; https://www.waterqualitydata.us). Kentucky Division of Water records are stored under organization code `21KY_WQX`.
- U.S. EPA National Lakes Assessment 2017 and 2022.
- Satellite imagery: Copernicus Sentinel-2 Level-2A (`COPERNICUS/S2_SR_HARMONIZED`) and USGS Landsat 8/9 Collection 2 Level-2, accessed through Google Earth Engine.
- Watersheds: USGS StreamStats.

## Earth Engine scripts

Upload the station-visit table as an Earth Engine asset and replace `projects/YOUR_PROJECT/assets/...` at the top of each script with your asset path.

- `regional_matchups.js`, `nla_matchups.js`: training matchups
- `eky_matchups.js`: EKY test matchups
- `eky_timeseries.js`, `eky_landsat_timeseries.js`: 2017–2026 station time series

The AquaSat benchmark (Section 3.4) uses the public AquaSat dataset (Ross et al., 2019, https://doi.org/10.1029/2019WR024883), which is not redistributed here.

## Licence

Code: MIT (see `LICENSE`). Data: CC BY 4.0 (see `LICENSE-DATA.md`).

## Contact

Amos Francis Akwemoh, Center for Geospatial Intelligence and Environmental Security, Kentucky State University.
