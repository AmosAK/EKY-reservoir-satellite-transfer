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
| `06_crosssensor_agreement.py` | Same-day Sentinel-2/Landsat agreement (Table 11) | `results/crosssensor_agreement.csv` |
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

## Revision analyses (v1.1.0)

Added in response to editorial comments on the matchup window and on sensor harmonization.

| Script (`scripts/revision/`) | Purpose | Output (`results/`) |
|---|---|---|
| `fetch_times.py` | Downloads field sampling times (ActivityStartTime) from the Water Quality Portal for the training stations | `data/wqp_sampling_times.csv` |
| `h3_analysis.py` | Exact sampling-to-overpass gaps in hours; paired window test scored on identical ±3 h samples (Table 8); EKY test by training window; error vs time gap (Figure S1) | `window_paired_3h.csv`, `eky_by_window.csv`, `error_vs_gap_v2.csv`, `gap_coverage.csv` |
| `hls_analysis.py` | HLS S30/L30 paired comparison, EKY test with HLS models, same-day reflectance agreement | `hls_paired_comparison.csv`, `hls_eky_test.csv`, `hls_sameday_reflectance.csv` |
| `paired_analysis.py` | Six models (standard and HLS) on identical station-dates (Table 9); reflectance agreement on identical same-day pairs | `six_model_identical_samples.csv`, `reflectance_identical.csv` |

HLS matchups (`data/hls/`) were extracted with `earth_engine/hls_matchups.js` (HLS v2.0, `NASA/HLS/HLSS30/v002` and `NASA/HLS/HLSL30/v002`) using the same stations, window and screening rules as the main analysis.

Example (from the repository root):

```bash
python scripts/revision/hls_analysis.py data/hls data results
python scripts/revision/paired_analysis.py data/hls data results
python scripts/revision/h3_analysis.py data data/hls data/wqp_sampling_times.csv path/to/EKY_resultphyschem.csv results
```

`h3_analysis.py` also needs the EKY Water Quality Portal result file (sampling times for the EKY test stations), downloadable from https://www.waterqualitydata.us with organization 21KY_WQX.

## Changes in v1.2.1

- `scripts/revision/h3_analysis.py`: EKY sampling times are now matched by column name (they were mislabelled, so EKY test gaps were not computed), and the error-by-gap table is written with the Figure S1 bins to `results/error_vs_gap_v2.csv`.
- `results/gap_coverage.csv`: regenerated with the corrected script (two Sentinel-2 and four Landsat EKY test matchups within ±3 h; a third Sentinel-2 matchup falls at 3.0 h). Adds a ±6 h column.
- README: cross-sensor agreement is Table 11 of the article.

All other results are unchanged; a full rerun of `run_all.py` and the revision scripts reproduces them exactly.

## Revision analyses (v1.2.0): atmospheric correction and new test data

Added in response to the editorial comment on atmospheric correction. The Level-1 scenes behind every EKY test matchup (26 Sentinel-2 L1C, 23 Landsat 8/9 Collection 2 Level-1) were reprocessed with ACOLITE (https://github.com/acolite/acolite, commit 929063a), and the final models (trained on standard products) were applied, unchanged, to both versions of the same samples (Tables S5 and S6).

ACOLITE settings: dark spectrum fitting with `dsf_aot_estimation=fixed`, `dsf_residual_glint_correction=True`, `ancillary_data=False`, `s2_target_res=10`, `l2w_parameters=rhow_*`, and `limit` set to the extent of each scene's test stations ±0.05°.

| File | Contents |
|---|---|
| `earth_engine/eky_scene_list.js` | Finds the Sentinel-2 L1C / Landsat C2 L1 product IDs behind each EKY test matchup |
| `data/acolite/eky_test_scene_list.csv` | Output of the script above (57 matchups, 49 products) |
| `scripts/revision/acolite_get_s2.py` | Downloads a Sentinel-2 L1C SAFE from the Google Cloud public archive |
| `scripts/revision/acolite_get_landsat.py` | Downloads a Landsat C2 L1 bundle through the USGS M2M API (needs your own M2M token in `~/.usgs_m2m`, format `username token`) |
| `scripts/revision/acolite_pipeline.py` | Download → ACOLITE → 3 × 3-pixel water median at each station → delete scene |
| `scripts/revision/acolite_compare.py` | Band agreement and model scores, standard vs ACOLITE |
| `data/acolite/acolite_extract_*.csv` | ACOLITE ρs (and ρw) station medians |
| `data/acolite/acolite_eky_test_matchups.csv` | EKY test matchups with standard and ACOLITE (`*_aco`) reflectance |
| `results/acolite_band_agreement.csv`, `results/acolite_model_scores.csv` | Tables S5 and S6 |

```bash
python scripts/revision/acolite_pipeline.py WORKDIR S2 path/to/acolite   # then L for Landsat
python scripts/revision/acolite_compare.py
```

The ACOLITE L2R/L2W NetCDF files (about 3 GB) are not included; the pipeline regenerates them.

Additional test data: the only public in-lake data from an EKY reservoir after 2023 are U.S. Army Corps of Engineers (Nashville District) samples from Martins Fork Lake in 2022 (Water Quality Portal, organization `USACOEND`). `earth_engine/martinsfork_matchups.js` extracts their matchups (`data/martins_fork/`); only three matchups from one sampling date resulted, too few to test.

## Supplementary materials

`supplementary/Supplementary_Materials.docx` contains Tables S1–S6 and Figure S1 of the article. Their source files:

| Item | Source |
|---|---|
| Table S1 (total iron) | `results/total_iron_cv.csv` |
| Table S2 (leave-one-reservoir-out) | `results/eky_test_leave_one_reservoir_out.csv` |
| Table S3 (alternative learners, independent analyst) | `results/alternative_learners_cv.csv` |
| Table S4 (coalfield holdout) | `results/robustness_checks.csv` |
| Figure S1 (error vs time gap) | `results/error_vs_gap_v2.csv`, `figures/FigureS1_error_vs_gap.png` |
