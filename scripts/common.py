"""Shared settings, feature definitions and metrics.

All scripts are run from the repository root, for example:
    python scripts/01_regional_cv.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

TARGETS = {
    "turb_ntu": "Turbidity",
    "secchi_m": "Secchi depth",
    "tss_mgl": "TSS",
    "chla_ugl": "Chlorophyll-a",
}
SENSORS = {"S2": "Sentinel-2", "L89": "Landsat 8/9"}

# Random forest settings used throughout the study (no further tuning).
RF_PARAMS = dict(min_samples_leaf=5, max_features=0.5, n_jobs=-1, random_state=1)
N_TREES_CV = 300      # trees per cross-validation fold
N_TREES_FINAL = 500   # trees in the final model applied to the EKY test set


# ----------------------------------------------------------------------------
# Features
# ----------------------------------------------------------------------------
def features_s2(d: pd.DataFrame) -> pd.DataFrame:
    """Sentinel-2 inputs: bands B1-B8A, NDCI, band ratios and NDTI."""
    X = d[["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A"]].copy()
    X["ndci"] = (d.B5 - d.B4) / (d.B5 + d.B4)
    X["b5b4"] = d.B5 / d.B4
    X["b3b2"] = d.B3 / d.B2
    X["b4b3"] = d.B4 / d.B3
    X["b2b4"] = d.B2 / d.B4
    X["b8b4"] = d.B8 / d.B4
    X["ndti"] = (d.B4 - d.B3) / (d.B4 + d.B3)
    X["b6b5"] = d.B6 / d.B5
    return X.replace([np.inf, -np.inf], np.nan)


def _landsat_like(blue, green, red, nir) -> pd.DataFrame:
    X = pd.DataFrame({"blue": blue, "green": green, "red": red, "nir": nir})
    X["br"] = blue / red
    X["gr"] = green / red
    X["bg"] = blue / green
    X["nr"] = nir / red
    X["ndti"] = (red - green) / (red + green)
    return X.replace([np.inf, -np.inf], np.nan)


def features_l89(d: pd.DataFrame) -> pd.DataFrame:
    """Landsat 8/9 inputs: blue, green, red, NIR, band ratios and NDTI."""
    return _landsat_like(d.blue, d.green, d.red, d.nir)


def features_s2_landsat_bands(d: pd.DataFrame) -> pd.DataFrame:
    """Sentinel-2 restricted to Landsat-equivalent bands (B2, B3, B4, B8A)."""
    return _landsat_like(d.B2, d.B3, d.B4, d.B8A)


FEATURES = {"S2": features_s2, "L89": features_l89}


# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------
def load_training(sensor: str) -> pd.DataFrame:
    name = {"S2": "training_sentinel2.csv", "L89": "training_landsat89.csv"}[sensor]
    return pd.read_csv(DATA / name)


def load_eky_test(sensor: str) -> pd.DataFrame:
    name = {"S2": "eky_test_sentinel2.csv", "L89": "eky_test_landsat89.csv"}[sensor]
    return pd.read_csv(DATA / name)


def model_table(d: pd.DataFrame, target: str, featurizer):
    """Rows with a positive target and complete features; returns (rows, X, log10 y)."""
    a = d[d[target] > 0]
    X = featurizer(a)
    ok = X.notna().all(axis=1)
    a, X = a[ok].reset_index(drop=True), X[ok].reset_index(drop=True)
    return a, X, np.log10(a[target].values)


# ----------------------------------------------------------------------------
# Metrics (computed on log10 values)
# ----------------------------------------------------------------------------
def metrics(y_log: np.ndarray, p_log: np.ndarray) -> dict:
    r = p_log - y_log
    return {
        "n": len(y_log),
        "R2": 1 - np.sum(r ** 2) / np.sum((y_log - y_log.mean()) ** 2),
        "MAPE": 100 * np.median(np.abs(10 ** r - 1)),
        "bias_ratio": np.median(10 ** r),
        "rho": spearmanr(y_log, p_log)[0],
    }
