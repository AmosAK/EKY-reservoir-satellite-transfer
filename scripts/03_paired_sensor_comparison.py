"""Sentinel-2 versus Landsat 8/9 on identical field samples (Tables 3 and 4).

Three models per parameter on the same station-grouped folds:
Landsat 8/9, Sentinel-2 restricted to Landsat-equivalent bands, Sentinel-2 with all
bands. Station-averaged absolute log errors are compared with the Wilcoxon
signed-rank test; R2 confidence intervals come from a station-block bootstrap.
The comparison is repeated for samples imaged by both sensors on the same day.
Outputs:
    results/paired_comparison_all.csv
    results/paired_comparison_same_day.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold

from common import (DATA, N_TREES_CV, RESULTS, RF_PARAMS, TARGETS, _landsat_like,
                    features_s2, features_s2_landsat_bands, metrics)

paired = pd.read_csv(DATA / "training_paired_s2_landsat.csv")
MODELS = {
    "Landsat 8/9": lambda d: _landsat_like(d.L_blue, d.L_green, d.L_red, d.L_nir),
    "S2, Landsat bands only": features_s2_landsat_bands,
    "S2, all bands": features_s2,
}


def compare(c: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target, name in TARGETS.items():
        a = c[c[target] > 0].copy()
        ok = np.ones(len(a), bool)
        for f in MODELS.values():
            ok &= f(a).notna().all(axis=1).values
        a = a[ok].reset_index(drop=True)
        y, g = np.log10(a[target].values), a.station.values
        folds = list(GroupKFold(5).split(a, y, g))
        errors = {}
        for model, f in MODELS.items():
            X = f(a).values
            p = np.zeros(len(y))
            for tr, te in folds:
                p[te] = RandomForestRegressor(N_TREES_CV, **RF_PARAMS).fit(X[tr], y[tr]).predict(X[te])
            m = metrics(y, p)
            errors[model] = np.abs(p - y)
            stations = np.unique(g)
            idx = {s: np.where(g == s)[0] for s in stations}
            rng = np.random.default_rng(0)
            boot = []
            for _ in range(500):
                ii = np.concatenate([idx[s] for s in rng.choice(stations, len(stations))])
                boot.append(metrics(y[ii], p[ii])["R2"])
            rows.append({"parameter": name, "model": model, "n": len(y), "stations": len(stations),
                         "R2": m["R2"], "R2_lo": np.percentile(boot, 2.5), "R2_hi": np.percentile(boot, 97.5),
                         "MAPE": m["MAPE"], "rho": m["rho"]})
        st = pd.DataFrame({"g": g, **errors}).groupby("g").mean()
        rows[-1]["p_vs_landsat"] = wilcoxon(st["S2, all bands"], st["Landsat 8/9"]).pvalue
        rows[-1]["pct_stations_s2_better"] = 100 * (st["S2, all bands"] < st["Landsat 8/9"]).mean()
        print(rows[-3:])
    return pd.DataFrame(rows)


compare(paired).round(3).to_csv(RESULTS / "paired_comparison_all.csv", index=False)
same_day = paired[paired.img_gap_days.abs() == 0].reset_index(drop=True)
compare(same_day).round(3).to_csv(RESULTS / "paired_comparison_same_day.csv", index=False)
