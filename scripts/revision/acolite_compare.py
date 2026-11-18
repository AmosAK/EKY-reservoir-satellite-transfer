"""Compare ACOLITE surface reflectance with the standard products on the EKY test matchups,
and score the regional models with ACOLITE inputs."""
import sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
import os
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# usage: python scripts/revision/acolite_compare.py   (reads data/acolite/, writes results/ and data/acolite/)
sys.path.insert(0, f"{REPO}/scripts")
from common import FEATURES, RF_PARAMS, N_TREES_FINAL, load_training, model_table, metrics
BANDS = {"S2": ["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11"],
         "L": ["coastal", "blue", "green", "red", "nir", "swir1"]}
TEST = {"S2": "eky_test_sentinel2.csv", "L": "eky_test_landsat89.csv"}
rows_band, rows_model, merged_all = [], [], []
for sen in ["S2", "L"]:
    try: a = pd.read_csv(f"{REPO}/data/acolite/acolite_extract_{'sentinel2' if sen == 'S2' else 'landsat89'}.csv").dropna(subset=["station"])
    except FileNotFoundError: continue
    rc = sorted([c for c in a if c.startswith("rhos_")], key=lambda c: int(c.split("_")[1]))
    # ACOLITE wavelengths differ slightly between S2A/S2B and L8/L9 -> map by band order
    a = a.copy()
    def to_std(r):
        w = sorted([c for c in r.index if c.startswith("rhos_") and pd.notna(r[c])], key=lambda c: int(c.split("_")[1]))
        return pd.Series({b: r[c] for b, c in zip(BANDS[sen], w)})
    A = pd.concat([a[["station", "img_date", "product", "n_px"]], a.apply(to_std, axis=1)], axis=1)
    A = A[A.n_px >= 5].groupby(["station", "img_date"])[BANDS[sen]].mean().reset_index()   # average overlapping tiles
    e = pd.read_csv(f"{REPO}/data/{TEST[sen]}")
    e["img_date"] = (pd.to_datetime(e.date) + pd.to_timedelta(np.floor(e.dt_days), "D")).dt.strftime("%Y-%m-%d")
    m = e.merge(A, on=["station", "img_date"], suffixes=("", "_aco"))
    m["sensor"] = sen; merged_all.append(m)
    for b in BANDS[sen]:
        if b + "_aco" not in m or b not in m: continue
        x, y = m[b].values, m[b + "_aco"].values; ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        x, y = x[ok], y[ok]
        rows_band.append(dict(sensor=sen, band=b, n=len(x), standard_median=np.median(x), acolite_median=np.median(y),
                              median_ratio=np.median(y / x), median_abs_diff=np.median(np.abs(y - x)),
                              r2=np.corrcoef(x, y)[0, 1] ** 2, rho=spearmanr(x, y)[0]))
    # model scoring: final models (trained on standard products) applied to both inputs on identical samples
    feat = FEATURES["S2" if sen == "S2" else "L89"]
    tr = load_training("S2" if sen == "S2" else "L89")
    Xs = feat(m); Xa = feat(m.assign(**{b: m[b + "_aco"] for b in BANDS[sen] if b + "_aco" in m and b in m}))
    for t in ["turb_ntu", "secchi_m", "tss_mgl", "chla_ugl"]:
        _, X, y = model_table(tr, t, feat)
        rf = RandomForestRegressor(N_TREES_FINAL, **RF_PARAMS).fit(X, y)
        ok = (m[t] > 0).values & Xs.notna().all(axis=1).values & Xa.notna().all(axis=1).values
        if ok.sum() < 5: continue
        yo = np.log10(m[t].values[ok])
        ms, ma = metrics(yo, rf.predict(Xs[ok])), metrics(yo, rf.predict(Xa[ok]))
        rows_model.append(dict(sensor=sen, target=t, n=int(ok.sum()), reservoirs=m.lake[ok].nunique(),
                               rho_standard=ms["rho"], rho_acolite=ma["rho"], MAPE_standard=ms["MAPE"], MAPE_acolite=ma["MAPE"],
                               bias_standard=ms["bias_ratio"], bias_acolite=ma["bias_ratio"]))
B = pd.DataFrame(rows_band).round(4); M = pd.DataFrame(rows_model).round(3)
B.to_csv(f"{REPO}/results/acolite_band_agreement.csv", index=False); M.to_csv(f"{REPO}/results/acolite_model_scores.csv", index=False)
pd.concat(merged_all).to_csv(f"{REPO}/data/acolite/acolite_eky_test_matchups.csv", index=False)
print(B.to_string(index=False)); print(M.to_string(index=False))
