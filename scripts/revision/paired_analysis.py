"""Paired versions of the window and HLS checks.
(1) Window test: same test samples (same-day matchups at held-out stations),
    models trained with +/-3 d, +/-1 d or same-day matchups.
(2) Six-model sensor comparison on IDENTICAL samples: original S2 L2A / Landsat C2
    and HLS S30 / L30, all on the station-dates present in both paired datasets.
(3) Same-day reflectance agreement on identical station-dates.
"""
import sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold

H, R, OUT = sys.argv[1:4]

RF = dict(min_samples_leaf=5, max_features=0.5, n_jobs=-1, random_state=1)
T = {"turb_ntu": "Turbidity", "secchi_m": "Secchi depth", "tss_mgl": "TSS", "chla_ugl": "Chlorophyll-a"}
import hls_common as hc
hc.H, hc.R = H, R
load, f_s30 = hc.load, hc.f_s30

def r2(y, p): return 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)
def mape(y, p): return 100 * np.median(np.abs(10 ** (p - y) - 1))

def s2feat(d):
    X = d[["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A"]].copy()
    X["ndci"] = (d.B5 - d.B4) / (d.B5 + d.B4); X["b5b4"] = d.B5 / d.B4; X["b3b2"] = d.B3 / d.B2
    X["b4b3"] = d.B4 / d.B3; X["b2b4"] = d.B2 / d.B4; X["b8b4"] = d.B8 / d.B4
    X["ndti"] = (d.B4 - d.B3) / (d.B4 + d.B3); X["b6b5"] = d.B6 / d.B5
    return X.replace([np.inf, -np.inf], np.nan)
def land(b, g, r, n):
    X = pd.DataFrame({"blue": b, "green": g, "red": r, "nir": n, "br": b / r, "gr": g / r, "bg": b / g, "nr": n / r,
                      "ndti": (r - g) / (r + g)})
    return X.replace([np.inf, -np.inf], np.nan)

# ---------------- (2) six models on identical samples ----------------
s30, _ = load("S30"); l30, _ = load("L30")
hp = s30.merge(l30, on=["station", "date"], suffixes=("", "_L"))
hp["hgap"] = (pd.to_datetime(hp.img_date) - pd.to_datetime(hp.img_date_L)).dt.days.abs()
op = pd.read_csv(f"{R}/training_paired_s2_landsat.csv")
both = op.merge(hp[["station", "date", "hgap", "coastal", "blue", "green", "red", "re1", "re2", "re3", "nir",
                    "blue_L", "green_L", "red_L", "nir_L"]], on=["station", "date"], how="inner")
print("identical station-dates in original and HLS paired sets:", len(both), "same-day in both:",
      ((both.img_gap_days == 0) & (both.hgap == 0)).sum())
S = []
for sub, B in [("all", both), ("same_day_both", both[(both.img_gap_days == 0) & (both.hgap == 0)])]:
    for t, name in T.items():
        d = B[B[t] > 0]
        F = {"Landsat C2 L2": land(d.L_blue, d.L_green, d.L_red, d.L_nir),
             "S2 L2A, Landsat bands": land(d.B2, d.B3, d.B4, d.B8A),
             "S2 L2A, all bands": s2feat(d),
             "HLS L30": land(d.blue_L, d.green_L, d.red_L, d.nir_L),
             "HLS S30, Landsat bands": land(d.blue, d.green, d.red, d.nir),
             "HLS S30, all bands": f_s30(d)}
        ok = np.ones(len(d), bool)
        for X in F.values(): ok &= X.notna().all(1).values
        d = d[ok]; y = np.log10(d[t].values); g = d.station.values
        if len(d) < 60: continue
        P = {}
        for k, X in F.items():
            X = X[ok].values; p = np.zeros(len(y))
            for a, b in GroupKFold(5).split(X, y, g):
                p[b] = RandomForestRegressor(300, **RF).fit(X[a], y[a]).predict(X[b])
            P[k] = p
        E = {k: pd.Series(np.abs(p - y)).groupby(g).mean() for k, p in P.items()}
        pv = lambda a, b: round(wilcoxon(E[a], E[b]).pvalue, 4)
        for k, p in P.items():
            S.append(dict(subset=sub, parameter=name, model=k, n=len(y), stations=len(set(g)),
                          R2=round(r2(y, p), 2), MAPE=round(mape(y, p), 1), rho=round(spearmanr(y, p)[0], 2)))
        S.append(dict(subset=sub, parameter=name, model="p: S2all vs Landsat | HLS S30all vs L30 | S2 LB vs L | S30 LB vs L30",
                      n=len(y), R2=f"{pv('S2 L2A, all bands','Landsat C2 L2')} | {pv('HLS S30, all bands','HLS L30')} | "
                                   f"{pv('S2 L2A, Landsat bands','Landsat C2 L2')} | {pv('HLS S30, Landsat bands','HLS L30')}"))
        print(sub, name, "done", flush=True)
S = pd.DataFrame(S); S.to_csv(f"{OUT}/six_model_identical_samples.csv", index=False)
print(S.to_string(index=False))

# ---------------- (3) reflectance agreement, identical station-dates ----------------
sd = both[(both.img_gap_days == 0) & (both.hgap == 0)]
C = []
for band, s2b in [("blue", "B2"), ("green", "B3"), ("red", "B4"), ("nir", "B8A")]:
    ok = (sd["L_" + band] > 0) & (sd[s2b] > 0) & (sd[band + "_L"] > 0) & (sd[band] > 0)
    ro = sd.loc[ok, "L_" + band] / sd.loc[ok, s2b]; rh = sd.loc[ok, band + "_L"] / sd.loc[ok, band]
    C.append(dict(band=band, n=int(ok.sum()), orig_median_ratio=round(ro.median(), 3),
                  orig_median_absdiff_pct=round(100 * np.median(np.abs(ro - 1)), 1),
                  hls_median_ratio=round(rh.median(), 3), hls_median_absdiff_pct=round(100 * np.median(np.abs(rh - 1)), 1),
                  p_absdiff=round(wilcoxon(np.abs(np.log(ro)), np.abs(np.log(rh))).pvalue, 4)))
C = pd.DataFrame(C); C.to_csv(f"{OUT}/reflectance_identical.csv", index=False); print(C.to_string(index=False))
