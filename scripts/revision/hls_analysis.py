"""HLS analysis: (A) paired S30 vs L30 comparison on identical samples,
(B) EKY withheld test with HLS-trained models, (C) same-day reflectance
agreement S30 vs L30 compared with the original S2 L2A vs Landsat C2 pairs."""
import sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold

H = sys.argv[1]           # folder with HLS csvs
R = sys.argv[2]           # repo data folder
OUT = sys.argv[3]
RF = dict(min_samples_leaf=5, max_features=0.5, n_jobs=-1, random_state=1)
T = {"turb_ntu": "Turbidity", "secchi_m": "Secchi depth", "tss_mgl": "TSS", "chla_ugl": "Chlorophyll-a"}
LAB = list(T) + ["fe_mgl"]

def clean(d):
    need = ["blue", "green", "red", "nir", "swir1"]
    d = d.dropna(subset=need)
    d = d[(d.swir1 <= 0.03) & (d.blue > 0) & (d.green > 0) & (d.red > 0) & (d.nir > -1)]
    d = d.assign(adt=d.dt_days.abs()).sort_values(["station", "date", "adt", "cloud_pct"])
    return d.drop_duplicates(["station", "date"]).drop(columns="adt")

def labels_regional():
    tr = pd.concat([pd.read_csv(f"{R}/training_sentinel2.csv"), pd.read_csv(f"{R}/training_landsat89.csv"),
                    pd.read_csv(f"{R}/training_paired_s2_landsat.csv")])
    return tr.groupby(["station", "date"])[LAB].first().reset_index()

def load(sensor):
    reg = pd.read_csv(f"{H}/HLS_{sensor}_REG_matchups.csv")
    reg = reg.drop(columns=[c for c in LAB if c in reg]).merge(labels_regional(), on=["station", "date"], how="inner")
    nla = pd.read_csv(f"{H}/HLS_{sensor}_NLA_matchups.csv")
    for c in LAB:
        if c not in nla: nla[c] = np.nan
    tr = clean(pd.concat([reg, nla], ignore_index=True))
    eky = pd.read_csv(f"{H}/HLS_{sensor}_EKY_matchups.csv")
    lab = pd.read_csv(f"{R}/eky_field_visits_screened.csv").groupby(["station", "date"])[LAB].first().reset_index()
    eky = clean(eky.drop(columns=[c for c in LAB if c in eky]).merge(lab, on=["station", "date"], how="inner"))
    return tr, eky

def f_land(d, p=""):
    b, g, r, n = d[p + "blue"], d[p + "green"], d[p + "red"], d[p + "nir"]
    X = pd.DataFrame({"blue": b, "green": g, "red": r, "nir": n, "br": b / r, "gr": g / r, "bg": b / g, "nr": n / r,
                      "ndti": (r - g) / (r + g)})
    return X.replace([np.inf, -np.inf], np.nan)

def f_s30(d, p=""):
    c = lambda k: d[p + k]
    X = pd.DataFrame({k: c(k) for k in ["coastal", "blue", "green", "red", "re1", "re2", "re3", "nir"]})
    X["ndci"] = (c("re1") - c("red")) / (c("re1") + c("red")); X["b5b4"] = c("re1") / c("red")
    X["b3b2"] = c("green") / c("blue"); X["b4b3"] = c("red") / c("green"); X["b2b4"] = c("blue") / c("red")
    X["b8b4"] = c("nir") / c("red"); X["ndti"] = (c("red") - c("green")) / (c("red") + c("green"))
    X["b6b5"] = c("re2") / c("re1")
    return X.replace([np.inf, -np.inf], np.nan)

def r2(y, p): return 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)
def mape(y, p): return 100 * np.median(np.abs(10 ** (p - y) - 1))
def bias(y, p): return 10 ** np.median(p - y)

s30, s30e = load("S30"); l30, l30e = load("L30")
print("HLS training matchups: S30", len(s30), "L30", len(l30), "| EKY: S30", len(s30e), "L30", len(l30e))
rows = []

# ---------------- (A) paired comparison on identical samples ----------------
pair = s30.merge(l30, on=["station", "date"], suffixes=("", "_L"))
pair["gap"] = (pd.to_datetime(pair.img_date) - pd.to_datetime(pair.img_date_L)).dt.days.abs()
print("paired samples:", len(pair), "same-day:", (pair.gap == 0).sum())
pairL = pair.copy()
for k in ["blue", "green", "red", "nir"]: pairL["L_" + k] = pair[k + "_L"]
A = []
for sub_name, P in [("all", pair), ("same_day", pair[pair.gap == 0])]:
    for t, name in T.items():
        d = P[P[t] > 0]
        feats = {"HLS L30": f_land(d.rename(columns={k + "_L": "L_" + k for k in ["blue", "green", "red", "nir"]}), "L_"),
                 "HLS S30, Landsat bands": f_land(d), "HLS S30, all bands": f_s30(d)}
        ok = np.ones(len(d), bool)
        for X in feats.values(): ok &= X.notna().all(1).values
        d = d[ok]; y = np.log10(d[t].values); g = d.station.values
        if len(d) < 50: continue
        preds = {}
        for m, X in feats.items():
            X = X[ok].values; p = np.zeros(len(y))
            for a, b in GroupKFold(5).split(X, y, g):
                p[b] = RandomForestRegressor(300, **RF).fit(X[a], y[a]).predict(X[b])
            preds[m] = p
        e = {m: pd.Series(np.abs(p - y)).groupby(g).mean() for m, p in preds.items()}
        pv = wilcoxon(e["HLS S30, all bands"], e["HLS L30"]).pvalue
        pv2 = wilcoxon(e["HLS S30, Landsat bands"], e["HLS L30"]).pvalue
        for m, p in preds.items():
            A.append(dict(subset=sub_name, parameter=name, model=m, n=len(y), stations=len(set(g)),
                          R2=round(r2(y, p), 2), MAPE=round(mape(y, p), 1), rho=round(spearmanr(y, p)[0], 2),
                          p_vs_L30=(round(pv, 4) if m == "HLS S30, all bands" else round(pv2, 4) if "Landsat bands" in m else None)))
A = pd.DataFrame(A); A.to_csv(f"{OUT}/hls_paired_comparison.csv", index=False)
print(A.to_string(index=False))

# ---------------- (B) EKY withheld test ----------------
B = []
for sensor, tr, te, fz in [("HLS S30", s30, s30e, f_s30), ("HLS L30", l30, l30e, f_land)]:
    for t, name in T.items():
        a = tr[tr[t] > 0]; Xa = fz(a); ok = Xa.notna().all(1); Xa, a = Xa[ok], a[ok]
        b = te[te[t] > 0]; Xb = fz(b); okb = Xb.notna().all(1); Xb, b = Xb[okb], b[okb]
        if len(b) < 5: B.append(dict(sensor=sensor, parameter=name, n_train=len(a), eky_n=len(b))); continue
        m = RandomForestRegressor(500, **RF).fit(Xa.values, np.log10(a[t].values))
        y = np.log10(b[t].values); p = m.predict(Xb.values)
        B.append(dict(sensor=sensor, parameter=name, n_train=len(a), eky_n=len(y), reservoirs=b.lake.nunique(),
                      rho=round(spearmanr(y, p)[0], 2), MAPE=round(mape(y, p), 1), bias=round(bias(y, p), 2)))
B = pd.DataFrame(B); B.to_csv(f"{OUT}/hls_eky_test.csv", index=False)
print(B.to_string(index=False))

# ---------------- (C) same-day reflectance agreement ----------------
orig = pd.read_csv(f"{R}/training_paired_s2_landsat.csv")
o = orig[orig.img_gap_days == 0]
C = []
for band, s2b in [("blue", "B2"), ("green", "B3"), ("red", "B4"), ("nir", "B8A")]:
    r_o = (o["L_" + band] / o[s2b]); r_o = r_o[(o["L_" + band] > 0) & (o[s2b] > 0)]
    sd = pair[pair.gap == 0]; r_h = sd[band + "_L"] / sd[band]; r_h = r_h[(sd[band + "_L"] > 0) & (sd[band] > 0)]
    C.append(dict(band=band, orig_n=len(r_o), orig_median_L_over_S2=round(r_o.median(), 3),
                  orig_median_abs_diff_pct=round(100 * np.median(np.abs(r_o - 1)), 1),
                  hls_n=len(r_h), hls_median_L30_over_S30=round(r_h.median(), 3),
                  hls_median_abs_diff_pct=round(100 * np.median(np.abs(r_h - 1)), 1)))
C = pd.DataFrame(C); C.to_csv(f"{OUT}/hls_sameday_reflectance.csv", index=False)
print(C.to_string(index=False))
