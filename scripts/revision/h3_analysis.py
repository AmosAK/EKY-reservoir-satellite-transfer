"""Exact sampling-to-overpass gaps (hours) and the paired window test.

Calendar offset of a matchup = floor(dt_days)  (dt_days = image time minus 00:00 UTC of the sampling date).
Overpass time: from HLS exports (station, img_date), else station median, else sensor median.
Sampling time: WQP ActivityStartTime converted to UTC; median over activities at a station-date.
"""
import sys, glob
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold

R, H, TIMES, EKYRAW, OUT = sys.argv[1:6]
RF = dict(min_samples_leaf=5, max_features=0.5, n_jobs=-1, random_state=1)
T = {"turb_ntu": "Turbidity", "secchi_m": "Secchi depth", "tss_mgl": "TSS", "chla_ugl": "Chlorophyll-a"}
TZ = {"EST": -5, "EDT": -4, "CST": -6, "CDT": -5, "MST": -7, "MDT": -6, "PST": -8, "PDT": -7, "AST": -4, "ADT": -3,
      "UTC": 0, "GMT": 0, "Z": 0}

def sampling_times():
    t = pd.read_csv(TIMES, dtype=str)
    t = t[t.get("ActivityMediaName", "Water").fillna("Water").str.contains("Water")]
    e = pd.read_csv(EKYRAW, dtype=str, usecols=lambda c: c in ["MonitoringLocationIdentifier", "ActivityStartDate",
                    "ActivityStartTime/Time", "ActivityStartTime/TimeZoneCode"])
    e = e.rename(columns={"MonitoringLocationIdentifier": "station", "ActivityStartDate": "date",
                          "ActivityStartTime/Time": "Time", "ActivityStartTime/TimeZoneCode": "TimeZoneCode"})
    t = pd.concat([t[["station", "date", "Time", "TimeZoneCode"]], e], ignore_index=True).dropna(subset=["Time"])
    off = t.TimeZoneCode.map(TZ)
    hh = pd.to_timedelta(t.Time, errors="coerce")
    t["samp_utc_h"] = hh.dt.total_seconds() / 3600 - off
    t = t.dropna(subset=["samp_utc_h"])
    return t.groupby(["station", "date"]).samp_utc_h.median().rename("samp_utc_h").reset_index()

def overpass_lookup(kind):
    rows = []
    for f in glob.glob(f"{H}/HLS_{kind}_*_matchups.csv"):
        d = pd.read_csv(f, usecols=["station", "img_date", "img_time_utc"])
        rows.append(d)
    d = pd.concat(rows).dropna()
    hm = d.img_time_utc.str.split(":", expand=True).astype(float)
    d["ov_h"] = hm[0] + hm[1] / 60
    return d.groupby(["station", "img_date"]).ov_h.median()

def add_gap(d, kind, ST):
    d = d.copy()
    d["cal"] = np.floor(d.dt_days).astype(int)
    if "img_date" not in d:
        d["img_date"] = (pd.to_datetime(d.date) + pd.to_timedelta(d.cal, "D")).dt.strftime("%Y-%m-%d")
    lk = overpass_lookup(kind)
    d = d.merge(lk.rename("ov_h").reset_index(), on=["station", "img_date"], how="left")
    st_med = d.groupby("station").ov_h.transform("median")
    d["ov_src"] = np.where(d.ov_h.notna(), "scene", np.where(st_med.notna(), "station", "regional"))
    d["ov_h"] = d.ov_h.fillna(st_med).fillna(lk.median())
    d = d.merge(ST, on=["station", "date"], how="left")
    d["gap_h"] = np.abs(d.cal * 24 + d.ov_h - d.samp_utc_h)
    return d

def r2(y, p): return 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)
def mape(y, p): return 100 * np.median(np.abs(10 ** (p - y) - 1))
def s2feat(d):
    X = d[["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A"]].copy()
    X["ndci"] = (d.B5 - d.B4) / (d.B5 + d.B4); X["b5b4"] = d.B5 / d.B4; X["b3b2"] = d.B3 / d.B2
    X["b4b3"] = d.B4 / d.B3; X["b2b4"] = d.B2 / d.B4; X["b8b4"] = d.B8 / d.B4
    X["ndti"] = (d.B4 - d.B3) / (d.B4 + d.B3); X["b6b5"] = d.B6 / d.B5
    return X.replace([np.inf, -np.inf], np.nan)
def lfeat(d):
    b, g, r, n = d.blue, d.green, d.red, d.nir
    X = pd.DataFrame({"blue": b, "green": g, "red": r, "nir": n, "br": b / r, "gr": g / r, "bg": b / g, "nr": n / r,
                      "ndti": (r - g) / (r + g)})
    return X.replace([np.inf, -np.inf], np.nan)

ST = sampling_times()
data = {"Sentinel-2": (add_gap(pd.read_csv(f"{R}/training_sentinel2.csv"), "S30", ST), s2feat),
        "Landsat 8/9": (add_gap(pd.read_csv(f"{R}/training_landsat89.csv"), "L30", ST), lfeat)}
eky = {"Sentinel-2": add_gap(pd.read_csv(f"{R}/eky_test_sentinel2.csv"), "S30", ST),
       "Landsat 8/9": add_gap(pd.read_csv(f"{R}/eky_test_landsat89.csv"), "L30", ST)}

# ---- coverage table
cov = []
for k, (d, _) in data.items():
    cov.append(dict(set="training " + k, n=len(d), with_sampling_time=int(d.samp_utc_h.notna().sum()),
                    overpass_from_scene=int((d.ov_src == "scene").sum()),
                    within_3h=int((d.gap_h <= 3).sum()), within_6h=int((d.gap_h <= 6).sum()),
                    same_calendar_day=int((d.cal == 0).sum()), within_1d=int((d.cal.abs() <= 1).sum()),
                    median_gap_h_sameday=round(d.loc[d.cal == 0, "gap_h"].median(), 2)))
for k, d in eky.items():
    cov.append(dict(set="EKY test " + k, n=len(d), with_sampling_time=int(d.samp_utc_h.notna().sum()),
                    within_3h=int((d.gap_h <= 3).sum()), within_6h=int((d.gap_h <= 6).sum()),
                    same_calendar_day=int((d.cal == 0).sum()), within_1d=int((d.cal.abs() <= 1).sum())))
cov = pd.DataFrame(cov); cov.to_csv(f"{OUT}/gap_coverage.csv", index=False); print(cov.to_string(index=False), flush=True)

# ---- paired window test: identical +/-3 h test samples at held-out stations
WIN = {"3h": lambda d: d.gap_h <= 3, "sameday": lambda d: d.cal == 0, "1d": lambda d: d.cal.abs() <= 1,
       "3d": lambda d: d.cal.abs() <= 3}
rows = []
for k, (d0, fz) in data.items():
    for t, name in T.items():
        d = d0[d0[t] > 0]; X = fz(d); ok = X.notna().all(axis=1).values; d, X = d[ok].reset_index(drop=True), X[ok].values
        y = np.log10(d[t].values); g = d.station.values
        test_mask = (d.gap_h <= 3).values
        masks = {w: f(d).values for w, f in WIN.items()}
        P = {w: np.full(len(y), np.nan) for w in WIN}
        for a, b in GroupKFold(5).split(X, y, g):
            bt = b[test_mask[b]]
            if len(bt) == 0: continue
            for w in WIN:
                aa = a[masks[w][a]]
                P[w][bt] = RandomForestRegressor(300, **RF).fit(X[aa], y[aa]).predict(X[bt])
        m = ~np.isnan(P["3h"]); yt, gt = y[m], g[m]
        E = {w: pd.Series(np.abs(p[m] - yt)).groupby(gt).mean() for w, p in P.items()}
        row = dict(sensor=k, parameter=name, test_n=int(m.sum()), test_stations=len(set(gt)),
                   train_n_3h=int(masks["3h"].sum()), train_n_3d=int(masks["3d"].sum()))
        for w, p in P.items():
            row[f"R2_{w}"] = round(r2(yt, p[m]), 2); row[f"MAPE_{w}"] = round(mape(yt, p[m]), 1)
            row[f"rho_{w}"] = round(spearmanr(yt, p[m])[0], 2)
        row["p_3d_vs_3h"] = round(wilcoxon(E["3d"], E["3h"]).pvalue, 3)
        row["p_3d_vs_sameday"] = round(wilcoxon(E["3d"], E["sameday"]).pvalue, 3)
        rows.append(row); print(row, flush=True)
W = pd.DataFrame(rows); W.to_csv(f"{OUT}/window_paired_3h.csv", index=False)

# ---- EKY test: models trained on each window, scored on the full EKY test set
rows = []
for k, (d0, fz) in data.items():
    te = eky[k]
    for t, name in T.items():
        b = te[te[t] > 0]; Xb = fz(b); okb = Xb.notna().all(axis=1).values; b, Xb = b[okb], Xb[okb].values
        yb = np.log10(b[t].values)
        row = dict(sensor=k, parameter=name, eky_n=len(yb))
        for w, f in WIN.items():
            a = d0[(d0[t] > 0) & f(d0)]; Xa = fz(a); ok = Xa.notna().all(axis=1).values
            m = RandomForestRegressor(500, **RF).fit(Xa[ok].values, np.log10(a[t].values[ok]))
            row[f"rho_{w}"] = round(spearmanr(yb, m.predict(Xb))[0], 2); row[f"ntrain_{w}"] = int(ok.sum())
        rows.append(row); print(row, flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/eky_by_window.csv", index=False)

# ---- error vs gap (for a supplementary figure): CV errors from the +/-3 d models binned by gap
rows = []
for k, (d0, fz) in data.items():
    for t, name in T.items():
        d = d0[(d0[t] > 0)]; X = fz(d); ok = X.notna().all(axis=1).values; d, X = d[ok].reset_index(drop=True), X[ok].values
        y = np.log10(d[t].values); g = d.station.values; p = np.zeros(len(y))
        for a, b in GroupKFold(5).split(X, y, g):
            p[b] = RandomForestRegressor(300, **RF).fit(X[a], y[a]).predict(X[b])
        d["abs_log_err"] = np.abs(p - y)
        bins = pd.cut(d.gap_h, [0, 3, 24, 48, 72, 100], include_lowest=True)
        s = d.groupby(bins, observed=True).abs_log_err.agg(["median", "count"])
        for iv, r_ in s.iterrows():
            rows.append(dict(sensor=k, parameter=name, gap_bin_h=str(iv), median_abs_log_err=round(r_["median"], 3),
                             n=int(r_["count"])))
        rho = spearmanr(d.gap_h, d.abs_log_err, nan_policy="omit")[0]
        rows.append(dict(sensor=k, parameter=name, gap_bin_h="spearman(gap, error)", median_abs_log_err=round(rho, 3), n=int(d.gap_h.notna().sum())))
E = pd.DataFrame(rows); E.to_csv(f"{OUT}/error_vs_gap_v2.csv", index=False); print(E.to_string(index=False))
