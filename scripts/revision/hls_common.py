import numpy as np, pandas as pd
H=R=None
LAB=["turb_ntu","secchi_m","tss_mgl","chla_ugl","fe_mgl"]
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

