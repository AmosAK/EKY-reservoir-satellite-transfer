"""ACOLITE on the EKY test scenes: download -> ACOLITE (small box per scene) -> 3x3 station medians -> clean up."""
import sys, os, glob, shutil, subprocess, json
import numpy as np, pandas as pd
from netCDF4 import Dataset

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# usage: python scripts/revision/acolite_pipeline.py WORKDIR SENSOR(S2|L) [ACOLITE_DIR]
A = os.path.abspath(sys.argv[1]); SENSOR = sys.argv[2]
ACO = sys.argv[3] if len(sys.argv) > 3 else f"{A}/acolite"   # clone of https://github.com/acolite/acolite
HERE = os.path.dirname(os.path.abspath(__file__)); WORK = f"{A}/work"; OUTNC = f"{A}/nc"; os.makedirs(OUTNC, exist_ok=True)
lst = pd.read_csv(f"{REPO}/data/acolite/eky_test_scene_list.csv")
st = pd.read_csv(f"{REPO}/data/eky_stations.csv").set_index("station")
lst = lst[lst.sensor == SENSOR]
prods = {}
for r in lst.itertuples():
    for p in r.product_ids.split(";"):
        prods.setdefault(p, set()).add((r.station, r.img_date))
done_f = f"{A}/extract_{SENSOR}.csv"
done = set(pd.read_csv(done_f).product) if os.path.exists(done_f) else set()
HALF = 1 if SENSOR == "L" else 1                     # 3x3 pixels either way

def window_median(nc, lat0, lon0):
    lat, lon = nc["lat"][:], nc["lon"][:]
    d = (lat - lat0) ** 2 + ((lon - lon0) * np.cos(np.radians(lat0))) ** 2
    i, j = np.unravel_index(np.argmin(d), d.shape)
    if np.sqrt(d[i, j]) * 111000 > 30: return None
    out = {}
    for v in nc.variables:
        if v.startswith(("rhos_", "rhow_")):
            a = np.ma.filled(nc[v][max(i - HALF, 0):i + HALF + 1, max(j - HALF, 0):j + HALF + 1].astype(float), np.nan)
            out[v] = a
    return out

for p, sts in sorted(prods.items()):
    if p in done: continue
    lats = [st.loc[s, "latitude"] for s, _ in sts]; lons = [st.loc[s, "longitude"] for s, _ in sts]
    box = [min(lats) - 0.05, min(lons) - 0.05, max(lats) + 0.05, max(lons) + 0.05]
    shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
    if SENSOR == "S2":
        r = subprocess.run([sys.executable, f"{HERE}/acolite_get_s2.py", p, WORK], capture_output=True, text=True)
        inp = f"{WORK}/{p}.SAFE"
    else:
        r = subprocess.run([sys.executable, f"{HERE}/acolite_get_landsat.py", p, WORK], capture_output=True, text=True)
        inp = f"{WORK}/{p}"
    if r.returncode: print("DOWNLOAD FAIL", p, r.stdout[-300:], r.stderr[-300:], flush=True); continue
    od = f"{WORK}/out"
    cfg = f"""inputfile={inp}
output={od}
limit={','.join(f'{x:.4f}' for x in box)}
s2_target_res=10
ancillary_data=False
dsf_aot_estimation=fixed
dsf_residual_glint_correction=True
l2w_parameters=rhow_*
l2r_export_geotiff=False
l1r_delete_netcdf=True
rgb_rhot=False
rgb_rhos=False
map_l2w=False
"""
    open(f"{WORK}/run.cfg", "w").write(cfg)
    r = subprocess.run([sys.executable, f"{ACO}/launch_acolite.py", "--cli", f"--settings={WORK}/run.cfg"],
                       capture_output=True, text=True)
    l2r = glob.glob(f"{od}/*_L2R.nc"); l2w = glob.glob(f"{od}/*_L2W.nc")
    if not l2r: print("ACOLITE FAIL", p, r.stdout[-500:], r.stderr[-500:], flush=True); continue
    rows = []
    with Dataset(l2r[0]) as a, (Dataset(l2w[0]) if l2w else Dataset(l2r[0])) as w:
        aot = getattr(a, "ac_aot_550", np.nan); model = getattr(a, "ac_model", "")
        for s, dte in sorted(sts):
            ra = window_median(a, st.loc[s, "latitude"], st.loc[s, "longitude"])
            if ra is None: continue
            rw = window_median(w, st.loc[s, "latitude"], st.loc[s, "longitude"]) if l2w else {}
            row = dict(product=p, station=s, img_date=dte, aot550=aot, model=model)
            # water mask from ACOLITE rhos (MNDWI > 0) as in the GEE extraction
            g = [k for k in ra if k.startswith("rhos_") and 540 <= int(k.split("_")[1]) <= 570]
            sw = [k for k in ra if k.startswith("rhos_") and 1550 <= int(k.split("_")[1]) <= 1660]
            mask = np.ones_like(next(iter(ra.values())), bool)
            if g and sw: mask = (ra[g[0]] - ra[sw[0]]) / (ra[g[0]] + ra[sw[0]]) > 0
            row["n_px"] = int(mask.sum())
            for k, v in {**ra, **(rw or {})}.items():
                vv = v[mask]; row[k] = float(np.nanmedian(vv)) if np.isfinite(vv).any() else np.nan
            rows.append(row)
    for f in l2r + l2w: shutil.copy(f, OUTNC)
    df = pd.DataFrame(rows) if rows else pd.DataFrame([dict(product=p)])
    df.to_csv(done_f, mode="a", header=not os.path.exists(done_f), index=False)
    print("DONE", p, len(rows), flush=True)
shutil.rmtree(WORK, ignore_errors=True)
print("ALL DONE")
