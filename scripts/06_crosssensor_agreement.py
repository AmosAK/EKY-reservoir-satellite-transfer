"""Same-day agreement between Sentinel-2 and Landsat 8/9 retrievals at EKY stations,
2017-2026 (Table 9). Only pairs in which both retrievals are in-domain are used.
Output: results/crosssensor_agreement.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import DATA, RESULTS, TARGETS

s2 = pd.read_csv(DATA / "eky_timeseries_sentinel2.csv")
l8 = pd.read_csv(DATA / "eky_timeseries_landsat89.csv")
for d in (s2, l8):
    d["img_date"] = pd.to_datetime(d.img_date).dt.normalize()
pairs = s2.merge(l8, on=["station", "img_date"], suffixes=("_s2", "_l"))

rows = []
for short, name in zip(["turb", "secchi", "tss", "chla"], TARGETS.values()):
    k = pairs[pairs[f"{short}_in_domain_s2"] & pairs[f"{short}_in_domain_l"]]
    r = np.log10(k[f"{short}_l"] / k[f"{short}_s2"])
    rows.append({"parameter": name, "same_day_pairs": len(k),
                 "rho": spearmanr(k[f"{short}_s2"], k[f"{short}_l"])[0],
                 "median_landsat_over_s2": 10 ** np.median(r),
                 "median_abs_difference_pct": 100 * np.median(np.abs(10 ** r - 1))})

out = pd.DataFrame(rows).round(3)
out.to_csv(RESULTS / "crosssensor_agreement.csv", index=False)
print(f"same-day station-dates: {len(pairs)}")
print(out.to_string(index=False))
