"""Robustness checks on the withheld EKY test (Section 4.4, Tables 7 and S4).

  * scene cloud limits (<80%, <50%, <20%)
  * training weighted to give each state equal influence
  * training restricted to Appalachian lakes (35-39.5 N, 79-85.5 W)
  * training without matchups inside the Eastern Kentucky Coal Field (optional;
    requires the Kentucky Geological Survey physiographic regions shapefile,
    https://www.uky.edu/KGS/gis/regions.htm, passed with --coalfield PATH)
Output: results/robustness_checks.csv
"""
import argparse

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

from common import (FEATURES, N_TREES_FINAL, RESULTS, RF_PARAMS, SENSORS, TARGETS,
                    load_eky_test, load_training, metrics, model_table)

ap = argparse.ArgumentParser()
ap.add_argument("--coalfield", help="KGS physiographic regions shapefile (optional)")
args = ap.parse_args()

coalfield = None
if args.coalfield:
    import geopandas as gpd
    reg = gpd.read_file(args.coalfield)
    coalfield = reg[reg.REGION == "Eastern Coal Field"].to_crs(3089).union_all()


def distance_to_coalfield_km(d):
    import geopandas as gpd
    pts = gpd.GeoSeries(gpd.points_from_xy(d.longitude, d.latitude), crs=4326).to_crs(3089)
    return pts.distance(coalfield).values * 0.3048 / 1000  # US survey feet to km


rows = []
for sensor, sensor_name in SENSORS.items():
    train0, test = load_training(sensor), load_eky_test(sensor)
    Xt_all = FEATURES[sensor](test)
    ok_t = Xt_all.notna().all(axis=1)
    train0["appalachian"] = train0.latitude.between(35, 39.5) & train0.longitude.between(-85.5, -79)
    if coalfield is not None:
        train0["dist_coalfield_km"] = distance_to_coalfield_km(train0)

    variants = {"as published": train0,
                "scene cloud < 80%": train0[train0.cloud_pct < 80],
                "scene cloud < 50%": train0[train0.cloud_pct < 50],
                "scene cloud < 20%": train0[train0.cloud_pct < 20],
                "Appalachian lakes only": train0[train0.appalachian],
                "state-balanced weights": train0}
    if coalfield is not None:
        variants["coalfield lakes removed"] = train0[train0.dist_coalfield_km > 0]
        variants["coalfield + 10 km removed"] = train0[train0.dist_coalfield_km > 10]

    for target, target_name in TARGETS.items():
        for variant, tr in variants.items():
            a, X, y = model_table(tr, target, FEATURES[sensor])
            w = None
            if variant == "state-balanced weights":
                st = a.state.astype(str)
                w = (1 / st.map(st.value_counts())).values
                w = w / w.mean()
            rf = RandomForestRegressor(N_TREES_FINAL, **RF_PARAMS).fit(X, y, sample_weight=w)
            t = test[ok_t & (test[target] > 0)]
            m = metrics(np.log10(t[target].values), rf.predict(Xt_all.loc[t.index]))
            rows.append({"sensor": sensor_name, "parameter": target_name, "variant": variant,
                         "n_train": len(y), "eky_n": m["n"], "eky_rho": m["rho"],
                         "eky_MAPE": m["MAPE"], "eky_bias_ratio": m["bias_ratio"]})
            print(rows[-1])

pd.DataFrame(rows).round(3).to_csv(RESULTS / "robustness_checks.csv", index=False)
