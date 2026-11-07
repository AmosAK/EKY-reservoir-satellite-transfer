"""Regional station-grouped five-fold cross-validation (Table 2, regional columns).

Also repeats the cross-validation with folds grouped by lake instead of station.
Outputs:
    results/regional_cv.csv
    results/regional_cv_predictions_S2.csv, results/regional_cv_predictions_L89.csv
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold

from common import (FEATURES, N_TREES_CV, RESULTS, RF_PARAMS, SENSORS, TARGETS,
                    load_training, metrics, model_table)


def grouped_cv(X, y, groups):
    pred = np.zeros(len(y))
    for train, test in GroupKFold(5).split(X, y, groups):
        rf = RandomForestRegressor(N_TREES_CV, **RF_PARAMS).fit(X.iloc[train], y[train])
        pred[test] = rf.predict(X.iloc[test])
    return pred


rows = []
for sensor, sensor_name in SENSORS.items():
    d = load_training(sensor)
    preds = []
    for target, target_name in TARGETS.items():
        a, X, y = model_table(d, target, FEATURES[sensor])
        p_station = grouped_cv(X, y, a.station)
        lake_key = a.lake.fillna(a.station).str.lower().str.replace(r"\s+", " ", regex=True).str.strip()
        p_lake = grouped_cv(X, y, lake_key)
        m, ml = metrics(y, p_station), metrics(y, p_lake)
        rows.append({"sensor": sensor_name, "parameter": target_name, "n_train": len(y),
                     "stations": a.station.nunique(), "R2": m["R2"], "MAPE": m["MAPE"],
                     "rho": m["rho"], "R2_lake_grouped": ml["R2"]})
        preds.append(pd.DataFrame({"sample_id": a.sample_id, "target": target,
                                   "observed": 10 ** y, "predicted": 10 ** p_station}))
        print(rows[-1])
    pd.concat(preds).to_csv(RESULTS / f"regional_cv_predictions_{sensor}.csv", index=False)

pd.DataFrame(rows).round(3).to_csv(RESULTS / "regional_cv.csv", index=False)
