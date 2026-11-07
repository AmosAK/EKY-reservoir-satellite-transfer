"""Permutation importance of model inputs on station-held-out folds (Figure 6).

Uses the fixed station folds stored in the 'fold' column of the training files
(five shuffles per fold, drop in R2 when an input is shuffled).
Output: results/permutation_importance.csv
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance

from common import (FEATURES, N_TREES_CV, RESULTS, RF_PARAMS, SENSORS, TARGETS,
                    load_training, model_table)

rows = []
for sensor, sensor_name in SENSORS.items():
    d = load_training(sensor)
    for target, target_name in TARGETS.items():
        a, X, y = model_table(d, target, FEATURES[sensor])
        per_fold = []
        for k in range(5):
            tr, te = a.fold.values != k, a.fold.values == k
            rf = RandomForestRegressor(N_TREES_CV, **RF_PARAMS).fit(X[tr], y[tr])
            r = permutation_importance(rf, X[te], y[te], n_repeats=5, random_state=k, scoring="r2", n_jobs=-1)
            per_fold.append(r.importances_mean)
        per_fold = np.array(per_fold)
        for j, col in enumerate(X.columns):
            rows.append({"sensor": sensor_name, "parameter": target_name, "feature": col,
                         "importance": per_fold[:, j].mean(), "sd": per_fold[:, j].std()})
        print(sensor_name, target_name, "done")

pd.DataFrame(rows).to_csv(RESULTS / "permutation_importance.csv", index=False)
