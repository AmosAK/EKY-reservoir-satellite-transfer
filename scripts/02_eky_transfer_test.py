"""Independent test on the withheld Eastern Kentucky (EKY) reservoirs.

Trains the final random forests on the full regional training set, predicts the
withheld EKY matchups and reports (Table 2, Table 5 and Table S2):
  * Spearman rho, MAPE and bias ratio
  * 95% bootstrap confidence intervals, resampling matchups and whole reservoirs
  * rho with each reservoir left out, and rho within reservoirs
Outputs:
    results/eky_test_predictions.csv
    results/eky_test_metrics.csv
    results/eky_test_bootstrap_ci.csv
    results/eky_test_leave_one_reservoir_out.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor

from common import (FEATURES, N_TREES_FINAL, RESULTS, RF_PARAMS, SENSORS, TARGETS,
                    load_eky_test, load_training, metrics, model_table)

N_BOOT = 5000
rng = np.random.default_rng(42)

# ---------------------------------------------------------------- predictions
pred_rows = []
for sensor in SENSORS:
    train, test = load_training(sensor), load_eky_test(sensor)
    Xt = FEATURES[sensor](test)
    ok = Xt.notna().all(axis=1)
    for target in TARGETS:
        _, X, y = model_table(train, target, FEATURES[sensor])
        rf = RandomForestRegressor(N_TREES_FINAL, **RF_PARAMS).fit(X, y)
        t = test[ok].copy()
        t["predicted"] = 10 ** rf.predict(Xt[ok])
        t = t[t[target] > 0]
        pred_rows.append(pd.DataFrame({"sample_id": t.sample_id, "sensor": sensor, "lake": t.lake,
                                       "target": target, "observed": t[target], "predicted": t.predicted}))
P = pd.concat(pred_rows)
P.to_csv(RESULTS / "eky_test_predictions.csv", index=False)


# ---------------------------------------------------------------- statistics
def basic(y, p):
    r = p - y
    return spearmanr(y, p)[0], 100 * np.median(np.abs(10 ** r - 1)), np.median(10 ** r)


def q(a, k):
    return np.nanpercentile(a[:, k], 2.5), np.nanpercentile(a[:, k], 97.5)


met_rows, ci_rows, loo_rows = [], [], []
for (sensor, target), g in P.groupby(["sensor", "target"], sort=False):
    y, p, lakes = np.log10(g.observed.values), np.log10(g.predicted.values), g.lake.values
    m = metrics(y, p)
    met_rows.append({"sensor": SENSORS[sensor], "parameter": TARGETS[target], **m,
                     "reservoirs": len(np.unique(lakes))})

    by_sample, by_lake = [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(y), len(y))
        if np.ptp(y[i]) > 0:
            by_sample.append(basic(y[i], p[i]))
    uniq = np.unique(lakes)
    for _ in range(N_BOOT):
        pick = rng.choice(uniq, len(uniq))
        i = np.concatenate([np.where(lakes == c)[0] for c in pick])
        if np.ptp(y[i]) > 0 and np.ptp(p[i]) > 0:
            by_lake.append(basic(y[i], p[i]))
    bs, bl = np.array(by_sample), np.array(by_lake)
    ci_rows.append({"sensor": SENSORS[sensor], "parameter": TARGETS[target], "n": len(y),
                    "reservoirs": len(uniq), "rho": m["rho"],
                    "rho_ci_matchups": "%.2f to %.2f" % q(bs, 0),
                    "rho_ci_reservoirs": "%.2f to %.2f" % q(bl, 0),
                    "MAPE": m["MAPE"], "MAPE_ci_matchups": "%.0f to %.0f" % q(bs, 1),
                    "bias_ratio": m["bias_ratio"], "bias_ci_matchups": "%.2f to %.2f" % q(bs, 2)})

    loo = {lk: spearmanr(y[lakes != lk], p[lakes != lk])[0] for lk in uniq}
    yo = y - pd.Series(y).groupby(lakes).transform("mean").values
    po = p - pd.Series(p).groupby(lakes).transform("mean").values
    loo_rows.append({"sensor": SENSORS[sensor], "parameter": TARGETS[target], "n": len(y),
                     "reservoirs": len(uniq), "rho_all": m["rho"],
                     "rho_loo_min": min(loo.values()), "rho_loo_max": max(loo.values()),
                     "most_influential_reservoir": max(loo, key=lambda k: abs(loo[k] - m["rho"])),
                     "rho_within_reservoirs": spearmanr(yo, po)[0]})

pd.DataFrame(met_rows).round(3).to_csv(RESULTS / "eky_test_metrics.csv", index=False)
pd.DataFrame(ci_rows).round(3).to_csv(RESULTS / "eky_test_bootstrap_ci.csv", index=False)
pd.DataFrame(loo_rows).round(3).to_csv(RESULTS / "eky_test_leave_one_reservoir_out.csv", index=False)
print(pd.DataFrame(met_rows).round(2).to_string(index=False))
