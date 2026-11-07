"""Figures 4-8 from the outputs of scripts 01-06 (600 dpi PNG in figures/)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
from scipy.stats import spearmanr

from common import DATA, FIGURES, RESULTS, TARGETS

plt.rcParams.update({"font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 9,
                     "axes.linewidth": 0.6, "axes.edgecolor": "#555", "xtick.color": "#333",
                     "ytick.color": "#333", "axes.labelcolor": "#222"})
BLUE, ORANGE, INK = "#2a78d6", "#eb6834", "#222222"
PARAMS = list(TARGETS.values())
plain = FuncFormatter(lambda v, _: f"{v:g}")


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=600, facecolor="white")
    plt.close(fig)


def log_axes(ax, lo, hi, ticks):
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(FixedLocator(ticks)); axis.set_minor_locator(NullLocator())
        axis.set_major_formatter(plain)
    ax.spines[["top", "right"]].set_visible(False)


# ------------------------------------------------------------ Figure 4: sampling dates
v = pd.read_csv(DATA / "eky_field_visits_screened.csv", parse_dates=["date"]).drop_duplicates(["lake", "date"])
order = ["Smoky Valley", "Grayson", "Yatesville", "Cave Run", "Paintsville", "Dewey", "Fishtrap",
         "Buckhorn", "Cranks Creek", "Martins Fork", "Cannon Creek", "Chenoa"]
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for i, lake in enumerate(order):
    y = len(order) - 1 - i
    d = v[v.lake == lake]
    ax.axhline(y, color="#eeeeee", lw=0.6, zorder=0)
    if len(d):
        ax.scatter(d.date, [y] * len(d), s=22, color=BLUE, edgecolor="white", linewidth=0.6, zorder=3)
    else:
        ax.text(pd.Timestamp("2019-07-01"), y, "no field data 2015–2023", va="center", ha="center",
                fontsize=7.5, color="#777", style="italic")
ax.axvline(pd.Timestamp("2017-03-28"), color="#555", lw=0.8, ls="--")
ax.text(pd.Timestamp("2017-04-20"), len(order) - 0.45, "Sentinel-2 L2A archive begins", fontsize=7.5, color="#444")
ax.set_yticks(range(len(order)), order[::-1]); ax.set_ylim(-0.6, len(order) - 0.1)
ax.set_xlim(pd.Timestamp("2015-01-01"), pd.Timestamp("2024-01-31"))
ax.xaxis.set_major_locator(mdates.YearLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax.set_xlabel("Sampling date"); ax.spines[["top", "right"]].set_visible(False); ax.tick_params(axis="y", length=0)
fig.tight_layout(); save(fig, "Figure4_sampling_dates")

# ------------------------------------------------------------ Figure 5: performance
cv = pd.read_csv(RESULTS / "regional_cv.csv")
ci = pd.read_csv(RESULTS / "eky_test_bootstrap_ci.csv")
x = np.arange(4); w = 0.36
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0))
for i, (sensor, col) in enumerate([("Sentinel-2", BLUE), ("Landsat 8/9", ORANGE)]):
    r2 = cv[cv.sensor == sensor].set_index("parameter").loc[PARAMS, "R2"].values
    xs = x + (i - 0.5) * (w + 0.02)
    a.bar(xs, r2, w, color=col, label=sensor)
    for xx, vv in zip(xs, r2):
        a.text(xx, vv + 0.015, f"{vv:.2f}", ha="center", fontsize=7, color="#333")
    c = ci[ci.sensor == sensor].set_index("parameter").loc[PARAMS]
    lo_hi = np.array([[float(s) for s in t.split(" to ")] for t in c.rho_ci_reservoirs])
    err = np.array([c.rho - lo_hi[:, 0], lo_hi[:, 1] - c.rho])
    b.errorbar(x + (i - 0.5) * 0.26, c.rho, yerr=err, fmt="o" if i == 0 else "s", ms=5, color=col,
               mec="white", mew=0.8, elinewidth=1, capsize=2.5, label=sensor)
n = {s: ci[ci.sensor == s].set_index("parameter").loc[PARAMS, "n"].values for s in ["Sentinel-2", "Landsat 8/9"]}
for ax in (a, b):
    ax.spines[["top", "right"]].set_visible(False); ax.yaxis.grid(True, color="#e6e6e6", lw=0.5)
    ax.set_axisbelow(True); ax.tick_params(axis="x", length=0)
a.set_xticks(x, PARAMS); a.set_ylim(0, 1); a.set_ylabel("Cross-validated R² (log₁₀ values)")
a.legend(frameon=False, fontsize=8)
b.set_xticks(x, [f"{p}\n(n = {s}/{l})" for p, s, l in zip(PARAMS, n["Sentinel-2"], n["Landsat 8/9"])])
b.axhline(0, color="#777", lw=0.6); b.set_ylim(-1.05, 1.05); b.set_ylabel("EKY test Spearman ρ")
b.legend(frameon=False, fontsize=8, loc="lower left")
for ax, letter in ((a, "a"), (b, "b")):
    ax.text(-0.16, 1.02, f"({letter})", transform=ax.transAxes, fontsize=10, fontweight="bold")
fig.tight_layout(w_pad=2.5); save(fig, "Figure5_performance")

# ------------------------------------------------------------ Figure 6: importance
imp = pd.read_csv(RESULTS / "permutation_importance.csv")
names = {"B5": "B5 red-edge 705", "B8A": "B8A NIR narrow", "B4": "B4 red", "B3": "B3 green",
         "ndci": "NDCI (B5, B4)", "ndti": "NDTI (red, green)", "b5b4": "B5/B4", "b6b5": "B6/B5",
         "b2b4": "B2/B4", "b4b3": "B4/B3", "b3b2": "B3/B2", "b8b4": "B8/B4", "blue": "Blue",
         "green": "Green", "red": "Red", "nir": "NIR", "br": "Blue/red", "gr": "Green/red",
         "bg": "Blue/green", "nr": "NIR/red"}
red_edge = lambda f: any(k in f.lower() for k in ["b5", "b6", "b7", "ndci"])
fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.6))
for r, sensor in enumerate(["Sentinel-2", "Landsat 8/9"]):
    for c, p in enumerate(PARAMS):
        ax = axes[r, c]
        d = imp[(imp.sensor == sensor) & (imp.parameter == p)].nlargest(6, "importance")[::-1]
        cols = [ORANGE if sensor == "Sentinel-2" and red_edge(f) else BLUE for f in d.feature]
        ax.barh(range(len(d)), d.importance, xerr=d["sd"], color=cols, height=0.7,
                error_kw=dict(ecolor="#555", elinewidth=0.6))
        ax.set_yticks(range(len(d)), [names.get(f, f) for f in d.feature], fontsize=7)
        ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False); ax.set_title(p, fontsize=8.5, loc="left")
        if r == 1:
            ax.set_xlabel("Drop in R²", fontsize=7.5)
h = [plt.Rectangle((0, 0), 1, 1, color=ORANGE), plt.Rectangle((0, 0), 1, 1, color=BLUE)]
fig.legend(h, ["Uses a Sentinel-2 red-edge band (B5–B7)", "Visible/NIR bands and ratios"],
           loc="lower center", ncol=2, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0.05, 1, 0.95), h_pad=3.2, w_pad=0.8)
for r, sensor in enumerate(["Sentinel-2", "Landsat 8/9"]):
    fig.text(0.01, axes[r, 0].get_position().y1 + 0.055, f"({'ab'[r]}) {sensor}", fontsize=9.5, fontweight="bold")
save(fig, "Figure6_permutation_importance")

# ------------------------------------------------------------ Figure 7: observed vs predicted
cvp = pd.read_csv(RESULTS / "regional_cv_predictions_S2.csv")
ek = pd.read_csv(RESULTS / "eky_test_predictions.csv")
blues = LinearSegmentedColormap.from_list("b", ["#cde2fb", "#86b6ef", "#2a78d6", "#184f95", "#0d366b"])
spec = {"turb_ntu": ("Turbidity (NTU)", (0.8, 100), [1, 2, 5, 10, 20, 50, 100]),
        "secchi_m": ("Secchi depth (m)", (0.2, 5), [0.2, 0.5, 1, 2, 5]),
        "tss_mgl": ("TSS (mg L$^{-1}$)", (1, 60), [1, 2, 5, 10, 20, 50]),
        "chla_ugl": ("Chlorophyll-a (µg L$^{-1}$)", (0.7, 200), [1, 2, 5, 10, 20, 50, 100, 200])}
fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.0))
for ax, (t, (label, (lo, hi), ticks)), letter in zip(axes.flat, spec.items(), "abcd"):
    c = cvp[cvp.target == t]
    e = ek[(ek.target == t) & (ek.sensor == "S2")]
    y, p = np.log10(c.observed), np.log10(c.predicted)
    r2 = 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)
    rho = spearmanr(e.observed, e.predicted)[0]
    ax.hexbin(c.observed, c.predicted, xscale="log", yscale="log", gridsize=38, mincnt=1, bins="log",
              cmap=blues, linewidths=0, extent=(np.log10(lo), np.log10(hi), np.log10(lo), np.log10(hi)))
    ax.plot([lo, hi], [lo, hi], ls="--", color="#333", lw=0.8)
    ax.scatter(e.observed, e.predicted, s=26, color=ORANGE, edgecolor="white", linewidth=0.8, zorder=3)
    log_axes(ax, lo, hi, ticks)
    ax.set_xlabel(f"Observed {label}"); ax.set_ylabel(f"Predicted {label}")
    ax.text(0.04, 0.96, f"Regional CV: R² = {r2:.2f}, n = {len(c):,}\nEKY test: ρ = {rho:.2f}, n = {len(e)}",
            transform=ax.transAxes, va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ccc", lw=0.5))
    ax.text(-0.2, 1.03, f"({letter})", transform=ax.transAxes, fontsize=10, fontweight="bold")
h = [plt.Line2D([], [], marker="h", ls="", color=BLUE, ms=8, label="Regional matchups (density)"),
     plt.Line2D([], [], marker="o", ls="", color=ORANGE, mec="white", ms=6, label="Withheld EKY reservoirs"),
     plt.Line2D([], [], ls="--", color="#333", lw=0.8, label="1:1")]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
fig.tight_layout(rect=(0, 0.04, 1, 1), h_pad=1.5, w_pad=2); save(fig, "Figure7_observed_vs_predicted")

# ------------------------------------------------------------ Figure 8: same-day agreement
s2 = pd.read_csv(DATA / "eky_timeseries_sentinel2.csv"); l8 = pd.read_csv(DATA / "eky_timeseries_landsat89.csv")
for d in (s2, l8):
    d["img_date"] = pd.to_datetime(d.img_date).dt.normalize()
j = s2.merge(l8, on=["station", "img_date"], suffixes=("_s2", "_l"))
spec8 = [("turb", "Turbidity (NTU)", (0.8, 100), [1, 2, 5, 10, 20, 50, 100]),
         ("secchi", "Secchi depth (m)", (0.2, 6), [0.2, 0.5, 1, 2, 5]),
         ("tss", "TSS (mg L$^{-1}$)", (1, 40), [1, 2, 5, 10, 20]),
         ("chla", "Chlorophyll-a (µg L$^{-1}$)", (0.7, 100), [1, 2, 5, 10, 20, 50, 100])]
fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.25))
for ax, (p, label, (lo, hi), ticks), letter in zip(axes, spec8, "abcd"):
    k = j[j[f"{p}_in_domain_s2"] & j[f"{p}_in_domain_l"]]
    rho = spearmanr(k[f"{p}_s2"], k[f"{p}_l"])[0]
    ratio = 10 ** np.median(np.log10(k[f"{p}_l"] / k[f"{p}_s2"]))
    ax.scatter(k[f"{p}_s2"], k[f"{p}_l"], s=5, color=BLUE, alpha=0.45, linewidth=0)
    ax.plot([lo, hi], [lo, hi], ls="--", color="#333", lw=0.7)
    log_axes(ax, lo, hi, ticks); ax.tick_params(labelsize=7)
    ax.set_xlabel("Sentinel-2", fontsize=8)
    if letter == "a":
        ax.set_ylabel("Landsat 8/9", fontsize=8)
    ax.set_title(f"({letter}) {label}", fontsize=8, loc="left")
    ax.text(0.97, 0.04, f"n = {len(k)}\nρ = {rho:.2f}\nL/S2 = {ratio:.2f}", transform=ax.transAxes,
            va="bottom", ha="right", fontsize=6.8)
fig.tight_layout(w_pad=1.0); save(fig, "Figure8_same_day_agreement")
print("figures written to", FIGURES)
