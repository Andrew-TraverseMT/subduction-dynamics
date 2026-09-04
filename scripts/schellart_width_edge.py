#!/usr/bin/env python3
"""Schellart trench-zone width / edge-distance test (Cenozoic, zone–time blocked).

Reads the existing UPS-proxy catalog. Does **not** re-extract trench kinematics.
Does **not** use Müller v_T vs v_OP as H1/H2 (they are collinear in this model).

Inference unit: (time_Ma, zone_id), not tessellated points.

Two UPS outcomes, kept separate
  A  bab_zone  = 1 if any sample in the block is bab_on (MOR back-arc)
  B  ext_zone  = 1 if median dilatation among dilat_in_network==1 samples is > 0;
                 0 if the block is all rigid / zero dilatation (no network samples,
                 or network median ≤ 0). Hellenic-style extension is B, not A.

Schellart et al. (2007) prediction under test: narrower zones and samples nearer
a slab/trench edge should be more extensional. Wide-zone interiors more
compressive / bab_off. ~1200 km is the near-edge length scale they associated
with retreat. This test is reconstruction-internal (same Müller2019 topologies
that encode back-arc ridges). Do not treat a significant ρ as independent geology.

Run:
    /workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
        /workspace/subduction-test/scripts/schellart_width_edge.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, mannwhitneyu, spearmanr
from sklearn.linear_model import LogisticRegression

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import cenozoic_ups_proxy as cup  # noqa: E402
import submap_m56_h1h2 as h12  # noqa: E402

ROOT = Path("/workspace/subduction-test")
RESULTS = ROOT / "results"
PROXY_CSV = RESULTS / "trench_kinematics_0-60Ma_with_ups_proxy.csv"
JOIN_CSV = RESULTS / "joined_submap_muller2019_0Ma.csv"
REPORT_MD = RESULTS / "schellart_width_edge.md"
FIG_WIDTH = RESULTS / "fig_width_vs_bab.png"
FIG_EDGE = RESULTS / "fig_edge_vs_bab.png"
FIG_DILAT = RESULTS / "fig_dilatation_vs_width_edge.png"
FIG_MAP = RESULTS / "fig_width_edge_map_0Ma.png"

EDGE_NEAR_KM = 1200.0
EXT_MEDIAN_THRESH = 0.0  # 10^{-15} s^{-1}; ext_zone = 1 iff median > this
CE_SCORE = h12.CE_SCORE
UPS3_ORDER = h12.UPS3_ORDER
UPS3_COLORS = h12.UPS3_COLORS


def _setup_style() -> None:
    cup._setup_style()
    plt.rcParams["axes.grid"] = False


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p <= 0:
        return "<1e-16"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_rho(r: float) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def fmt_num(x: float, digits: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def _finite_pair(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


def spearman_safe(x, y):
    x, y = _finite_pair(x, y)
    n = int(len(x))
    if n < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, np.nan, n
    r, p = spearmanr(x, y)
    return float(r), float(p), n


def logit_univariate(x, y):
    return cup.logit_univariate(x, y)


def logit_with_time(df, x_col, y_col, time_col="time_Ma"):
    """Logistic y ~ x + time dummies. LRT vs time-only."""
    d = df[[x_col, y_col, time_col]].dropna()
    if len(d) < 30 or d[y_col].nunique() < 2:
        return None
    y = d[y_col].to_numpy(dtype=float)
    x = d[x_col].to_numpy(dtype=float)
    times = d[time_col].to_numpy(dtype=float)
    uniq = np.sort(np.unique(times))
    if uniq.size < 2:
        uni = logit_univariate(x, y)
        if uni is None:
            return None
        return {
            "n": uni["n"],
            "n_on": uni["n_on"],
            "coef": uni["coef"],
            "odds_ratio_per_unit": uni["odds_ratio_per_unit"],
            "lrt_vs_time_only": uni["lrt"],
            "p_lrt": uni["p_lrt"],
        }
    dummies = np.column_stack([(times == t).astype(float) for t in uniq[1:]])
    X_full = np.column_stack([x, dummies])
    X_red = dummies
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = LogisticRegression(
            penalty=None, solver="lbfgs", max_iter=4000
        ).fit(X_full, y)
        red = LogisticRegression(
            penalty=None, solver="lbfgs", max_iter=4000
        ).fit(X_red, y)
    coef = float(full.coef_[0, 0])
    p_full = np.clip(full.predict_proba(X_full)[:, 1], 1e-12, 1 - 1e-12)
    p_red = np.clip(red.predict_proba(X_red)[:, 1], 1e-12, 1 - 1e-12)
    ll_full = float(np.sum(y * np.log(p_full) + (1 - y) * np.log(1 - p_full)))
    ll_red = float(np.sum(y * np.log(p_red) + (1 - y) * np.log(1 - p_red)))
    lrt = 2.0 * (ll_full - ll_red)
    p_lrt = float(chi2.sf(max(lrt, 0.0), 1))
    return {
        "n": int(len(y)),
        "n_on": int(y.sum()),
        "coef": coef,
        "odds_ratio_per_unit": float(np.exp(coef)),
        "lrt_vs_time_only": lrt,
        "p_lrt": p_lrt,
    }


def logit_multi_with_time(df, x_cols, y_col, time_col="time_Ma"):
    """Logistic y ~ X + time dummies. LRT vs time-only (df = len(X))."""
    cols = list(x_cols) + [y_col, time_col]
    d = df[cols].dropna()
    if len(d) < 30 or d[y_col].nunique() < 2:
        return None
    y = d[y_col].to_numpy(dtype=float)
    Xpred = np.column_stack([d[c].to_numpy(dtype=float) for c in x_cols])
    times = d[time_col].to_numpy(dtype=float)
    uniq = np.sort(np.unique(times))
    if uniq.size < 2:
        return None
    dummies = np.column_stack([(times == t).astype(float) for t in uniq[1:]])
    X_full = np.column_stack([Xpred, dummies])
    X_red = dummies
    k = Xpred.shape[1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = LogisticRegression(
            penalty=None, solver="lbfgs", max_iter=4000
        ).fit(X_full, y)
        red = LogisticRegression(
            penalty=None, solver="lbfgs", max_iter=4000
        ).fit(X_red, y)
    coefs = [float(full.coef_[0, i]) for i in range(k)]
    p_full = np.clip(full.predict_proba(X_full)[:, 1], 1e-12, 1 - 1e-12)
    p_red = np.clip(red.predict_proba(X_red)[:, 1], 1e-12, 1 - 1e-12)
    ll_full = float(np.sum(y * np.log(p_full) + (1 - y) * np.log(1 - p_full)))
    ll_red = float(np.sum(y * np.log(p_red) + (1 - y) * np.log(1 - p_red)))
    lrt = 2.0 * (ll_full - ll_red)
    p_lrt = float(chi2.sf(max(lrt, 0.0), k))
    return {
        "n": int(len(y)),
        "n_on": int(y.sum()),
        "coefs": coefs,
        "odds_ratios": [float(np.exp(c)) for c in coefs],
        "lrt_vs_time_only": lrt,
        "p_lrt": p_lrt,
        "df": k,
    }


def mantel_haenszel(tables):
    return cup.mantel_haenszel(tables)


def tertile_code(x: pd.Series) -> pd.Series:
    """Pooled rank tertiles: 0 = low, 1 = mid, 2 = high. Ties use average rank."""
    r = x.rank(method="average", pct=True)
    out = pd.Series(np.ones(len(x), dtype=int), index=x.index)
    out[r <= (1.0 / 3.0)] = 0
    out[r > (2.0 / 3.0)] = 2
    out[~np.isfinite(x.to_numpy(dtype=float))] = -1
    return out


def mh_binary_exposure(df, y_col, exposed_mask, time_col="time_Ma"):
    """MH OR for y vs binary exposure, stratified by time.

    tables are (a,b,c,d) = on|exposed, off|exposed, on|unexp, off|unexp.
    OR > 1 means exposed blocks are more likely y=1.
    """
    tables = []
    n_exp = n_un = n_on_exp = n_on_un = 0
    for _, s in df.groupby(time_col):
        y = s[y_col].to_numpy(dtype=int)
        exp = exposed_mask.loc[s.index].to_numpy(dtype=bool)
        a = int(np.sum(exp & (y == 1)))
        b = int(np.sum(exp & (y == 0)))
        c = int(np.sum(~exp & (y == 1)))
        d = int(np.sum(~exp & (y == 0)))
        if (a + b + c + d) <= 1:
            continue
        tables.append((a, b, c, d))
        n_exp += a + b
        n_un += c + d
        n_on_exp += a
        n_on_un += c
    mh = mantel_haenszel(tables) if tables else None
    p_exp = (n_on_exp / n_exp) if n_exp else np.nan
    p_un = (n_on_un / n_un) if n_un else np.nan
    return {
        "mh": mh,
        "n_exposed": n_exp,
        "n_unexposed": n_un,
        "n_on_exposed": n_on_exp,
        "n_on_unexposed": n_on_un,
        "p_exposed": p_exp,
        "p_unexposed": p_un,
        "n_strata": 0 if mh is None else mh["n_strata"],
    }


def zone_blocks(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (time_Ma, zone_id)."""
    rows = []
    for (t, z), s in df.groupby(["time_Ma", "zone_id"], sort=False):
        n_points = int(len(s))
        n_bab = int(s["bab_on"].sum())
        bab_zone = int(n_bab > 0)
        net = s.loc[s["dilat_in_network"] == 1, "dilat_rate_1e-15_s"].to_numpy(
            dtype=float
        )
        n_net = int(np.isfinite(net).sum())
        if n_net == 0:
            med_dilat_net = np.nan
            ext_zone = 0  # all rigid / zero
            all_rigid = 1
        else:
            med_dilat_net = float(np.nanmedian(net))
            ext_zone = int(np.isfinite(med_dilat_net) and med_dilat_net > EXT_MEDIAN_THRESH)
            all_rigid = 0
        med_dilat_all = float(np.nanmedian(s["dilat_rate_1e-15_s"].to_numpy(dtype=float)))
        rows.append(
            {
                "time_Ma": float(t),
                "zone_id": int(z),
                "n_points": n_points,
                "n_bab_on": n_bab,
                "bab_zone": bab_zone,
                "bab_frac": n_bab / n_points if n_points else np.nan,
                "n_in_network": n_net,
                "all_rigid": all_rigid,
                "med_dilat_network": med_dilat_net,
                "med_dilat_all": med_dilat_all,
                "ext_zone": ext_zone,
                "seafloor_age_Ma": float(s["seafloor_age_Ma"].median()),
                "trench_zone_width_km": float(s["trench_zone_width_km"].median()),
                "dist_nearest_trench_edge_km": float(
                    s["dist_nearest_trench_edge_km"].median()
                ),
                "lon": float(s["lon"].median()),
                "lat": float(s["lat"].median()),
            }
        )
    z = pd.DataFrame(rows)
    z["width_tertile"] = tertile_code(z["trench_zone_width_km"])
    z["edge_tertile"] = tertile_code(z["dist_nearest_trench_edge_km"])
    z["near_edge"] = (
        z["dist_nearest_trench_edge_km"] < EDGE_NEAR_KM
    ).astype(int)
    return z


def univariate_table(z, y_col, level):
    rows = []
    for pred, pretty in [
        ("trench_zone_width_km", "trench-zone width (km)"),
        ("dist_nearest_trench_edge_km", "distance to trench edge (km)"),
        ("seafloor_age_Ma", "seafloor age (Ma)"),
    ]:
        x, y = _finite_pair(z[pred], z[y_col])
        rho, p_s, n = spearman_safe(x, y)
        logit = logit_univariate(x, y)
        timed = logit_with_time(z, pred, y_col)
        rows.append(
            {
                "level": level,
                "predictor": pretty,
                "col": pred,
                "n": n,
                "n_on": int(y.sum()) if n else 0,
                "spearman_rho": rho,
                "spearman_p": p_s,
                "logit_OR": logit["odds_ratio_per_unit"] if logit else np.nan,
                "logit_p_lrt": logit["p_lrt"] if logit else np.nan,
                "timed_n": timed["n"] if timed else 0,
                "timed_n_on": timed["n_on"] if timed else 0,
                "timed_coef": timed["coef"] if timed else np.nan,
                "timed_OR": timed["odds_ratio_per_unit"] if timed else np.nan,
                "timed_lrt": timed["lrt_vs_time_only"] if timed else np.nan,
                "timed_p": timed["p_lrt"] if timed else np.nan,
            }
        )
    return pd.DataFrame(rows)


def quantile_split_table(z, y_col):
    """Width tertiles (narrow vs wide) and edge <1200 vs far, plus MH."""
    out = {}
    w = z["width_tertile"]
    y = z[y_col]
    for k, lab in [(0, "narrow"), (1, "mid"), (2, "wide")]:
        m = w == k
        n = int(m.sum())
        n_on = int(y[m].sum()) if n else 0
        cuts = z.loc[m, "trench_zone_width_km"]
        out[f"width_{lab}"] = {
            "n": n,
            "n_on": n_on,
            "p_on": (n_on / n) if n else np.nan,
            "w_min": float(cuts.min()) if n else np.nan,
            "w_max": float(cuts.max()) if n else np.nan,
            "w_med": float(cuts.median()) if n else np.nan,
        }
    # MH narrow (T1) vs wide (T3); drop mid
    sub = z.loc[z["width_tertile"].isin([0, 2])].copy()
    exposed = sub["width_tertile"] == 0  # narrow = Schellart "should be more extensional"
    out["mh_width_narrow_vs_wide"] = mh_binary_exposure(sub, y_col, exposed)

    near = z["near_edge"] == 1
    far = z["near_edge"] == 0
    for lab, m in [("near", near), ("far", far)]:
        n = int(m.sum())
        n_on = int(y[m].sum()) if n else 0
        out[f"edge_{lab}"] = {
            "n": n,
            "n_on": n_on,
            "p_on": (n_on / n) if n else np.nan,
        }
    out["mh_edge_near_vs_far"] = mh_binary_exposure(z, y_col, z["near_edge"] == 1)

    e = z["edge_tertile"]
    for k, lab in [(0, "nearT"), (1, "midT"), (2, "farT")]:
        m = e == k
        n = int(m.sum())
        n_on = int(y[m].sum()) if n else 0
        cuts = z.loc[m, "dist_nearest_trench_edge_km"]
        out[f"edge_{lab}"] = {
            "n": n,
            "n_on": n_on,
            "p_on": (n_on / n) if n else np.nan,
            "e_min": float(cuts.min()) if n else np.nan,
            "e_max": float(cuts.max()) if n else np.nan,
            "e_med": float(cuts.median()) if n else np.nan,
        }
    return out


def mw_on_off(z, pred, y_col):
    on = z.loc[z[y_col] == 1, pred].to_numpy(dtype=float)
    off = z.loc[z[y_col] == 0, pred].to_numpy(dtype=float)
    on = on[np.isfinite(on)]
    off = off[np.isfinite(off)]
    if len(on) < 3 or len(off) < 3:
        return {
            "n_on": len(on),
            "n_off": len(off),
            "med_on": np.nan,
            "med_off": np.nan,
            "U": np.nan,
            "p": np.nan,
        }
    U, p = mannwhitneyu(on, off, alternative="two-sided")
    return {
        "n_on": int(len(on)),
        "n_off": int(len(off)),
        "med_on": float(np.median(on)),
        "med_off": float(np.median(off)),
        "U": float(U),
        "p": float(p),
    }


def region_zone_table(df0, z0):
    """0 Ma boxes: point-level vs the zone-time flags of overlapping zones."""
    rows = []
    for name, box in cup.SANITY_REGIONS.items():
        m = cup.region_mask(df0, box)
        s = df0.loc[m]
        n = int(m.sum())
        n_bab = int(s["bab_on"].sum()) if n else 0
        zids = sorted(s["zone_id"].unique().tolist()) if n else []
        zb = z0.loc[z0["zone_id"].isin(zids)] if zids else z0.iloc[0:0]
        rows.append(
            {
                "region": name,
                "n_points": n,
                "n_bab_on_pts": n_bab,
                "frac_bab_pts": (n_bab / n) if n else np.nan,
                "n_zones": int(len(zids)),
                "zone_ids": ",".join(str(int(i)) for i in zids),
                "bab_zone": int(zb["bab_zone"].max()) if len(zb) else np.nan,
                "ext_zone": int(zb["ext_zone"].max()) if len(zb) else np.nan,
                "med_width": float(zb["trench_zone_width_km"].median())
                if len(zb)
                else np.nan,
                "med_edge": float(zb["dist_nearest_trench_edge_km"].median())
                if len(zb)
                else np.nan,
                "med_dilat_net": float(zb["med_dilat_network"].median())
                if len(zb)
                else np.nan,
                "n_in_network_pts": int((s["dilat_in_network"] == 1).sum()) if n else 0,
            }
        )
    return pd.DataFrame(rows)


def submap_width_edge(join_csv: Path):
    """0 Ma SubMap UPS_3 vs Müller width/edge. Independent of bab_on."""
    if not join_csv.is_file():
        return None
    j = pd.read_csv(join_csv)
    if "UPS_3" not in j.columns:
        j["UPS_3"] = h12.recode_ups3(j["Geodetic_UPS"])
    j = j.loc[j["UPS_3"].isin(["C", "N", "E", "SS"])].copy()
    ce = j.loc[j["UPS_3"].isin(["C", "E"])].copy()
    ce["score"] = ce["UPS_3"].map(CE_SCORE)
    rows = []
    for col, pretty in [
        ("trench_zone_width_km", "trench-zone width (km)"),
        ("dist_nearest_trench_edge_km", "distance to trench edge (km)"),
        ("seafloor_age_Ma", "seafloor age (Ma)"),
    ]:
        r, p, n = spearman_safe(ce[col], ce["score"])
        rows.append({"predictor": pretty, "rho": r, "p": p, "n": n})
    meds = {}
    for cls in UPS3_ORDER:
        meds[cls] = {
            "n": int((j["UPS_3"] == cls).sum()),
            "width": h12.median_n(j, "trench_zone_width_km", cls),
            "edge": h12.median_n(j, "dist_nearest_trench_edge_km", cls),
            "age": h12.median_n(
                j.loc[np.isfinite(j["seafloor_age_Ma"])], "seafloor_age_Ma", cls
            ),
        }
    n_c = int((ce["UPS_3"] == "C").sum())
    n_e = int((ce["UPS_3"] == "E").sum())
    return {
        "n_defined": int(len(j)),
        "n_ce": int(len(ce)),
        "n_c": n_c,
        "n_e": n_e,
        "spearman": rows,
        "meds": meds,
        "join": j,
        "ce": ce,
    }


def savefig(fig, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"Wrote {path}")


def _strip_box(ax, data, colors, labels, ylabel, title, hline=None):
    rng = np.random.default_rng(20260903)
    bp = ax.boxplot(
        data,
        positions=np.arange(len(data)),
        widths=0.58,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="0.1", linewidth=1.6),
        whiskerprops=dict(color="0.3", linewidth=1.0),
        capprops=dict(color="0.3", linewidth=1.0),
        boxprops=dict(linewidth=0.8, edgecolor="0.25"),
        zorder=2,
    )
    for patch, color, vals in zip(bp["boxes"], colors, data):
        patch.set_facecolor(color)
        patch.set_alpha(0.40 if len(vals) else 0.08)
    for i, (vals, color) in enumerate(zip(data, colors)):
        if len(vals) == 0:
            continue
        x = i + rng.uniform(-0.18, 0.18, size=len(vals))
        ax.scatter(
            x,
            vals,
            s=14,
            c=color,
            alpha=0.55,
            edgecolors="0.15",
            linewidths=0.25,
            zorder=3,
        )
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=8)
    if hline is not None:
        ax.axhline(hline, color="0.55", linewidth=0.8, linestyle="--", zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _bar_p(ax, labels, p_on, ns, colors, ylabel, title, expected_note=None):
    x = np.arange(len(labels))
    ax.bar(x, p_on, color=colors, edgecolor="0.2", linewidth=0.7, width=0.7, zorder=2)
    ymax = max([p for p in p_on if np.isfinite(p)] + [0.05])
    for i, (p, n) in enumerate(zip(p_on, ns)):
        if not np.isfinite(p):
            continue
        ax.text(
            i,
            p + 0.03 * ymax,
            f"{p:.3f}\nn={n}",
            ha="center",
            va="bottom",
            fontsize=8,
            color="0.15",
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=8)
    ax.set_ylim(0.0, max(0.25, ymax * 1.35))
    if expected_note:
        ax.text(
            0.02,
            0.98,
            expected_note,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            color="0.3",
        )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fig_width_vs_bab(z, qA, path):
    _setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.0))
    on = z.loc[z["bab_zone"] == 1, "trench_zone_width_km"].to_numpy(dtype=float)
    off = z.loc[z["bab_zone"] == 0, "trench_zone_width_km"].to_numpy(dtype=float)
    on = on[np.isfinite(on)]
    off = off[np.isfinite(off)]
    _strip_box(
        axes[0],
        [off, on],
        ["0.65", "#0A9396"],
        [
            f"bab_zone=0\nn={len(off)}",
            f"bab_zone=1\nn={len(on)}",
        ],
        "Median trench-zone width (km)",
        "Zone–time width vs reconstructed back-arc (A)",
    )
    labs = ["narrow\nT1", "mid\nT2", "wide\nT3"]
    p_on = [qA[f"width_{k}"]["p_on"] for k in ("narrow", "mid", "wide")]
    ns = [qA[f"width_{k}"]["n"] for k in ("narrow", "mid", "wide")]
    _bar_p(
        axes[1],
        labs,
        p_on,
        ns,
        ["#94D2BD", "#E9D8A6", "#BB3E03"],
        r"$P(\mathrm{bab\_zone}=1)$",
        "Width tertiles (pooled zone–time)",
        expected_note="Schellart 2007: T1 (narrow) should be highest",
    )
    fig.suptitle(
        "Layer A — trench-zone width vs MOR back-arc (blocked zone × time)\n"
        "Müller et al. (2019) topologies; reconstruction-internal",
        y=1.03,
        fontsize=12,
    )
    fig.tight_layout()
    savefig(fig, path)


def fig_edge_vs_bab(z, qA, path):
    _setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.0))
    on = z.loc[z["bab_zone"] == 1, "dist_nearest_trench_edge_km"].to_numpy(dtype=float)
    off = z.loc[z["bab_zone"] == 0, "dist_nearest_trench_edge_km"].to_numpy(dtype=float)
    on = on[np.isfinite(on)]
    off = off[np.isfinite(off)]
    _strip_box(
        axes[0],
        [off, on],
        ["0.65", "#0A9396"],
        [
            f"bab_zone=0\nn={len(off)}",
            f"bab_zone=1\nn={len(on)}",
        ],
        "Median distance to nearest trench edge (km)",
        "Zone–time edge distance vs reconstructed back-arc (A)",
    )
    axes[0].axhline(EDGE_NEAR_KM, color="#9B2226", ls="--", lw=1.0, zorder=1)
    labs = ["near T1", "mid T2", "far T3"]
    p_on = [qA[k]["p_on"] for k in ("edge_nearT", "edge_midT", "edge_farT")]
    ns = [qA[k]["n"] for k in ("edge_nearT", "edge_midT", "edge_farT")]
    _bar_p(
        axes[1],
        labs,
        p_on,
        ns,
        ["#94D2BD", "#E9D8A6", "#BB3E03"],
        r"$P(\mathrm{bab\_zone}=1)$",
        "Edge-distance tertiles (pooled zone–time)",
        expected_note="Schellart 2007: T1 (near-edge) should be highest",
    )
    fig.suptitle(
        "Layer A — distance to trench edge vs MOR back-arc (blocked zone × time)\n"
        "Median edge distance at zone level is partly a width proxy",
        y=1.03,
        fontsize=12,
    )
    fig.tight_layout()
    savefig(fig, path)


def fig_dilatation(z, qB, path):
    _setup_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.0))
    net = z.loc[z["all_rigid"] == 0].copy()
    wcode = net["width_tertile"]
    data = []
    cols = []
    labs = []
    for k, lab, c in [
        (0, "narrow T1", "#94D2BD"),
        (1, "mid T2", "#E9D8A6"),
        (2, "wide T3", "#BB3E03"),
    ]:
        vals = net.loc[wcode == k, "med_dilat_network"].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        data.append(vals)
        cols.append(c)
        labs.append(f"{lab}\nn={len(vals)}")
    _strip_box(
        axes[0],
        data,
        cols,
        labs,
        r"Median dilatation ($10^{-15}$ s$^{-1}$)",
        "Layer B — dilatation vs width tertile\n(network blocks only)",
        hline=0.0,
    )

    near = net.loc[net["near_edge"] == 1, "med_dilat_network"].to_numpy(dtype=float)
    far = net.loc[net["near_edge"] == 0, "med_dilat_network"].to_numpy(dtype=float)
    near = near[np.isfinite(near)]
    far = far[np.isfinite(far)]
    _strip_box(
        axes[1],
        [near, far],
        ["#0A9396", "#BB3E03"],
        [
            f"near <{int(EDGE_NEAR_KM)}\nn={len(near)}",
            f"far ≥{int(EDGE_NEAR_KM)}\nn={len(far)}",
        ],
        r"Median dilatation ($10^{-15}$ s$^{-1}$)",
        "Layer B — dilatation vs 1200 km edge split\n(network blocks only)",
        hline=0.0,
    )

    labs = ["narrow T1", "mid T2", "wide T3"]
    p_on = [qB[f"width_{k}"]["p_on"] for k in ("narrow", "mid", "wide")]
    ns = [qB[f"width_{k}"]["n"] for k in ("narrow", "mid", "wide")]
    _bar_p(
        axes[2],
        labs,
        p_on,
        ns,
        ["#94D2BD", "#E9D8A6", "#BB3E03"],
        r"$P(\mathrm{ext\_zone}=1)$",
        "Layer B — ext_zone by width tertile\n(all blocks; rigid → 0)",
        expected_note="Schellart: narrow should be highest",
    )
    fig.suptitle(
        "Layer B — deforming-network dilatation (blocked zone × time)\n"
        r"ext_zone = 1 iff median dilatation among network samples $> 0$; "
        "Hellenic-style extension is this layer, not MOR bab_on",
        y=1.04,
        fontsize=12,
    )
    fig.tight_layout()
    savefig(fig, path)


def fig_map_0ma(df0, path):
    _setup_style()
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        proj = ccrs.Robinson()
        data_crs = ccrs.PlateCarree()
        fig, axes = plt.subplots(
            1, 2, figsize=(13.4, 5.6), subplot_kw={"projection": proj}
        )
        for ax in axes:
            ax.set_global()
            try:
                ax.add_feature(
                    cfeature.LAND, facecolor="0.88", edgecolor="none", zorder=0
                )
                ax.add_feature(
                    cfeature.COASTLINE, linewidth=0.4, edgecolor="0.35", zorder=1
                )
            except Exception as err:
                print(f"WARNING: NaturalEarth features unavailable ({err})")
            ax.gridlines(draw_labels=False, linewidth=0.3, color="0.7", linestyle=":")
        transform = data_crs
    except Exception as err:
        print(f"WARNING: cartopy unavailable ({err}); plain lon/lat map.")
        fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.6))
        for ax in axes:
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.set_aspect("equal", adjustable="box")
        transform = None

    off = df0[df0["bab_on"] == 0]
    on = df0[df0["bab_on"] == 1]
    fields = [
        (
            axes[0],
            "trench_zone_width_km",
            "Trench-zone width (km)",
            "cividis",
        ),
        (
            axes[1],
            "dist_nearest_trench_edge_km",
            "Distance to nearest trench edge (km)",
            "plasma",
        ),
    ]
    for ax, col, label, cmap in fields:
        kw = dict(s=9, linewidths=0, zorder=3)
        if transform is not None:
            sc = ax.scatter(
                df0["lon"],
                df0["lat"],
                c=df0[col],
                cmap=cmap,
                transform=transform,
                **kw,
            )
            ax.scatter(
                on["lon"],
                on["lat"],
                facecolors="none",
                edgecolors="#0A9396",
                transform=transform,
                s=28,
                linewidths=0.9,
                zorder=4,
            )
        else:
            sc = ax.scatter(df0["lon"], df0["lat"], c=df0[col], cmap=cmap, **kw)
            ax.scatter(
                on["lon"],
                on["lat"],
                facecolors="none",
                edgecolors="#0A9396",
                s=28,
                linewidths=0.9,
                zorder=4,
            )
        cb = fig.colorbar(sc, ax=ax, shrink=0.72, pad=0.03)
        cb.set_label(label)
        ax.set_title(label + "\nteal rings = bab_on")
    fig.suptitle(
        "0 Ma trench samples — width / edge with MOR back-arc overlay",
        y=1.02,
        fontsize=12,
    )
    fig.tight_layout()
    savefig(fig, path)


def mh_line(mh_pack, exposed_lab, unexp_lab):
    mh = mh_pack["mh"]
    if mh is None:
        return (
            f"Mantel–Haenszel {exposed_lab} vs {unexp_lab}: not defined "
            f"(n_exp={mh_pack['n_exposed']}, n_un={mh_pack['n_unexposed']})"
        )
    return (
        f"Mantel–Haenszel OR ({exposed_lab} vs {unexp_lab}, y=1 more likely if exposed) "
        f"= {fmt_num(mh['or_mh'])} (χ²={fmt_num(mh['chi2'], 2)}, p={fmt_p(mh['p'])}, "
        f"{mh['n_strata']} time strata); "
        f"P(y=1|{exposed_lab})={fmt_num(mh_pack['p_exposed'])} "
        f"(n={mh_pack['n_exposed']}, n_on={mh_pack['n_on_exposed']}) vs "
        f"P(y=1|{unexp_lab})={fmt_num(mh_pack['p_unexposed'])} "
        f"(n={mh_pack['n_unexposed']}, n_on={mh_pack['n_on_unexposed']})"
    )


def spearman_md(tab: pd.DataFrame) -> str:
    lines = [
        "| predictor | n | n_on | Spearman ρ | p | logistic OR (no time) | LRT p |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in tab.iterrows():
        lines.append(
            f"| {r['predictor']} | {int(r['n'])} | {int(r['n_on'])} | "
            f"{fmt_rho(r['spearman_rho'])} | {fmt_p(r['spearman_p'])} | "
            f"{fmt_num(r['logit_OR'], 4)} | {fmt_p(r['logit_p_lrt'])} |"
        )
    return "\n".join(lines)


def timed_md(tab: pd.DataFrame) -> str:
    lines = [
        "| predictor | n | n_on | coef | OR per unit | LRT vs time-only χ² | p |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in tab.iterrows():
        lines.append(
            f"| {r['predictor']} | {int(r['timed_n'])} | {int(r['timed_n_on'])} | "
            f"{fmt_num(r['timed_coef'], 5)} | {fmt_num(r['timed_OR'], 5)} | "
            f"{fmt_num(r['timed_lrt'], 2)} | {fmt_p(r['timed_p'])} |"
        )
    return "\n".join(lines)


def sign_vs_schellart(rho, kind):
    """kind: 'width' (Schellart wants ρ < 0 vs extension) or 'edge' (ρ < 0)."""
    if not np.isfinite(rho):
        return "NA"
    if rho < 0:
        return (
            "sign matches Schellart (narrower / nearer-edge more extensional)"
            if kind in ("width", "edge")
            else "negative"
        )
    if rho > 0:
        return (
            "sign is opposite Schellart (wider / farther-from-edge more extensional)"
        )
    return "zero"


def write_report(
    *,
    n_rows,
    n_times,
    z,
    tabA,
    tabB,
    tabBnet,
    qA,
    qB,
    qBnet,
    mwA_w,
    mwA_e,
    mwB_w,
    mwB_e,
    jointA,
    jointB,
    rho_we,
    rho_we_p,
    rho_we_n,
    contB,
    sub,
    region_tab,
    z0,
):
    n_blocks = int(len(z))
    n_bab = int(z["bab_zone"].sum())
    n_ext = int(z["ext_zone"].sum())
    n_net = int((z["all_rigid"] == 0).sum())
    n_age = int(np.isfinite(z["seafloor_age_Ma"]).sum())

    rA_w = tabA.loc[tabA["col"] == "trench_zone_width_km"].iloc[0]
    rA_e = tabA.loc[tabA["col"] == "dist_nearest_trench_edge_km"].iloc[0]
    rA_a = tabA.loc[tabA["col"] == "seafloor_age_Ma"].iloc[0]
    rB_w = tabB.loc[tabB["col"] == "trench_zone_width_km"].iloc[0]
    rB_e = tabB.loc[tabB["col"] == "dist_nearest_trench_edge_km"].iloc[0]

    # Winner language: only if blocked numbers are strong AND sign matches.
    def strong(row):
        return (
            np.isfinite(row["spearman_rho"])
            and abs(row["spearman_rho"]) >= 0.25
            and row["spearman_p"] < 1e-6
            and np.isfinite(row["timed_p"])
            and row["timed_p"] < 1e-3
        )

    def claim(row, kind, layer):
        if not np.isfinite(row["spearman_rho"]):
            return f"{layer}: no finite ρ."
        match = row["spearman_rho"] < 0
        sig = np.isfinite(row["spearman_p"]) and row["spearman_p"] < 0.05
        if strong(row) and match:
            return (
                f"{layer}: blocked association is strong and in the Schellart direction "
                f"(ρ={fmt_rho(row['spearman_rho'])}, p={fmt_p(row['spearman_p'])}). "
                "Still reconstruction-internal; not independent geology."
            )
        if strong(row) and not match:
            return (
                f"{layer}: blocked association is strong but the **sign is opposite** "
                f"Schellart (ρ={fmt_rho(row['spearman_rho'])}, p={fmt_p(row['spearman_p'])}). "
                "Not a Schellart confirmation."
            )
        extra = (
            sign_vs_schellart(row["spearman_rho"], kind) + "; "
            if sig
            else "not distinguishable from null at 0.05; "
        )
        return (
            f"{layer}: blocked numbers are not a Schellart confirmation "
            f"(ρ={fmt_rho(row['spearman_rho'])}, p={fmt_p(row['spearman_p'])}; "
            f"{extra}). No winner claimed."
        )

    region_md = [
        "| region | n pts | n bab_on | bab_zone | ext_zone | med width (km) | med edge (km) | med dilat net |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in region_tab.iterrows():
        region_md.append(
            f"| {r['region']} | {int(r['n_points'])} | {int(r['n_bab_on_pts'])} | "
            f"{r['bab_zone']} | {r['ext_zone']} | {fmt_num(r['med_width'], 0)} | "
            f"{fmt_num(r['med_edge'], 0)} | {fmt_num(r['med_dilat_net'], 3)} |"
        )

    def tertile_md(q, prefix="width"):
        if prefix == "width":
            lines = [
                "| tertile | n blocks | n on | P(on) | width min–max (km) | median width (km) |",
                "|---|---:|---:|---:|---|---:|",
            ]
            for lab in ("narrow", "mid", "wide"):
                d = q[f"width_{lab}"]
                lines.append(
                    f"| {lab} | {d['n']} | {d['n_on']} | {fmt_num(d['p_on'])} | "
                    f"{fmt_num(d['w_min'], 0)}–{fmt_num(d['w_max'], 0)} | "
                    f"{fmt_num(d['w_med'], 0)} |"
                )
            return "\n".join(lines)
        lines = [
            "| tertile | n blocks | n on | P(on) | edge min–max (km) | median edge (km) |",
            "|---|---:|---:|---:|---|---:|",
        ]
        for lab, key in (("near T1", "edge_nearT"), ("mid T2", "edge_midT"), ("far T3", "edge_farT")):
            d = q[key]
            lines.append(
                f"| {lab} | {d['n']} | {d['n_on']} | {fmt_num(d['p_on'])} | "
                f"{fmt_num(d['e_min'], 0)}–{fmt_num(d['e_max'], 0)} | "
                f"{fmt_num(d['e_med'], 0)} |"
            )
        return "\n".join(lines)

    jointA_md = "joint logistic not defined"
    if jointA is not None:
        names = ["width", "edge", "age"]
        bits = ", ".join(
            f"{nm} coef={fmt_num(c, 5)} (OR={fmt_num(o, 5)})"
            for nm, c, o in zip(names, jointA["coefs"], jointA["odds_ratios"])
        )
        jointA_md = (
            f"n={jointA['n']} (n_on={jointA['n_on']}); {bits}; "
            f"LRT vs time-only χ²={fmt_num(jointA['lrt_vs_time_only'], 2)} "
            f"(df={jointA['df']}, p={fmt_p(jointA['p_lrt'])})"
        )
    jointB_md = "joint logistic not defined"
    if jointB is not None:
        names = ["width", "edge", "age"]
        bits = ", ".join(
            f"{nm} coef={fmt_num(c, 5)} (OR={fmt_num(o, 5)})"
            for nm, c, o in zip(names, jointB["coefs"], jointB["odds_ratios"])
        )
        jointB_md = (
            f"n={jointB['n']} (n_on={jointB['n_on']}); {bits}; "
            f"LRT vs time-only χ²={fmt_num(jointB['lrt_vs_time_only'], 2)} "
            f"(df={jointB['df']}, p={fmt_p(jointB['p_lrt'])})"
        )

    if sub is None:
        sub_section = "Joined SubMap file not found; 0 Ma UPS_3 sanity skipped.\n"
    else:
        spear_lines = [
            "| predictor | Spearman ρ vs C vs E (C=+1, E=−1) | p | n |",
            "|---|---:|---:|---:|",
        ]
        for row in sub["spearman"]:
            spear_lines.append(
                f"| {row['predictor']} | {fmt_rho(row['rho'])} | {fmt_p(row['p'])} | {row['n']} |"
            )
        med_lines = [
            "| UPS_3 | n | median width (km) | median dist-to-edge (km) | median age (Ma) |",
            "|---|---:|---|---|---|",
        ]
        for cls in UPS3_ORDER:
            m = sub["meds"][cls]
            med_lines.append(
                f"| {cls} | {m['n']} | {h12.med_cell(m['width'], 'km')} | "
                f"{h12.med_cell(m['edge'], 'km')} | {h12.med_cell(m['age'], 'Ma')} |"
            )
        # Schellart vs SubMap: width ρ vs C=+1 should be >0 if wide → C
        wrow = sub["spearman"][0]
        erow = sub["spearman"][1]
        sub_width_note = (
            "Schellart at 0 Ma vs independent UPS: ρ(width, C=+1/E=−1) "
            f"= {fmt_rho(wrow['rho'])} (p={fmt_p(wrow['p'])}, n={wrow['n']}). "
            "A Schellart width effect (wide → more C) would be **positive**. "
        )
        if np.isfinite(wrow["rho"]) and wrow["p"] >= 0.05:
            sub_width_note += "This 0 Ma width result is weak (not distinguishable from null at 0.05)."
        elif np.isfinite(wrow["rho"]) and wrow["rho"] > 0:
            sub_width_note += "Sign matches wide → C, but see p and n before any claim."
        elif np.isfinite(wrow["rho"]):
            sub_width_note += "Sign is opposite wide → C (wider associated with E)."
        sub_section = f"""Join: `{JOIN_CSV}` (nearest Müller trench sample, already computed).
Width and edge come from Müller2019; UPS_3 is SubMap geodetic strain.
This table does **not** use `bab_on`. SS held out of C vs E. Undefined dropped.
C = C+C-SS, E = E+E-SS. Score C=+1, E=−1.

- n defined UPS_3 after join: **{sub['n_defined']}**
- C vs E n: **{sub['n_ce']}** (C n={sub['n_c']}, E n={sub['n_e']})

{chr(10).join(spear_lines)}

{chr(10).join(med_lines)}

{sub_width_note}
ρ(edge, C/E) = {fmt_rho(erow['rho'])} (p={fmt_p(erow['p'])}, n={erow['n']}).
A Schellart edge effect (near-edge → more E) would be **positive** ρ vs C=+1
(farther from edge → more C). Do not treat this as a Cenozoic confirmation.
"""

    n0 = int(len(z0))
    n0_bab = int(z0["bab_zone"].sum())
    n0_ext = int(z0["ext_zone"].sum())

    md = f"""# Schellart width / edge test — Cenozoic, zone–time blocked

Built by `scripts/schellart_width_edge.py` from
`{PROXY_CSV}`. Trench kinematics were **not** re-extracted.
Müller `v_T` vs `v_OP` is **not** used (collinear in this model).

- **n trench samples:** {n_rows}
- **n times:** {n_times} (0–60 Ma, 1 Myr)
- **n zone–time blocks:** {n_blocks}
- **n bab_zone=1:** {n_bab} ({100.0 * n_bab / n_blocks:.1f}%)
- **n ext_zone=1:** {n_ext} ({100.0 * n_ext / n_blocks:.1f}%)
- **n blocks with a deforming-network sample:** {n_net}
- **n blocks with finite median age:** {n_age}

## Circularity (read this first)

Layer A (`bab_on` / `bab_zone`) is a back-arc MOR flag taken from the **same**
Müller et al. (2019) topologies that supply the trench catalog. Layer B is
dilatation inside those deforming networks. Both are reconstruction-internal.
The Cenozoic test asks whether *this model's* back-arc / extension geometry
follows Schellart et al. (2007) width and edge patterns, not whether nature
does. Present-day SubMap geodetic UPS remains the only UPS here that is
independent of the plate model; that 0 Ma table is a sanity check.

## Schellart prediction under test

Schellart, Freeman, Stegman, Moresi & May (2007, *Nature*): rollback and
back-arc opening concentrate on **narrow** slabs and **near a lateral slab
edge**; interiors of wide slabs (≳4000 km) are slow or advancing / more
compressive. They associated near-edge retreat with a length scale of order
**1200 km**.

Operationally, at zone–time:

- narrower `trench_zone_width_km` → higher P(bab_zone) and higher P(ext_zone)
  → Spearman ρ(**width**, extension flag) should be **negative**
- smaller `dist_nearest_trench_edge_km` → more extension
  → Spearman ρ(**edge distance**, extension flag) should be **negative**
- near-edge split: P(on | median edge < 1200 km) > P(on | ≥ 1200 km)

`trench_zone_width_km` is the along-trench length of a spatially connected
polyline sharing a subducting-plate ID (see `extract_trench_kinematics_0Ma.py`).
It is a trench-zone width, not a tomographic slab width. `zone_id` is unique
only within a time slice.

Zone-level median edge distance is mechanically related to width (a uniform
sample on a zone of length W has median distance-to-edge ≈ W/4).
Spearman ρ(width, edge) at zone–time = {fmt_rho(rho_we)} (p={fmt_p(rho_we_p)}, n={rho_we_n}).
Treat the two predictors as overlapping geometric information.

## Inference unit and outcomes

**Unit:** `groupby(time_Ma, zone_id)` — **{n_blocks}** blocks. Point-level n≈76k
is spatially autocorrelated along the trench and is not used for inference.

**Layer A — `bab_zone`.** 1 if any tessellated sample in the block has `bab_on=1`
(qualifying overriding-plate MOR; see `cenozoic_age_gate.md`). 0 otherwise.

**Layer B — `ext_zone`.** Among samples with `dilat_in_network==1`, take the
median `dilat_rate_1e-15_s`. `ext_zone = 1` if that median **> {EXT_MEDIAN_THRESH:g}**
(threshold: strictly positive dilatation, units 10⁻¹⁵ s⁻¹). If the block has
no network samples (all rigid plates, dilatation recorded as zero /
`dilat_in_network==0`), `ext_zone = 0`. Negative median dilatation (network
shortening) is also `ext_zone = 0`. Hellenic / Aegean extension is a deforming
network in Müller2019, not a MOR: it should be A=0, B=1.

A network-only sensitivity drops `all_rigid` blocks so rigid zeros do not
dilute layer B.

Age is a covariate already known to be weak at this blocking
(`cenozoic_age_gate.md`); it is not the headline.

## 0 Ma region check (A vs B)

{chr(10).join(region_md)}

Expected: Mariana / Tonga / Scotia / Andaman **bab_zone=1**. Andes / Cascadia /
MAT_Mexico **bab_zone=0**. Hellenic **bab_zone=0** and **ext_zone=1**.
The `> 0` cut is literal: Andes and MAT_Mexico can land on `ext_zone=1` with
median network dilatation of order 10⁻³ (10⁻¹⁵ s⁻¹), which is not
Hellenic-style extension (Hellenic median is larger; see table). Continuous
median dilatation (below) does not use that binary cut.

0 Ma zone–time blocks: n={n0}, bab_zone=1 in {n0_bab}, ext_zone=1 in {n0_ext}.

## Layer A — MOR back-arc (`bab_zone`)

n_blocks = **{n_blocks}**, n_on = **{n_bab}**.

Median width bab_zone=1 vs 0: {fmt_num(mwA_w['med_on'], 0)} vs {fmt_num(mwA_w['med_off'], 0)} km
(U={fmt_num(mwA_w['U'], 1)}, p={fmt_p(mwA_w['p'])}, n_on={mwA_w['n_on']}, n_off={mwA_w['n_off']}).

Median edge distance bab_zone=1 vs 0: {fmt_num(mwA_e['med_on'], 0)} vs {fmt_num(mwA_e['med_off'], 0)} km
(U={fmt_num(mwA_e['U'], 1)}, p={fmt_p(mwA_e['p'])}).

### Spearman / univariate logistic (zone–time)

{spearman_md(tabA)}

{claim(rA_w, 'width', 'Width vs bab_zone')}
{claim(rA_e, 'edge', 'Edge vs bab_zone')}
Age covariate: ρ={fmt_rho(rA_a['spearman_rho'])}, p={fmt_p(rA_a['spearman_p'])} (n={int(rA_a['n'])}).

### Logistic with time dummies

`y ~ predictor + time factor`. Coefficient is after absorbing time.

{timed_md(tabA)}

Joint `bab_zone ~ width + edge + age + time dummies`: {jointA_md}.
Width and edge are collinear; joint coefficients are descriptive only.

### Quantile splits (Schellart geometry)

Width tertiles are pooled rank tertiles of zone–time median width (average-rank
ties). Narrow vs wide MH drops the middle tertile. Edge split is median
distance < {int(EDGE_NEAR_KM)} km vs ≥ {int(EDGE_NEAR_KM)} km.

{tertile_md(qA, 'width')}

{mh_line(qA['mh_width_narrow_vs_wide'], 'narrow T1', 'wide T3')}

Edge <{int(EDGE_NEAR_KM)} km: P(bab_zone)={fmt_num(qA['edge_near']['p_on'])} (n={qA['edge_near']['n']}, n_on={qA['edge_near']['n_on']});
≥{int(EDGE_NEAR_KM)} km: P(bab_zone)={fmt_num(qA['edge_far']['p_on'])} (n={qA['edge_far']['n']}, n_on={qA['edge_far']['n_on']}).

{mh_line(qA['mh_edge_near_vs_far'], 'near <1200 km', 'far ≥1200 km')}

{tertile_md(qA, 'edge')}

An MH OR > 1 for *narrow* vs *wide* (or *near* vs *far*) would match Schellart
(exposed class more extensional).

## Layer B — dilatation (`ext_zone`)

Threshold: median network dilatation **> {EXT_MEDIAN_THRESH:g}** × 10⁻¹⁵ s⁻¹.
n_blocks = **{n_blocks}**, n_on = **{n_ext}**, of which network blocks = {n_net}.

Median width ext_zone=1 vs 0: {fmt_num(mwB_w['med_on'], 0)} vs {fmt_num(mwB_w['med_off'], 0)} km
(U={fmt_num(mwB_w['U'], 1)}, p={fmt_p(mwB_w['p'])}, n_on={mwB_w['n_on']}, n_off={mwB_w['n_off']}).

Median edge distance ext_zone=1 vs 0: {fmt_num(mwB_e['med_on'], 0)} vs {fmt_num(mwB_e['med_off'], 0)} km
(U={fmt_num(mwB_e['U'], 1)}, p={fmt_p(mwB_e['p'])}).

### Spearman / univariate logistic (all blocks; rigid → ext_zone=0)

{spearman_md(tabB)}

{claim(rB_w, 'width', 'Width vs ext_zone')}
{claim(rB_e, 'edge', 'Edge vs ext_zone')}

### Logistic with time dummies

{timed_md(tabB)}

Joint `ext_zone ~ width + edge + age + time dummies`: {jointB_md}.

### Quantile splits

{tertile_md(qB, 'width')}

{mh_line(qB['mh_width_narrow_vs_wide'], 'narrow T1', 'wide T3')}

Edge <{int(EDGE_NEAR_KM)} km: P(ext_zone)={fmt_num(qB['edge_near']['p_on'])} (n={qB['edge_near']['n']});
≥{int(EDGE_NEAR_KM)} km: P(ext_zone)={fmt_num(qB['edge_far']['p_on'])} (n={qB['edge_far']['n']}).

{mh_line(qB['mh_edge_near_vs_far'], 'near <1200 km', 'far ≥1200 km')}

### Continuous dilatation (network blocks; no binary cut)

Spearman of median network dilatation (signed, 10⁻¹⁵ s⁻¹) vs predictors.
Positive ρ means larger width / edge distance goes with more extension.

| predictor | n | Spearman ρ | p |
|---|---:|---:|---:|
| trench-zone width (km) | {contB['width'][2]} | {fmt_rho(contB['width'][0])} | {fmt_p(contB['width'][1])} |
| distance to trench edge (km) | {contB['edge'][2]} | {fmt_rho(contB['edge'][0])} | {fmt_p(contB['edge'][1])} |
| seafloor age (Ma) | {contB['age'][2]} | {fmt_rho(contB['age'][0])} | {fmt_p(contB['age'][1])} |

Schellart direction for this signed rate is still negative ρ (narrower /
nearer-edge → more positive dilatation).

### Sensitivity: network blocks only (drop all-rigid)

n={int(len(z[z['all_rigid']==0]))}, n_on={int(z.loc[z['all_rigid']==0, 'ext_zone'].sum())}.

{spearman_md(tabBnet)}

{tertile_md(qBnet, 'width')}

{mh_line(qBnet['mh_width_narrow_vs_wide'], 'narrow T1', 'wide T3')}
{mh_line(qBnet['mh_edge_near_vs_far'], 'near <1200 km', 'far ≥1200 km')}

## 0 Ma SubMap sanity (width / edge vs UPS_3, independent of bab_on)

{sub_section}
## Do width / edge survive blocking?

Blocked n = {n_blocks} zone–time units, still reconstruction-internal.

**Layer A width does survive blocking as an association, but not as a Schellart test.** ρ={fmt_rho(rA_w['spearman_rho'])}, p={fmt_p(rA_w['spearman_p'])}; time-dummy LRT p={fmt_p(rA_w['timed_p'])}. The sign is **opposite** Schellart: wide T3 P(bab_zone)={fmt_num(qA['width_wide']['p_on'])} vs narrow T1 {fmt_num(qA['width_narrow']['p_on'])}; MH OR (narrow vs wide) = {fmt_num(qA['mh_width_narrow_vs_wide']['mh']['or_mh']) if qA['mh_width_narrow_vs_wide']['mh'] else np.nan} (p={fmt_p(qA['mh_width_narrow_vs_wide']['mh']['p']) if qA['mh_width_narrow_vs_wide']['mh'] else np.nan}). Wider reconstructed zones are more likely to host a back-arc MOR in this model.

**Layer A edge association also survives blocking, also opposite Schellart** (ρ={fmt_rho(rA_e['spearman_rho'])}, time-dummy p={fmt_p(rA_e['timed_p'])}; far tertile P(bab_zone) higher than near). The Schellart 1200 km cut is **degenerate at zone-median** (n_far={qA['edge_far']['n']} of {n_blocks}); almost every zone has median edge distance < 1200 km, so that split is not an interior-vs-edge test.

**Layer B width/edge do not survive time dummies** (width LRT p={fmt_p(rB_w['timed_p'])}; edge LRT p={fmt_p(rB_e['timed_p'])}). Continuous median dilatation vs width among network blocks is null (ρ={fmt_rho(contB['width'][0])}, p={fmt_p(contB['width'][1])}).

0 Ma SubMap C vs E (independent of bab_on): width ρ={fmt_rho(sub['spearman'][0]['rho']) if sub else np.nan} (p={fmt_p(sub['spearman'][0]['p']) if sub else 'NA'}, n={sub['spearman'][0]['n'] if sub else 0}), matching the weak 0 Ma width result in `cenozoic_age_gate.md`.

No hypothesis winner. The Schellart *direction* does not survive this blocking.

## Figures

- `{FIG_WIDTH}`
- `{FIG_EDGE}`
- `{FIG_DILAT}`
- `{FIG_MAP}` (0 Ma map; teal rings = `bab_on`)

## How to rerun

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
    /workspace/subduction-test/scripts/schellart_width_edge.py
```

## Caveats

- Reconstruction-internal. Opening histories of back-arc basins are encoded
  in the topologies that define `bab_on`.
- `bab_zone=1` if *any* point in the zone is `bab_on`, so a long connected
  trench with a local MOR (e.g. Mariana Trough on an Izu–Bonin–Mariana
  polyline) flags the whole zone.
- Median edge distance at zone level is not a within-zone interior-vs-edge
  contrast; it covaries with width.
- Hellenic / Okinawa-style continental extension is layer B, not A.
- Absolute frame / H1 vs H2 is not tested here. See
  `results/submap_h1h2_frames.md` for SubMap Vupn/Vtn across frames.
"""
    REPORT_MD.write_text(md)
    print(f"Wrote {REPORT_MD}")


def main() -> int:
    _setup_style()
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not PROXY_CSV.is_file():
        raise SystemExit(f"Missing {PROXY_CSV}")

    need = [
        "time_Ma",
        "lon",
        "lat",
        "zone_id",
        "seafloor_age_Ma",
        "trench_zone_width_km",
        "dist_nearest_trench_edge_km",
        "bab_on",
        "dist_to_backarc_ridge_km",
        "dilat_rate_1e-15_s",
        "dilat_in_network",
    ]
    df = pd.read_csv(PROXY_CSV, usecols=lambda c: c in set(need) or True)
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit(f"Proxy CSV missing columns: {missing}")

    n_rows = int(len(df))
    n_times = int(df["time_Ma"].nunique())
    print(f"Loaded {n_rows} samples, {n_times} times from {PROXY_CSV}")

    z = zone_blocks(df)
    print(
        f"zone-time blocks n={len(z)} bab_zone=1 n={int(z.bab_zone.sum())} "
        f"ext_zone=1 n={int(z.ext_zone.sum())} network n={int((z.all_rigid==0).sum())}"
    )

    rho_we, rho_we_p, rho_we_n = spearman_safe(
        z["trench_zone_width_km"], z["dist_nearest_trench_edge_km"]
    )

    tabA = univariate_table(z, "bab_zone", "A")
    tabB = univariate_table(z, "ext_zone", "B")
    znet = z.loc[z["all_rigid"] == 0].copy()
    tabBnet = univariate_table(znet, "ext_zone", "B_network")

    qA = quantile_split_table(z, "bab_zone")
    qB = quantile_split_table(z, "ext_zone")
    qBnet = quantile_split_table(znet, "ext_zone")

    mwA_w = mw_on_off(z, "trench_zone_width_km", "bab_zone")
    mwA_e = mw_on_off(z, "dist_nearest_trench_edge_km", "bab_zone")
    mwB_w = mw_on_off(z, "trench_zone_width_km", "ext_zone")
    mwB_e = mw_on_off(z, "dist_nearest_trench_edge_km", "ext_zone")

    jointA = logit_multi_with_time(
        z,
        ["trench_zone_width_km", "dist_nearest_trench_edge_km", "seafloor_age_Ma"],
        "bab_zone",
    )
    jointB = logit_multi_with_time(
        z,
        ["trench_zone_width_km", "dist_nearest_trench_edge_km", "seafloor_age_Ma"],
        "ext_zone",
    )
    contB = {
        "width": spearman_safe(znet["trench_zone_width_km"], znet["med_dilat_network"]),
        "edge": spearman_safe(
            znet["dist_nearest_trench_edge_km"], znet["med_dilat_network"]
        ),
        "age": spearman_safe(znet["seafloor_age_Ma"], znet["med_dilat_network"]),
    }

    df0 = df.loc[df["time_Ma"] == 0].copy()
    z0 = z.loc[z["time_Ma"] == 0].copy()
    region_tab = region_zone_table(df0, z0)
    print("\n=== 0 Ma region zone flags ===")
    print(region_tab.to_string(index=False))

    sub = submap_width_edge(JOIN_CSV)

    fig_width_vs_bab(z, qA, FIG_WIDTH)
    fig_edge_vs_bab(z, qA, FIG_EDGE)
    fig_dilatation(z, qB, FIG_DILAT)
    fig_map_0ma(df0, FIG_MAP)

    write_report(
        n_rows=n_rows,
        n_times=n_times,
        z=z,
        tabA=tabA,
        tabB=tabB,
        tabBnet=tabBnet,
        qA=qA,
        qB=qB,
        qBnet=qBnet,
        mwA_w=mwA_w,
        mwA_e=mwA_e,
        mwB_w=mwB_w,
        mwB_e=mwB_e,
        jointA=jointA,
        jointB=jointB,
        rho_we=rho_we,
        rho_we_p=rho_we_p,
        rho_we_n=rho_we_n,
        contB=contB,
        sub=sub,
        region_tab=region_tab,
        z0=z0,
    )

    print("\nLayer A Spearman:")
    print(tabA[["predictor", "n", "spearman_rho", "spearman_p", "timed_p"]].to_string(index=False))
    print("\nLayer B Spearman:")
    print(tabB[["predictor", "n", "spearman_rho", "spearman_p", "timed_p"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
