#!/usr/bin/env python3
"""0 Ma H1 vs H2 kinematic test from SubMap's own M56 velocities.

Uses SubMap M56+ (Cerpa, Lallemand & Heuret 2025 / Sub-DATA v7), independent
of the Müller 2019 trench=OP collinearity. Numbers are computed from the
CSVs; nothing is invented.

The compact UPS table does not store absolute V_UP / V_T. Those trench-normal
scalars (Vupn1, Vtn1; MORVEL56-NNR) are merged from the full Sub-DATA sheet
in the same dump, by Short_Name.

Run:
    /workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
        /workspace/subduction-test/scripts/submap_m56_h1h2.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr, t as student_t

ROOT = Path("/workspace/subduction-test")
SUB_CSV = ROOT / "data" / "submap" / "submap_geodetic_UPS.csv"
FULL_CSV = ROOT / "data" / "submap" / "submap_transects_full.csv"
JOIN_CSV = ROOT / "results" / "joined_submap_muller2019_0Ma.csv"
STATS_MD = ROOT / "results" / "submap_m56_h1h2.md"
FIG_VCN = ROOT / "results" / "fig_submap_vcn_vs_UPS.png"
FIG_VSN = ROOT / "results" / "fig_submap_vsn_vs_UPS.png"
FIG_VDN = ROOT / "results" / "fig_submap_vdn_vs_UPS.png"
FIG_AGE = ROOT / "results" / "fig_age_width_vs_UPS.png"
FIG_ABS = ROOT / "results" / "fig_submap_Vupn_Vtn_vs_UPS.png"

SENTINEL = -999.0

UPS3_MAP = {
    "C": "C",
    "C-SS": "C",
    "E": "E",
    "E-SS": "E",
    "Neutral": "N",
    "SS": "SS",
    "SS-C": "SS",
    "SS-E": "SS",
    "Undefined": "Undefined",
}

UPS3_ORDER = ["C", "N", "E", "SS"]
UPS3_COLORS = {
    "C": "#9B2226",
    "N": "#6C757D",
    "E": "#005F73",
    "SS": "#C9A227",
}

# C=1, N=0, E=-1 aligns with SubMap vdn (shortening > 0).
CNE_SCORE = {"C": 1.0, "N": 0.0, "E": -1.0}
CE_SCORE = {"C": 1.0, "E": -1.0}

M56_VEL_MM = ["M56_vc", "M56_vcn", "M56_vs", "M56_vsn", "M56_vd", "M56_vdn"]
ABS_VEL_MM = ["Vupn1", "Vtn1", "Vsubn1"]


def _setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "none",
            "axes.edgecolor": "0.2",
            "axes.labelcolor": "0.1",
            "text.color": "0.1",
            "xtick.color": "0.2",
            "ytick.color": "0.2",
            "axes.grid": False,
            "font.size": 11,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "figure.dpi": 140,
            "savefig.dpi": 200,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def recode_ups3(geodetic_ups: pd.Series) -> pd.Series:
    mapped = geodetic_ups.map(UPS3_MAP)
    unknown = mapped.isna() & geodetic_ups.notna()
    if unknown.any():
        extras = sorted(geodetic_ups[unknown].unique().tolist())
        raise ValueError(f"Unexpected Geodetic_UPS values: {extras}")
    return mapped


def sentinel_to_nan(s: pd.Series) -> pd.Series:
    out = pd.to_numeric(s, errors="coerce")
    return out.mask(out == SENTINEL)


def box_strip(ax, df, ycol, order, colors, ylabel, title, annotate_n=True, zero_line=True):
    rng = np.random.default_rng(20260902)
    data = []
    plotted_colors = []
    ns = []
    for cat in order:
        vals = df.loc[df["UPS_3"] == cat, ycol].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        data.append(vals)
        plotted_colors.append(colors[cat])
        ns.append(len(vals))

    bp = ax.boxplot(
        data,
        positions=np.arange(len(order)),
        widths=0.58,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="0.1", linewidth=1.6),
        whiskerprops=dict(color="0.3", linewidth=1.0),
        capprops=dict(color="0.3", linewidth=1.0),
        boxprops=dict(linewidth=0.8, edgecolor="0.25"),
        zorder=2,
    )
    for patch, color, vals in zip(bp["boxes"], plotted_colors, data):
        patch.set_facecolor(color)
        patch.set_alpha(0.40 if len(vals) else 0.08)

    for i, (vals, color) in enumerate(zip(data, plotted_colors)):
        if len(vals) == 0:
            continue
        x = i + rng.uniform(-0.18, 0.18, size=len(vals))
        ax.scatter(
            x,
            vals,
            s=22,
            c=color,
            alpha=0.75,
            edgecolors="0.15",
            linewidths=0.35,
            zorder=3,
        )

    ax.set_xticks(np.arange(len(order)))
    if annotate_n:
        labels = [f"{cat}\nn={n}" for cat, n in zip(order, ns)]
    else:
        labels = list(order)
    ax.set_xticklabels(labels, rotation=0)
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=8)
    if zero_line:
        ax.axhline(0.0, color="0.55", linewidth=0.8, linestyle="--", zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return ns


def savefig(fig, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"Wrote {path}")


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_rho(r: float) -> str:
    if not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def fmt_num(x: float, digits: int = 3) -> str:
    if not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def spearman_safe(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 3:
        return np.nan, np.nan, n
    if np.unique(x[m]).size < 2 or np.unique(y[m]).size < 2:
        return np.nan, np.nan, n
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p), n


def ols_simple(x, y):
    """Simple OLS y ~ 1 + x. Returns slope, intercept, R2, slope p, n."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    n = int(len(x))
    out = {"n": n, "slope": np.nan, "intercept": np.nan, "r2": np.nan, "p": np.nan, "se": np.nan}
    if n < 3:
        return out
    X = np.column_stack([np.ones(n), x])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    out["intercept"] = float(beta[0])
    out["slope"] = float(beta[1])
    out["r2"] = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    if n > 2 and ss_res >= 0:
        mse = ss_res / (n - 2)
        try:
            xtx_inv = np.linalg.inv(X.T @ X)
            se = float(np.sqrt(mse * xtx_inv[1, 1]))
            out["se"] = se
            if se > 0 and np.isfinite(se):
                tstat = out["slope"] / se
                out["p"] = float(2.0 * student_t.sf(abs(tstat), n - 2))
        except np.linalg.LinAlgError:
            pass
    return out


def ols_multi(Xcols, y):
    """OLS y ~ 1 + Xcols (list of 1d arrays). Returns r2, n, betas."""
    cols = [np.asarray(c, dtype=float) for c in Xcols]
    y = np.asarray(y, dtype=float)
    m = np.isfinite(y)
    for c in cols:
        m &= np.isfinite(c)
    n = int(m.sum())
    if n < 3 + len(cols):
        return {"n": n, "r2": np.nan, "betas": None}
    Y = y[m]
    X = np.column_stack([np.ones(n)] + [c[m] for c in cols])
    beta, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    yhat = X @ beta
    ss_res = float(np.sum((Y - yhat) ** 2))
    ss_tot = float(np.sum((Y - Y.mean()) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    return {"n": n, "r2": r2, "betas": beta}


def logistic_mle(Xcols, y):
    """Bernoulli logistic y ~ 1 + Xcols. McFadden R2, AIC, betas."""
    cols = [np.asarray(c, dtype=float) for c in Xcols]
    y = np.asarray(y, dtype=float)
    m = np.isfinite(y)
    for c in cols:
        m &= np.isfinite(c)
    y = y[m]
    n = int(len(y))
    k = 1 + len(cols)
    out = {
        "n": n,
        "ll": np.nan,
        "ll0": np.nan,
        "r2_mcf": np.nan,
        "aic": np.nan,
        "betas": None,
        "converged": False,
        "separation": False,
    }
    if n < 8 or y.min() == y.max():
        return out
    X = np.column_stack([np.ones(n)] + [c[m] for c in cols])

    def nll(beta):
        z = np.clip(X @ beta, -30.0, 30.0)
        p = 1.0 / (1.0 + np.exp(-z))
        p = np.clip(p, 1e-12, 1.0 - 1e-12)
        return -np.sum(y * np.log(p) + (1.0 - y) * np.log(1.0 - p))

    pbar = float(np.clip(y.mean(), 1e-12, 1.0 - 1e-12))
    ll0 = float(np.sum(y * np.log(pbar) + (1.0 - y) * np.log(1.0 - pbar)))
    res = minimize(nll, np.zeros(k), method="BFGS")
    ll = float(-res.fun)
    out["ll"] = ll
    out["ll0"] = ll0
    out["r2_mcf"] = float(1.0 - ll / ll0) if ll0 != 0 else np.nan
    out["aic"] = float(2 * k - 2 * ll)
    out["betas"] = res.x
    out["converged"] = bool(res.success)
    z = np.clip(X @ res.x, -30.0, 30.0)
    p = 1.0 / (1.0 + np.exp(-z))
    out["separation"] = bool(np.any((p < 1e-5) | (p > 1.0 - 1e-5)))
    return out


def median_n(df, col, cls, cls_col="UPS_3"):
    vals = df.loc[df[cls_col] == cls, col].to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.nan, 0
    return float(np.median(vals)), int(len(vals))


def med_cell(pair, unit="cm/yr"):
    med, n = pair
    if n == 0 or not np.isfinite(med):
        return "NA (n=0)"
    return f"{med:.3f} {unit} (n={n})"


def main() -> int:
    _setup_style()
    RESULTS = ROOT / "results"
    RESULTS.mkdir(parents=True, exist_ok=True)

    sub = pd.read_csv(SUB_CSV)
    n_raw = len(sub)
    if n_raw != 260:
        print(f"WARNING: compact SubMap n={n_raw} (expected 260)")

    need = [
        "Short_Name",
        "Zone",
        "Trench_Name",
        "Lon",
        "Lat",
        "UPN",
        "Geodetic_UPS",
        "M56_vcn",
        "M56_vsn",
        "M56_vdn",
        "M56_vd",
    ]
    missing = [c for c in need if c not in sub.columns]
    if missing:
        raise SystemExit(f"Compact CSV missing columns: {missing}")

    full = pd.read_csv(FULL_CSV)
    abs_keep = ["Short_Name"] + ABS_VEL_MM
    missing_abs = [c for c in abs_keep if c not in full.columns]
    if missing_abs:
        raise SystemExit(f"Full Sub-DATA sheet missing absolute-velocity columns: {missing_abs}")
    sub = sub.merge(full[abs_keep], on="Short_Name", how="left", validate="1:1")

    for col in M56_VEL_MM + ABS_VEL_MM:
        if col not in sub.columns:
            continue
        mm = sentinel_to_nan(sub[col])
        sub[col] = mm
        sub[f"{col}_cm_yr"] = mm / 10.0

    sub["UPS_3"] = recode_ups3(sub["Geodetic_UPS"])

    # Drop Undefined / vd sentinel. Neutral (vd=0) is kept.
    vd_ok = sub["M56_vd"].notna()
    defined = sub.loc[vd_ok & (sub["Geodetic_UPS"] != "Undefined")].copy()
    n_defined = len(defined)
    n_dropped = n_raw - n_defined

    # Identity diagnostics (mm/yr, then also cm/yr residuals).
    ident = defined["M56_vcn"] - (defined["M56_vsn"] + defined["M56_vdn"])
    ident_abs = ident.abs()
    n_ident = int(np.isfinite(ident).sum())
    n_ident_exact = int((ident == 0).sum())
    n_ident_le1 = int((ident_abs <= 1).sum())
    n_ident_gt1 = int((ident_abs > 1).sum())
    ident_med = float(np.nanmedian(ident_abs))
    ident_max = float(np.nanmax(ident_abs)) if n_ident else np.nan
    outlier_rows = defined.loc[ident_abs > 1, ["Short_Name", "Trench_Name", "Geodetic_UPS", "M56_vcn", "M56_vsn", "M56_vdn"]].copy()
    outlier_rows["vcn-(vsn+vdn)"] = ident.loc[ident_abs > 1].to_numpy()

    vdn_from_abs = defined["Vupn1"] - defined["Vtn1"]
    abs_resid = (defined["M56_vdn"] - vdn_from_abs).abs()
    n_abs = int((defined["Vupn1"].notna() & defined["Vtn1"].notna() & defined["M56_vdn"].notna()).sum())
    n_abs_le1 = int((abs_resid <= 1).sum())
    abs_med = float(np.nanmedian(abs_resid))
    r_vup_vt, p_vup_vt, n_vup_vt = spearman_safe(defined["Vupn1"], defined["Vtn1"])

    ups3_counts = (
        defined["UPS_3"].value_counts().reindex(["C", "N", "E", "SS"]).fillna(0).astype(int)
    )
    geo_counts = defined["Geodetic_UPS"].value_counts()

    # ---- figures: vcn / vsn / vdn vs UPS_3 ----
    plot_df = defined.copy()

    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    box_strip(
        ax,
        plot_df,
        "M56_vcn_cm_yr",
        UPS3_ORDER,
        UPS3_COLORS,
        ylabel=r"SubMap $v_{cn}$ (cm/yr)",
        title=(
            "SubMap M56 trench-normal convergence vs UPS$_3$\n"
            r"$v_{cn}$: far-field SP–UP convergence; shortening-positive $v_{dn}$ identity $v_{cn}\approx v_{sn}+v_{dn}$"
        ),
    )
    ax.set_xlabel("UPS$_3$ (C = C+C-SS; E = E+E-SS; N = Neutral; SS = SS/SS-C/SS-E)")
    fig.tight_layout()
    savefig(fig, FIG_VCN)

    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    box_strip(
        ax,
        plot_df,
        "M56_vsn_cm_yr",
        UPS3_ORDER,
        UPS3_COLORS,
        ylabel=r"SubMap $v_{sn}$ (cm/yr)",
        title=(
            "SubMap M56 trench-normal subduction rate vs UPS$_3$\n"
            r"$v_{sn}$: SP consumption at trench (hinge-related); $v_{sn}\approx v_{cn}-v_{dn}$"
        ),
    )
    ax.set_xlabel("UPS$_3$ (C = C+C-SS; E = E+E-SS; N = Neutral; SS = SS/SS-C/SS-E)")
    fig.tight_layout()
    savefig(fig, FIG_VSN)

    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    box_strip(
        ax,
        plot_df,
        "M56_vdn_cm_yr",
        UPS3_ORDER,
        UPS3_COLORS,
        ylabel=r"SubMap $v_{dn}$ (cm/yr)",
        title=(
            "SubMap M56 trench-normal OP deformation vs UPS$_3$  (sanity; circular)\n"
            r"shortening $\Rightarrow v_{dn}>0$; extension $\Rightarrow v_{dn}<0$; UPS is classified from $O_{vd}$/$v_d$"
        ),
    )
    ax.set_xlabel("UPS$_3$ (C = C+C-SS; E = E+E-SS; N = Neutral; SS = SS/SS-C/SS-E)")
    fig.tight_layout()
    savefig(fig, FIG_VDN)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.2), sharex=True)
    box_strip(
        axes[0],
        plot_df,
        "Vupn1_cm_yr",
        UPS3_ORDER,
        UPS3_COLORS,
        ylabel=r"SubMap $V_{\mathrm{UPn}1}$ (cm/yr)",
        title="H1 proxy: trench-normal absolute upper-plate velocity\nMORVEL56-NNR (frame 1); + = oceanward (toward trench)",
    )
    box_strip(
        axes[1],
        plot_df,
        "Vtn1_cm_yr",
        UPS3_ORDER,
        UPS3_COLORS,
        ylabel=r"SubMap $V_{Tn}1$ (cm/yr)",
        title="H2 proxy: trench-normal absolute trench velocity\nMORVEL56-NNR (frame 1); + = oceanward (rollback)",
    )
    axes[0].set_xlabel("UPS$_3$")
    axes[1].set_xlabel("UPS$_3$")
    fig.suptitle("SubMap own absolute velocities (not in the compact UPS CSV)", fontsize=12, y=1.02)
    fig.tight_layout()
    savefig(fig, FIG_ABS)

    # ---- Müller join covariates (independent of SubMap class) ----
    join_ok = JOIN_CSV.is_file()
    join = None
    age_df = None
    if join_ok:
        join = pd.read_csv(JOIN_CSV)
        if "UPS_3" not in join.columns:
            join["UPS_3"] = recode_ups3(join["Geodetic_UPS"])
        join_def = join.loc[join["UPS_3"].isin(["C", "N", "E", "SS"])].copy()
        fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.0), sharex=True)
        age_df = join_def.loc[np.isfinite(join_def["seafloor_age_Ma"])].copy()
        box_strip(
            axes[0],
            age_df,
            "seafloor_age_Ma",
            UPS3_ORDER,
            UPS3_COLORS,
            ylabel="Seafloor age (Ma)",
            title="Müller 2019 seafloor age at trench",
            zero_line=False,
        )
        box_strip(
            axes[1],
            join_def,
            "trench_zone_width_km",
            UPS3_ORDER,
            UPS3_COLORS,
            ylabel="Trench-zone width (km)",
            title="Along-trench connected-zone width",
            zero_line=False,
        )
        box_strip(
            axes[2],
            join_def,
            "dist_nearest_trench_edge_km",
            UPS3_ORDER,
            UPS3_COLORS,
            ylabel="Distance to nearest trench edge (km)",
            title="Distance along trench to nearest edge",
            zero_line=False,
        )
        for ax in axes:
            ax.set_xlabel("UPS$_3$")
        fig.suptitle(
            "Müller 2019 H2/H3 covariates vs SubMap UPS$_3$ (0 Ma join ≤ 200 km)\n"
            "Independent of SubMap class and of SubMap M56 velocities",
            fontsize=12,
            y=1.04,
        )
        fig.tight_layout()
        savefig(fig, FIG_AGE)

    # ---- stats: C/N/E and C vs E ----
    cne = defined.loc[defined["UPS_3"].isin(["C", "N", "E"])].copy()
    cne["score"] = cne["UPS_3"].map(CNE_SCORE)
    ce = defined.loc[defined["UPS_3"].isin(["C", "E"])].copy()
    ce["score"] = ce["UPS_3"].map(CE_SCORE)
    ce["yC"] = (ce["UPS_3"] == "C").astype(float)

    predictors = [
        ("M56_vcn_cm_yr", "vcn"),
        ("M56_vsn_cm_yr", "vsn"),
        ("M56_vdn_cm_yr", "vdn"),
        ("Vupn1_cm_yr", "Vupn1"),
        ("Vtn1_cm_yr", "Vtn1"),
    ]

    rows_cne = []
    rows_ce = []
    ols_ce = {}
    log_ce = {}
    for col, name in predictors:
        r, p, n = spearman_safe(cne[col], cne["score"])
        rows_cne.append((name, r, p, n))
        r2, p2, n2 = spearman_safe(ce[col], ce["score"])
        rows_ce.append((name, r2, p2, n2))
        ols_ce[name] = ols_simple(ce[col].to_numpy(), ce["yC"].to_numpy())
        log_ce[name] = logistic_mle([ce[col].to_numpy()], ce["yC"].to_numpy())

    ols_both = ols_multi(
        [ce["M56_vcn_cm_yr"].to_numpy(), ce["M56_vsn_cm_yr"].to_numpy()],
        ce["yC"].to_numpy(),
    )
    log_both = logistic_mle(
        [ce["M56_vcn_cm_yr"].to_numpy(), ce["M56_vsn_cm_yr"].to_numpy()],
        ce["yC"].to_numpy(),
    )
    log_null_n = int(ce["yC"].notna().sum())

    # Medians
    meds = {}
    for cls in UPS3_ORDER:
        meds[cls] = {
            "vcn": median_n(defined, "M56_vcn_cm_yr", cls),
            "vsn": median_n(defined, "M56_vsn_cm_yr", cls),
            "vdn": median_n(defined, "M56_vdn_cm_yr", cls),
            "Vupn1": median_n(defined, "Vupn1_cm_yr", cls),
            "Vtn1": median_n(defined, "Vtn1_cm_yr", cls),
        }

    # Müller covariate stats
    muller_lines = []
    muller_med_lines = []
    if join_ok:
        jdef = join.loc[join["UPS_3"].isin(["C", "N", "E", "SS"])].copy()
        jcne = jdef.loc[jdef["UPS_3"].isin(["C", "N", "E"])].copy()
        jcne["score"] = jcne["UPS_3"].map(CNE_SCORE)
        jce = jdef.loc[jdef["UPS_3"].isin(["C", "E"])].copy()
        jce["score"] = jce["UPS_3"].map(CE_SCORE)
        covs = [
            ("seafloor_age_Ma", "seafloor_age_Ma", "Ma"),
            ("trench_zone_width_km", "trench_zone_width_km", "km"),
            ("dist_nearest_trench_edge_km", "dist_nearest_trench_edge_km", "km"),
        ]
        muller_lines.append(
            "| predictor | sample | Spearman ρ vs UPS score | p | n |"
        )
        muller_lines.append("|---|---|---:|---:|---:|")
        for col, name, unit in covs:
            r, p, n = spearman_safe(jcne[col], jcne["score"])
            muller_lines.append(
                f"| `{name}` | C/N/E (C=1,N=0,E=−1) | {fmt_rho(r)} | {fmt_p(p)} | {n} |"
            )
            r, p, n = spearman_safe(jce[col], jce["score"])
            muller_lines.append(
                f"| `{name}` | C vs E (C=1,E=−1) | {fmt_rho(r)} | {fmt_p(p)} | {n} |"
            )
        muller_med_lines.append("| UPS_3 | median age (Ma) | median zone width (km) | median dist-to-edge (km) |")
        muller_med_lines.append("|---|---|---|---|")
        for cls in UPS3_ORDER:
            a = median_n(jdef.loc[np.isfinite(jdef["seafloor_age_Ma"])], "seafloor_age_Ma", cls)
            w = median_n(jdef, "trench_zone_width_km", cls)
            d = median_n(jdef, "dist_nearest_trench_edge_km", cls)
            muller_med_lines.append(
                f"| {cls} | {med_cell(a, 'Ma')} | {med_cell(w, 'km')} | {med_cell(d, 'km')} |"
            )
        n_join_def = len(jdef)
        n_join_age = int(np.isfinite(jdef["seafloor_age_Ma"]).sum())
        n_join_raw = len(join)
    else:
        n_join_def = 0
        n_join_age = 0
        n_join_raw = 0

    # UPN is not a velocity
    upn_tab = pd.crosstab(defined["UPS_3"], defined["UPN"])

    # Plain-language comparison of vcn vs vsn on C vs E (non-circular).
    rho_vcn_ce = rows_ce[0][1]
    rho_vsn_ce = rows_ce[1][1]
    rho_vdn_ce = rows_ce[2][1]
    r2_vcn = ols_ce["vcn"]["r2"]
    r2_vsn = ols_ce["vsn"]["r2"]
    mcf_vcn = log_ce["vcn"]["r2_mcf"]
    mcf_vsn = log_ce["vsn"]["r2_mcf"]
    abs_vcn = abs(rho_vcn_ce) if np.isfinite(rho_vcn_ce) else np.nan
    abs_vsn = abs(rho_vsn_ce) if np.isfinite(rho_vsn_ce) else np.nan

    stronger = "neither stands out"
    if np.isfinite(abs_vcn) and np.isfinite(abs_vsn):
        gap = abs(abs_vcn - abs_vsn)
        if gap < 0.08:
            stronger = (
                f"|ρ| differs by only {gap:.3f}, so vcn and vsn are similar in rank strength "
                "on this C vs E sample"
            )
        elif abs_vcn > abs_vsn:
            stronger = (
                f"|ρ(vcn, C/E)| = {abs_vcn:.3f} exceeds |ρ(vsn, C/E)| = {abs_vsn:.3f} "
                f"(gap {gap:.3f})"
            )
        else:
            stronger = (
                f"|ρ(vsn, C/E)| = {abs_vsn:.3f} exceeds |ρ(vcn, C/E)| = {abs_vcn:.3f} "
                f"(gap {gap:.3f})"
            )

    overwhelming = (
        np.isfinite(abs_vcn)
        and np.isfinite(abs_vsn)
        and max(abs_vcn, abs_vsn) >= 0.50
        and abs(abs_vcn - abs_vsn) >= 0.25
    )

    rho_vup_ce = rows_ce[3][1]
    rho_vt_ce = rows_ce[4][1]
    abs_vup = abs(rho_vup_ce) if np.isfinite(rho_vup_ce) else np.nan
    abs_vt = abs(rho_vt_ce) if np.isfinite(rho_vt_ce) else np.nan
    h12_note = "SubMap Vupn1 and Vtn1 are only moderately rank-correlated, unlike Müller v_OP ≈ −v_T."
    if np.isfinite(abs_vup) and np.isfinite(abs_vt):
        h12_note += (
            f" On C vs E, |ρ| is {abs_vup:.3f} for Vupn1 (H1) and {abs_vt:.3f} for Vtn1 (H2)."
        )
        if max(abs_vup, abs_vt) >= 0.50 and abs(abs_vup - abs_vt) >= 0.25:
            h12_note += " That gap is large enough to note, but still not treated as a hypothesis winner here."
        else:
            h12_note += " That is not an overwhelming H1-vs-H2 separation."

    outlier_md = ""
    if len(outlier_rows):
        lines = [
            "| Short_Name | Trench_Name | Geodetic_UPS | vcn | vsn | vdn | vcn−(vsn+vdn) |",
            "|---|---|---|---:|---:|---:|---:|",
        ]
        for _, r in outlier_rows.iterrows():
            lines.append(
                f"| {r['Short_Name']} | {r['Trench_Name']} | {r['Geodetic_UPS']} | "
                f"{r['M56_vcn']:.1f} | {r['M56_vsn']:.1f} | {r['M56_vdn']:.1f} | "
                f"{r['vcn-(vsn+vdn)']:.1f} |"
            )
        outlier_md = "\n".join(lines)

    spear_cne_md = [
        "| predictor | Spearman ρ vs C/N/E (C=1, N=0, E=−1) | p | n |",
        "|---|---:|---:|---:|",
    ]
    for name, r, p, n in rows_cne:
        spear_cne_md.append(f"| {name} | {fmt_rho(r)} | {fmt_p(p)} | {n} |")

    spear_ce_md = [
        "| predictor | Spearman ρ vs C vs E (C=1, E=−1) | p | n |",
        "|---|---:|---:|---:|",
    ]
    for name, r, p, n in rows_ce:
        spear_ce_md.append(f"| {name} | {fmt_rho(r)} | {fmt_p(p)} | {n} |")

    ols_md = [
        "| predictor | OLS slope (y=1 if C, 0 if E) | R² | p (slope) | n |",
        "|---|---:|---:|---:|---:|",
    ]
    for name in ["vcn", "vsn", "vdn", "Vupn1", "Vtn1"]:
        o = ols_ce[name]
        ols_md.append(
            f"| {name} | {fmt_num(o['slope'])} | {fmt_num(o['r2'])} | {fmt_p(o['p'])} | {o['n']} |"
        )
    ols_md.append(
        f"| vcn + vsn (nested) | — | {fmt_num(ols_both['r2'])} | — | {ols_both['n']} |"
    )

    def log_row(name, L):
        sep = " yes" if L["separation"] else " no"
        return (
            f"| {name} | {fmt_num(L['r2_mcf'])} | {fmt_num(L['aic'], 1)} | "
            f"{'yes' if L['converged'] else 'no'} |{sep} | {L['n']} |"
        )

    log_md = [
        "| predictor | McFadden R² | AIC | converged | quasi-separation | n |",
        "|---|---:|---:|---|---|---:|",
        log_row("vcn", log_ce["vcn"]),
        log_row("vsn", log_ce["vsn"]),
        log_row("vdn (circular sanity)", log_ce["vdn"]),
        log_row("Vupn1", log_ce["Vupn1"]),
        log_row("Vtn1", log_ce["Vtn1"]),
        log_row("vcn + vsn (nested)", log_both),
    ]

    med_md = [
        "| UPS_3 | n | median vcn | median vsn | median vdn |",
        "|---|---:|---|---|---|",
    ]
    for cls in UPS3_ORDER:
        ncls = int(ups3_counts[cls])
        med_md.append(
            f"| {cls} | {ncls} | {med_cell(meds[cls]['vcn'])} | "
            f"{med_cell(meds[cls]['vsn'])} | {med_cell(meds[cls]['vdn'])} |"
        )

    abs_med_md = [
        "| UPS_3 | median Vupn1 (H1) | median Vtn1 (H2) |",
        "|---|---|---|",
    ]
    for cls in UPS3_ORDER:
        abs_med_md.append(
            f"| {cls} | {med_cell(meds[cls]['Vupn1'])} | {med_cell(meds[cls]['Vtn1'])} |"
        )

    geo_md = ["| Geodetic_UPS | UPS_3 | n |", "|---|---|---:|"]
    geo_to_3 = [
        ("C", "C"),
        ("C-SS", "C"),
        ("Neutral", "N"),
        ("E", "E"),
        ("E-SS", "E"),
        ("SS", "SS"),
        ("SS-C", "SS"),
        ("SS-E", "SS"),
    ]
    for g, u in geo_to_3:
        geo_md.append(f"| {g} | {u} | {int(geo_counts.get(g, 0))} |")

    upn_md_lines = ["| UPS_3 | UPN=1 | UPN=2 |", "|---|---:|---:|"]
    for cls in UPS3_ORDER:
        n1 = int(upn_tab.loc[cls, 1]) if cls in upn_tab.index and 1 in upn_tab.columns else 0
        n2 = int(upn_tab.loc[cls, 2]) if cls in upn_tab.index and 2 in upn_tab.columns else 0
        upn_md_lines.append(f"| {cls} | {n1} | {n2} |")
    upn_md = "\n".join(upn_md_lines)

    winner_line = (
        "No hypothesis is declared a winner: the vcn-vs-vsn contrast is not overwhelming "
        "on this sample."
        if not overwhelming
        else (
            "The vcn-vs-vsn contrast is numerically large on this sample; still treat it as "
            "a description of these rows, not a general H1/H2 verdict."
        )
    )

    muller_section = ""
    if join_ok:
        muller_section = f"""## Müller 2019 covariates (H2/H3; independent of SubMap class)

Joined table: `{JOIN_CSV}` (nearest trench sample, ≤ 200 km; n matched = {n_join_raw}).
These fields are **not** derived from SubMap UPS or M56 velocities.

- n with defined UPS_3 after join: **{n_join_def}**
- n with finite seafloor age: **{n_join_age}**

{chr(10).join(muller_lines)}

{chr(10).join(muller_med_lines)}

Age, trench-zone width, and distance-to-edge are the independent H2/H3-style covariates
(slab age; slab/trench width; proximity to a slab edge). They are reported here against
SubMap UPS_3; they are not a test of SubMap M56 kinematics.
"""
    else:
        muller_section = "## Müller 2019 covariates\n\nJoined file not found; age/width figure and covariate stats skipped.\n"

    stats = f"""# SubMap M56 H1 vs H2 kinematic test (0 Ma)

Values are computed from the SubMap dump and, where used, the existing Müller join.
Nothing is invented. UPS ~ vdn is circular by construction and is reported only as a sanity check.

- **Compact UPS table:** `{SUB_CSV}` (n={n_raw})
- **Full Sub-DATA sheet (absolute velocities only):** `{FULL_CSV}`
- **Script:** `/workspace/subduction-test/scripts/submap_m56_h1h2.py`
- **Units:** SubMap stores mm/yr; all velocities below are converted to **cm/yr** (÷ 10). Sentinel −999 → NaN.

## Sample

- n compact rows: **{n_raw}**
- n dropped (Undefined and/or `M56_vd` missing/−999): **{n_dropped}**
- n kept: **{n_defined}**
- C vs E test n: **{len(ce)}** (C n={int((ce['UPS_3']=='C').sum())}, E n={int((ce['UPS_3']=='E').sum())})
- C/N/E ordinal n: **{len(cne)}** (N n={int((cne['UPS_3']=='N').sum())})
- SS held out of C/N/E tests, shown separately in box plots and median tables (n={int(ups3_counts['SS'])})

### UPS_3 recode

C = C, C-SS; E = E, E-SS; N = Neutral; SS = SS, SS-C, SS-E (kept separate). Undefined dropped.

{chr(10).join(geo_md)}

`UPN` is a 1/2 upper-plate-nature code, **not** a velocity. It is not used as an H1/H2 proxy.
Among kept rows, UPN counts by UPS_3 (not used as a predictor):

{upn_md}

## Kinematic identity (documented, not assumed)

Cerpa, Lallemand & Heuret (2025, §II.1) define **relative** velocities from MORVEL56 + an arc block:

- `vc` = convergence of undeformed SP and UP interiors
- `vd` = arc-block motion relative to UP (forearc-to-back-arc deformation)
- `vs` = subduction / SP consumption at the trench (`vs = V_SP + V_trench` as vectors)
- Vector statement in the paper: `vs = vc + vd`

SubMap **signed trench-normal scalars** (shortening ⇒ `vdn > 0`, extension ⇒ `vdn < 0`) obey

**`vcn ≈ vsn + vdn`**  (equivalently **`vsn ≈ vcn − vdn`**).

On the kept rows this is a data identity, not a fit:

- n finite: **{n_ident}**
- median |vcn − (vsn + vdn)|: **{ident_med:.3f} mm/yr**
- n with exact 0 (integer mm/yr): **{n_ident_exact}**
- n with |residual| ≤ 1 mm/yr: **{n_ident_le1}**
- n with |residual| > 1 mm/yr: **{n_ident_gt1}** (max {ident_max:.1f} mm/yr)

Rows with |residual| > 1 mm/yr (left as in the database; not corrected):

{outlier_md}

Absolute trench-normal scalars in the full sheet (frame 1 = MORVEL56-NNR) satisfy

**`vdn ≈ Vupn1 − Vtn1`** (median |residual| = {abs_med:.3f} mm/yr; {n_abs_le1}/{n_abs} within 1 mm/yr)

and **`vcn ≈ Vsubn1 + Vupn1`**, **`vsn ≈ Vsubn1 + Vtn1`**.

Sign of the absolute trench-normal components implied by that identity and by SubMap’s
`vdn` convention: **positive Vupn / Vtn = oceanward** (toward the subducting plate). Then
`vdn > 0` (shortening) when the upper plate moves oceanward faster than the trench
(OP advance relative to the hinge), and `vdn < 0` (extension) when the trench moves
oceanward faster than the upper plate (rollback relative to OP).

## Which SubMap fields map to H1 / H2

The compact UPS CSV does **not** contain absolute upper-plate or trench velocities.
`UPN` is not a proxy. `vc` and `vs` are relative rates; they do not uniquely recover
`V_UP` and `V_T` without an absolute frame.

| role | compact CSV | what we use | notes |
|---|---|---|---|
| H1 (upper-plate trench-normal absolute) | not present | **`Vupn1`** from full Sub-DATA sheet, frame 1 (MORVEL56-NNR) | + = oceanward (toward trench) |
| H2 (trench-normal trench velocity / rollback) | not present | **`Vtn1`** from the same sheet / frame | + = oceanward (rollback) |
| relative convergence | `M56_vcn` | `M56_vcn` | SP vs undeformed UP |
| relative subduction / hinge consumption | `M56_vsn` | `M56_vsn` | SP vs trench |
| OP deformation (circular with UPS) | `M56_vdn` | `M56_vdn` | class is from `Ovd` / `vd` |

Spearman ρ(`Vupn1`, `Vtn1`) = {fmt_rho(r_vup_vt)} (p = {fmt_p(p_vup_vt)}, n = {n_vup_vt}).
That is **not** Müller-style collinearity (`v_T ≈ −v_OP`). SubMap’s own absolute
velocities therefore *can* separate H1 from H2, but only after merging `Vupn1`/`Vtn1`
from the full sheet. The compact CSV alone cannot.

If those absolute columns are ignored, the compact table still supports a
**non-circular relative-velocity test**: does `vcn` or `vsn` predict UPS_3 (C vs E)
better? That is not a clean H1-vs-H2 split (`vcn` mixes both plates; `vsn` mixes SP
and trench), but it does not use the class-defining `vdn`.

## Median M56 velocities by UPS_3 (cm/yr)

{chr(10).join(med_md)}

Absolute-velocity medians (same rows; cm/yr):

{chr(10).join(abs_med_md)}

## Spearman: C/N/E ordinal and C vs E

Score sign: C = +1, E = −1 (and N = 0), matching SubMap `vdn` (shortening positive).
`vdn` vs class is **circular**. `vcn` and `vsn` are the non-circular relative predictors.
`Vupn1` / `Vtn1` are the H1 / H2 absolute predictors.

{chr(10).join(spear_cne_md)}

{chr(10).join(spear_ce_md)}

## OLS and logistic (C vs E only)

Binary y = 1 if C, 0 if E. Predictors in cm/yr. Nested model is `y ~ vcn + vsn`.
`vdn` is listed as a sanity check (quasi-complete separation is expected).

{chr(10).join(ols_md)}

{chr(10).join(log_md)}

Because `vdn ≈ vcn − vsn`, the two-predictor nested model reconstructs the
class-defining deformation rate. Nested logistic McFadden R² = 1 and nested OLS R²
near the univariate `vdn` R² are therefore **circular**, not evidence that `vcn` and
`vsn` jointly explain UPS independently of `vd`.

## Which predictor is stronger (plain language)

Non-circular relative velocities, C vs E (n={len(ce)}):
Spearman ρ(vcn) = {fmt_rho(rho_vcn_ce)}, ρ(vsn) = {fmt_rho(rho_vsn_ce)};
OLS R²(vcn) = {fmt_num(r2_vcn)}, R²(vsn) = {fmt_num(r2_vsn)};
logistic McFadden R²(vcn) = {fmt_num(mcf_vcn)}, R²(vsn) = {fmt_num(mcf_vsn)}.
{stronger}. Univariate comparison only: nested `vcn+vsn` reconstructs `vdn` via the identity, so nested R² values are not an independent test (OLS R²={fmt_num(ols_both['r2'])}; logistic McFadden R²={fmt_num(log_both['r2_mcf'])}).

`vdn` vs C/E is circular (ρ = {fmt_rho(rho_vdn_ce)}) and is not a test.

H1 vs H2 from SubMap absolute velocities: {h12_note}

{winner_line}

{muller_section}
## Figures

- `{FIG_VCN}`
- `{FIG_VSN}`
- `{FIG_VDN}` (sanity; circular)
- `{FIG_ABS}` (H1/H2 absolute velocities from the full sheet)
- `{FIG_AGE if join_ok else "(age/width figure skipped — join file missing)"}`

## Caveats

- Geodetic UPS is SubMap’s classification from `Ovd` / `vd`. Any statistic that uses `vdn`
  as a predictor of UPS is recovering that definition.
- Compact-CSV columns cannot split `V_T` from `V_UP`. Absolute H1/H2 uses frame-1
  MORVEL56-NNR scalars from the full sheet of the **same** dump (`Vupn1`, `Vtn1`).
  Other frames (2–4) exist and would shift the absolute numbers; they were not mixed in.
- M56 magnitudes are stored as integers mm/yr; 1 mm/yr identity residuals are rounding.
- Pure C and pure E are small (each n=7 before recode). Recoded C and E are larger
  (C-SS and E-SS folded in) but still modest.
- The Müller join is a nearest-point geographic match, not a common segmentation.
"""

    STATS_MD.write_text(stats)
    print(f"Wrote {STATS_MD}")
    print(f"n_raw={n_raw} n_defined={n_defined} n_CE={len(ce)}")
    print("UPS_3", ups3_counts.to_dict())
    print("medians vcn", {k: meds[k]["vcn"] for k in UPS3_ORDER})
    print("medians vsn", {k: meds[k]["vsn"] for k in UPS3_ORDER})
    print("medians vdn", {k: meds[k]["vdn"] for k in UPS3_ORDER})
    print("spearman CE", rows_ce)
    return 0


if __name__ == "__main__":
    sys.exit(main())
