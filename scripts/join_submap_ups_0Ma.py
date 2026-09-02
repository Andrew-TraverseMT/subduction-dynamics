#!/usr/bin/env python3
"""Join 0 Ma Müller2019 trench kinematics to SubMap geodetic UPS.

Nearest trench sample to each SubMap transect (great-circle, lon wrapping),
drop joins > 200 km, recode UPS to a 3-class scheme, convert M56 mm/yr → cm/yr,
write the joined table, first-milestone replication figures, and join stats.

Run:
    /workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
        /workspace/subduction-test/scripts/join_submap_ups_0Ma.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path("/workspace/subduction-test")
KIN_CSV = ROOT / "results" / "trench_kinematics_0Ma_Muller2019.csv"
SUB_CSV = ROOT / "data" / "submap" / "submap_geodetic_UPS.csv"
OUT_CSV = ROOT / "results" / "joined_submap_muller2019_0Ma.csv"
STATS_MD = ROOT / "results" / "join_stats_0Ma.md"
FIG_VOP = ROOT / "results" / "fig_vOP_vs_UPS.png"
FIG_VT = ROOT / "results" / "fig_vT_vs_UPS.png"
FIG_CLASSIC = ROOT / "results" / "fig_classic_replication.png"
FIG_AGE = ROOT / "results" / "fig_age_vs_UPS.png"
FIG_SCATTER = ROOT / "results" / "fig_vOP_vs_submap_vdn.png"

EARTH_RADIUS_KM = 6371.009  # pygplates.Earth.mean_radius_in_kms
MAX_JOIN_KM = 200.0
SENTINEL = -999

# Full Geodetic_UPS display order used in box plots (Undefined excluded:
# those transects have missing M56 kinematics). This is the 9-class scheme
# minus Undefined.
UPS_ORDER = ["C", "C-SS", "SS-C", "SS", "SS-E", "E-SS", "E", "Neutral"]

# Compression (red) → strike-slip (gold/grey) → extension (blue) → Neutral.
UPS_COLORS = {
    "C": "#9B2226",
    "C-SS": "#BB3E03",
    "SS-C": "#CA6702",
    "SS": "#C9A227",
    "SS-E": "#94D2BD",
    "E-SS": "#0A9396",
    "E": "#005F73",
    "Neutral": "#6C757D",
}

# 3-class recode. Mixed strike-slip kept as SS so they can be shown
# separately; Undefined is not part of the C/N/E test.
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

M56_VEL_MM = ["M56_vc", "M56_vcn", "M56_vs", "M56_vsn", "M56_vd", "M56_vdn"]

KIN_COLS = [
    "lon",
    "lat",
    "subducting_plate_id",
    "trench_plate_id",
    "trench_normal_azimuth_deg",
    "arc_segment_length_km",
    "conv_mag_cm_yr",
    "conv_obliquity_deg",
    "conv_perp_cm_yr",
    "conv_parallel_cm_yr",
    "v_T_perp",
    "v_T_parallel",
    "v_OP_perp",
    "v_OP_parallel",
    "seafloor_age_Ma",
    "dist_nearest_trench_edge_km",
    "zone_id",
    "trench_zone_width_km",
    "upper_plate_nature",
]


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


def haversine_km(lon1, lat1, lon2, lat2):
    """Great-circle distance (km). Inputs in degrees; lon wrapping is inherent.

    Broadcasting: lon1/lat1 shape (n,), lon2/lat2 shape (m,) → (n, m).
    """
    lon1 = np.radians(np.asarray(lon1, dtype=float))[:, None]
    lat1 = np.radians(np.asarray(lat1, dtype=float))[:, None]
    lon2 = np.radians(np.asarray(lon2, dtype=float))[None, :]
    lat2 = np.radians(np.asarray(lat2, dtype=float))[None, :]
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    a = np.clip(a, 0.0, 1.0)
    return EARTH_RADIUS_KM * 2.0 * np.arcsin(np.sqrt(a))


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


def nearest_join(sub: pd.DataFrame, kin: pd.DataFrame) -> pd.DataFrame:
    """Nearest Müller trench sample to each SubMap (Lon, Lat)."""
    dist = haversine_km(
        sub["Lon"].to_numpy(),
        sub["Lat"].to_numpy(),
        kin["lon"].to_numpy(),
        kin["lat"].to_numpy(),
    )
    idx = dist.argmin(axis=1)
    dmin = dist[np.arange(dist.shape[0]), idx]

    kin_hit = kin.iloc[idx].reset_index(drop=True)
    kin_hit = kin_hit.rename(columns={"lon": "muller_lon", "lat": "muller_lat"})

    out = sub.reset_index(drop=True).copy()
    out["join_distance_km"] = dmin
    out["muller_row"] = kin.index.to_numpy()[idx]
    for col in kin_hit.columns:
        out[col] = kin_hit[col].to_numpy()
    return out


def box_strip(ax, df, ycol, order, colors, ylabel, title, annotate_n=True, zero_line=True):
    """Box plot + jittered strip. Empty classes still occupy the x slot."""
    rng = np.random.default_rng(20260902)
    data = []
    plotted_colors = []
    ns = []
    for cat in order:
        vals = df.loc[df["Geodetic_UPS"] == cat, ycol].to_numpy(dtype=float)
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


def spearman_safe(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 3:
        return np.nan, np.nan, n
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p), n


def main() -> int:
    _setup_style()
    RESULTS = ROOT / "results"
    RESULTS.mkdir(parents=True, exist_ok=True)

    print("Reading kinematics:", KIN_CSV)
    kin = pd.read_csv(KIN_CSV)
    print("Reading SubMap:", SUB_CSV)
    sub = pd.read_csv(SUB_CSV)

    n_sub = len(sub)
    if n_sub != 260:
        print(f"WARNING: SubMap n={n_sub} (expected 260)")
    missing_kin = [c for c in KIN_COLS if c not in kin.columns]
    if missing_kin:
        raise SystemExit(f"Kinematics CSV missing columns: {missing_kin}")
    for need in ["Short_Name", "Lon", "Lat", "Geodetic_UPS", "M56_vdn"]:
        if need not in sub.columns:
            raise SystemExit(f"SubMap CSV missing column: {need}")

    print(f"Kinematics n={len(kin)}  SubMap n={n_sub}")
    print(
        f"lon ranges  SubMap [{sub['Lon'].min():.1f}, {sub['Lon'].max():.1f}]  "
        f"Müller [{kin['lon'].min():.1f}, {kin['lon'].max():.1f}]"
    )

    joined_all = nearest_join(sub, kin)
    n_far = int((joined_all["join_distance_km"] > MAX_JOIN_KM).sum())
    n_ok = int((joined_all["join_distance_km"] <= MAX_JOIN_KM).sum())
    dropped = joined_all.loc[
        joined_all["join_distance_km"] > MAX_JOIN_KM,
        ["Short_Name", "Zone", "Trench_Name", "Lon", "Lat", "join_distance_km", "Geodetic_UPS"],
    ].copy()

    joined = joined_all.loc[joined_all["join_distance_km"] <= MAX_JOIN_KM].copy()
    joined = joined.reset_index(drop=True)

    joined["UPS_3"] = recode_ups3(joined["Geodetic_UPS"])

    # Sentinel → NaN, then mm/yr → cm/yr for plot-side-by-side with Müller.
    for col in M56_VEL_MM:
        if col not in joined.columns:
            continue
        mm = sentinel_to_nan(joined[col])
        joined[col] = mm  # keep original column, but without -999
        joined[f"{col}_cm_yr"] = mm / 10.0

    # Seismic_UPS is all sentinel; leave as-is (ignored).
    joined.to_csv(OUT_CSV, index=False, float_format="%.6f")
    print(f"Wrote {len(joined)} rows -> {OUT_CSV}")
    print(f"n matched (≤ {MAX_JOIN_KM:.0f} km): {n_ok} / {n_sub}")
    print(f"n dropped  (> {MAX_JOIN_KM:.0f} km): {n_far}")
    med_dist = float(np.median(joined["join_distance_km"]))
    print(f"median join distance: {med_dist:.3f} km")
    if n_far:
        print("Dropped transects:")
        print(dropped.to_string(index=False))

    # Plotting frame: matched rows, 8-class UPS (no Undefined).
    plot_df = joined.loc[joined["Geodetic_UPS"].isin(UPS_ORDER)].copy()

    # ------------------------------------------------------------------
    # Fig 1 — v_OP_perp vs Geodetic_UPS
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    box_strip(
        ax,
        plot_df,
        "v_OP_perp",
        UPS_ORDER,
        UPS_COLORS,
        ylabel=r"$v_{\mathrm{OP}\perp}$ (cm/yr), positive = upper-plate retreat",
        title=(
            "Müller et al. (2019) overriding-plate motion vs SubMap geodetic UPS\n"
            "0 Ma nearest-trench join (≤ 200 km)"
        ),
    )
    ax.set_xlabel("Geodetic UPS (SubMap)")
    fig.tight_layout()
    savefig(fig, FIG_VOP)

    # ------------------------------------------------------------------
    # Fig 2 — v_T_perp vs Geodetic_UPS
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    box_strip(
        ax,
        plot_df,
        "v_T_perp",
        UPS_ORDER,
        UPS_COLORS,
        ylabel=r"$v_{T\perp}$ (cm/yr), positive = oceanward trench retreat",
        title=(
            "Müller et al. (2019) trench-normal trench motion vs SubMap geodetic UPS\n"
            "0 Ma nearest-trench join (≤ 200 km)"
        ),
    )
    ax.set_xlabel("Geodetic UPS (SubMap)")
    fig.tight_layout()
    savefig(fig, FIG_VT)

    # ------------------------------------------------------------------
    # Fig 3 — classic replication: SubMap M56_vdn vs Müller v_OP
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 8.6), sharex=True)
    box_strip(
        axes[0],
        plot_df,
        "M56_vdn_cm_yr",
        UPS_ORDER,
        UPS_COLORS,
        ylabel=r"SubMap $v_{dn}$ (cm/yr)",
        title=(
            "Heuret-style replication with SubMap’s own M56 kinematics\n"
            r"shortening $\Rightarrow v_{dn}>0$; extension $\Rightarrow v_{dn}<0$"
        ),
    )
    box_strip(
        axes[1],
        plot_df,
        "v_OP_perp",
        UPS_ORDER,
        UPS_COLORS,
        ylabel=r"Müller $v_{\mathrm{OP}\perp}$ (cm/yr), positive = upper-plate retreat",
        title=(
            "Same UPS classes, Müller et al. (2019) overriding-plate motion\n"
            "Does reconstruction kinematics recover the SubMap pattern?"
        ),
    )
    axes[1].set_xlabel("Geodetic UPS (SubMap)")
    fig.suptitle(
        "Classic UPS replication, 0 Ma (nearest-trench join ≤ 200 km)",
        fontsize=13,
        y=1.01,
    )
    fig.tight_layout()
    savefig(fig, FIG_CLASSIC)

    # ------------------------------------------------------------------
    # Fig 4 — seafloor age vs UPS (drop NaN ages)
    # ------------------------------------------------------------------
    age_df = plot_df.loc[np.isfinite(plot_df["seafloor_age_Ma"])].copy()
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    box_strip(
        ax,
        age_df,
        "seafloor_age_Ma",
        UPS_ORDER,
        UPS_COLORS,
        ylabel="Seafloor age (Ma)",
        title=(
            "Seafloor age at trench vs SubMap geodetic UPS\n"
            "Müller 2019 age grid, 0 Ma; NaN ages dropped"
        ),
        zero_line=False,
    )
    ax.set_xlabel("Geodetic UPS (SubMap)")
    fig.tight_layout()
    savefig(fig, FIG_AGE)

    # ------------------------------------------------------------------
    # Fig 5 — scatter v_OP vs M56_vdn coloured by UPS_3
    # ------------------------------------------------------------------
    sc = joined.loc[
        joined["UPS_3"].isin(["C", "N", "E", "SS"])
        & np.isfinite(joined["v_OP_perp"])
        & np.isfinite(joined["M56_vdn_cm_yr"])
    ].copy()
    r_sc, p_sc, n_sc = spearman_safe(sc["v_OP_perp"], sc["M56_vdn_cm_yr"])

    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    for cls in UPS3_ORDER:
        subc = sc.loc[sc["UPS_3"] == cls]
        if subc.empty:
            continue
        ax.scatter(
            subc["M56_vdn_cm_yr"],
            subc["v_OP_perp"],
            s=32,
            c=UPS3_COLORS[cls],
            alpha=0.8,
            edgecolors="0.15",
            linewidths=0.35,
            label=f"{cls} (n={len(subc)})",
            zorder=3,
        )
    lims_x = ax.get_xlim()
    lims_y = ax.get_ylim()
    lo = min(lims_x[0], lims_y[0], -1)
    hi = max(lims_x[1], lims_y[1], 1)
    ax.plot([lo, hi], [lo, hi], color="0.65", lw=0.9, ls="--", zorder=1, label="1:1")
    ax.axhline(0.0, color="0.7", lw=0.6, zorder=1)
    ax.axvline(0.0, color="0.7", lw=0.6, zorder=1)
    ax.set_xlim(lims_x)
    ax.set_ylim(lims_y)
    ax.set_xlabel(r"SubMap M56 $v_{dn}$ (cm/yr)")
    ax.set_ylabel(r"Müller $v_{\mathrm{OP}\perp}$ (cm/yr)")
    ax.set_title(
        "Müller $v_{OP\\perp}$ vs SubMap $v_{dn}$ (0 Ma join)\n"
        "coloured by UPS$_3$  (C = C+C-SS; E = E+E-SS; N = Neutral; SS = mixed SS)"
    )
    ax.legend(frameon=False, loc="best", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(
        0.04,
        0.04,
        f"Spearman ρ = {fmt_rho(r_sc)}\n"
        f"p = {fmt_p(p_sc)}\nn = {n_sc}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.8", alpha=0.9),
    )
    fig.tight_layout()
    savefig(fig, FIG_SCATTER)

    # ------------------------------------------------------------------
    # Stats
    # ------------------------------------------------------------------
    defined = joined.loc[joined["Geodetic_UPS"] != "Undefined"].copy()

    r_op, p_op, n_op = spearman_safe(defined["v_OP_perp"], defined["M56_vdn_cm_yr"])
    r_t, p_t, n_t = spearman_safe(defined["v_T_perp"], defined["M56_vdn_cm_yr"])

    age_ok = defined.loc[
        np.isfinite(defined["seafloor_age_Ma"]) & np.isfinite(defined["M56_vdn_cm_yr"])
    ]
    r_age, p_age, n_age = spearman_safe(age_ok["seafloor_age_Ma"], age_ok["M56_vdn_cm_yr"])

    def med_by(col, cls):
        vals = joined.loc[joined["UPS_3"] == cls, col].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            return np.nan, 0
        return float(np.median(vals)), int(len(vals))

    meds = {}
    for cls in ["C", "N", "E", "SS"]:
        meds[cls] = {
            "v_OP": med_by("v_OP_perp", cls),
            "v_T": med_by("v_T_perp", cls),
        }

    vt_vs_vop = spearman_safe(joined["v_T_perp"], joined["v_OP_perp"])
    vt_plus_vop = joined["v_T_perp"] + joined["v_OP_perp"]
    med_sum = float(np.nanmedian(np.abs(vt_plus_vop)))

    ups_counts = (
        joined["Geodetic_UPS"].value_counts().reindex(UPS_ORDER + ["Undefined"]).fillna(0).astype(int)
    )
    ups3_counts = joined["UPS_3"].value_counts().reindex(["C", "N", "E", "SS", "Undefined"]).fillna(0).astype(int)

    dropped_lines = ""
    if n_far:
        rows = []
        for _, r in dropped.iterrows():
            rows.append(
                f"| {r['Short_Name']} | {r['Zone']} | {r['Trench_Name']} | "
                f"{r['Lon']:.2f} | {r['Lat']:.2f} | {r['join_distance_km']:.1f} | {r['Geodetic_UPS']} |"
            )
        dropped_lines = (
            "\nDropped transects (join > 200 km):\n\n"
            "| Short_Name | Zone | Trench_Name | Lon | Lat | dist_km | Geodetic_UPS |\n"
            "|---|---|---|---:|---:|---:|---|\n"
            + "\n".join(rows)
            + "\n"
        )

    def med_fmt(pair):
        med, n = pair
        if n == 0 or not np.isfinite(med):
            return f"NA (n=0)"
        return f"{med:.3f} cm/yr (n={n})"

    stats = f"""# 0 Ma SubMap UPS ↔ Müller 2019 trench-kinematics join

First-milestone look only. Values are computed from the joined table; nothing is invented.

- **Kinematics:** `{KIN_CSV}`
- **SubMap:** `{SUB_CSV}` (n={n_sub})
- **Joined table:** `{OUT_CSV}`
- **Join:** nearest Müller 2019 trench sample to each SubMap (Lon, Lat), great-circle distance with lon wrapping (SubMap lon is 0–360 in some zones; Müller is −180–180). Earth radius {EARTH_RADIUS_KM} km.
- **Cutoff:** drop joins farther than {MAX_JOIN_KM:.0f} km.

## Match rates

- n SubMap transects: **{n_sub}**
- n matched (≤ {MAX_JOIN_KM:.0f} km): **{n_ok}**
- n dropped (> {MAX_JOIN_KM:.0f} km): **{n_far}**
- median join distance (matched): **{med_dist:.2f} km**
- mean join distance (matched): **{float(joined['join_distance_km'].mean()):.2f} km**
- max join distance (matched): **{float(joined['join_distance_km'].max()):.2f} km**
{dropped_lines}
## UPS recode (`UPS_3`)

| Geodetic_UPS | UPS_3 | n matched |
|---|---|---:|
| C, C-SS | C | {int(ups3_counts['C'])} |
| Neutral | N | {int(ups3_counts['N'])} |
| E, E-SS | E | {int(ups3_counts['E'])} |
| SS, SS-C, SS-E | SS | {int(ups3_counts['SS'])} |
| Undefined | Undefined (out of C/N/E test) | {int(ups3_counts.get('Undefined', 0))} |

Geodetic_UPS counts among matched rows:

| class | n |
|---|---:|
"""
    for cat in UPS_ORDER + ["Undefined"]:
        stats += f"| {cat} | {int(ups_counts.get(cat, 0))} |\n"

    stats += f"""
M56 velocities were converted mm/yr → cm/yr (`*_cm_yr` columns) so they sit next to Müller cm/yr. SubMap sentinel −999 was set to NaN before conversion. `Seismic_UPS` is all −999 and is ignored.

## Spearman rank correlations (defined UPS only)

Defined UPS = `Geodetic_UPS != Undefined`, plus finite values of both series. Mixed strike-slip (UPS_3 = SS) **is included** here (it is a defined class); Undefined is not.

| pair | Spearman ρ | p | n |
|---|---:|---:|---:|
| `v_OP_perp` vs SubMap `M56_vdn` (cm/yr) | {fmt_rho(r_op)} | {fmt_p(p_op)} | {n_op} |
| `v_T_perp` vs SubMap `M56_vdn` (cm/yr) | {fmt_rho(r_t)} | {fmt_p(p_t)} | {n_t} |
| `seafloor_age_Ma` vs SubMap `M56_vdn` (cm/yr) | {fmt_rho(r_age)} | {fmt_p(p_age)} | {n_age} |

Sign conventions (do not mix them up):

- SubMap: shortening ⇒ 0 ≤ Ovd < 90°, `vdn > 0`; extension ⇒ 90° < Ovd ≤ 180°, `vdn < 0`.
- Müller `v_OP_perp > 0` = overriding plate moving **away from the trench** (upper-plate retreat).
- Müller `v_T_perp > 0` = **oceanward trench retreat**.

Because the signs differ, ρ(`v_OP_perp`, `M56_vdn`) need not be positive even if the two fields describe the same deformation sense (OP retreat ↔ extension ↔ `vdn < 0`). The observed sign of ρ is reported as computed; it is not a hypothesis verdict.

## Median kinematics by UPS_3

| UPS_3 | median `v_OP_perp` | median `v_T_perp` |
|---|---|---|
| C (C + C-SS) | {med_fmt(meds['C']['v_OP'])} | {med_fmt(meds['C']['v_T'])} |
| N (Neutral) | {med_fmt(meds['N']['v_OP'])} | {med_fmt(meds['N']['v_T'])} |
| E (E + E-SS) | {med_fmt(meds['E']['v_OP'])} | {med_fmt(meds['E']['v_T'])} |
| SS (SS + SS-C + SS-E), shown separately | {med_fmt(meds['SS']['v_OP'])} | {med_fmt(meds['SS']['v_T'])} |

C / N / E are the classes for a first H1/H2 look. SS is reported so it is not silently folded in.

## `v_T_perp` and `v_OP_perp` are not independent in Müller 2019

In this rigid/deforming-network model the trench line is typically attached to `trench_plate_id`, so **`v_T_perp ≈ −v_OP_perp`**. Among matched rows:

- Spearman ρ(`v_T_perp`, `v_OP_perp`) = {fmt_rho(vt_vs_vop[0])} (n={vt_vs_vop[2]})
- median |`v_T_perp` + `v_OP_perp`| = {med_sum:.4f} cm/yr

Correlations of `v_T_perp` vs `M56_vdn` and of `v_OP_perp` vs `M56_vdn` are therefore **not independent tests of H1 (trench rollback) vs H2 (upper-plate motion)**. They are approximately sign-flipped versions of the same kinematic field. Do not treat a “winner” from this pair.

## First-look caveat (do not overclaim)

This is a nearest-point geographic join of two independently built 0 Ma products, not a common segmentation. UPS class is SubMap’s geodetic classification (from their `Ovd` / `vd`); the top panel of `fig_classic_replication.png` therefore largely recovers SubMap’s own definition of UPS, and is a sanity check rather than an independent hypothesis test. The Müller panels ask whether reconstruction kinematics line up with that classification. Sample sizes in pure C and pure E are small. No hypothesis is declared a winner here.

## Figures

- `{FIG_VOP}`
- `{FIG_VT}`
- `{FIG_CLASSIC}`
- `{FIG_AGE}`
- `{FIG_SCATTER}`
"""

    STATS_MD.write_text(stats)
    print(f"Wrote {STATS_MD}")

    print("\n=== join summary ===")
    print(f"n matched: {n_ok}   n dropped: {n_far}   median dist: {med_dist:.2f} km")
    print(f"Spearman v_OP vs M56_vdn: rho={fmt_rho(r_op)}  p={fmt_p(p_op)}  n={n_op}")
    print(f"Spearman v_T  vs M56_vdn: rho={fmt_rho(r_t)}  p={fmt_p(p_t)}  n={n_t}")
    print(f"Spearman age  vs M56_vdn: rho={fmt_rho(r_age)}  p={fmt_p(p_age)}  n={n_age}")
    print(
        "median v_OP_perp  C={:.3f}  N={:.3f}  E={:.3f}".format(
            meds["C"]["v_OP"][0], meds["N"]["v_OP"][0], meds["E"]["v_OP"][0]
        )
    )
    print(
        "median v_T_perp   C={:.3f}  N={:.3f}  E={:.3f}".format(
            meds["C"]["v_T"][0], meds["N"]["v_T"][0], meds["E"]["v_T"][0]
        )
    )
    print("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
