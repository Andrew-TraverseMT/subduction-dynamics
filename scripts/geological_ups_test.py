#!/usr/bin/env python3
"""Join a geological UPS time series to Müller2019 trench kinematics and retest
the Sdrolias age gate and Schellart width/edge — without using bab_on or
dilatation as the outcome.

Geological UPS comes from published basin opening/cessation ages and Andean
shortening intervals in data/geological_ups/cenozoic_ups_intervals.csv.
It is independent of Müller et al. (2019) topologies.

Inference unit: named geological segment × time, not Müller zone_id, and not
tessellated points (those are spatially autocorrelated along the trench).

Does NOT re-extract trench kinematics.
Does NOT read bab_on or dilat_rate.

Run:
    /workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
        /workspace/subduction-test/scripts/geological_ups_test.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

ROOT = Path("/workspace/subduction-test")
INTERVALS_CSV = ROOT / "data/geological_ups/cenozoic_ups_intervals.csv"
KIN_CSV = ROOT / "results/trench_kinematics_0-60Ma_Muller2019.csv"
JOIN_CSV = ROOT / "results/joined_submap_muller2019_0Ma.csv"
RESULTS = ROOT / "results"
OUT_SUMMARY = RESULTS / "geological_ups_joined_summary.csv"
OUT_MD = RESULTS / "geological_ups_test.md"
FIG_AGE = RESULTS / "fig_geo_ups_age.png"
FIG_WIDTH = RESULTS / "fig_geo_ups_width_edge.png"
FIG_TIMELINE = RESULTS / "fig_geo_ups_timeline.png"

AGE_GATE_MA = 55.0
EDGE_NEAR_KM = 1200.0
CLASS_ORDER = ("E", "C", "N")
CLASS_COLORS = {"E": "#005F73", "C": "#9B2226", "N": "#6C757D"}
# Prefer E over C over N when two intervals hit the same sample.
CLASS_PRIORITY = {"E": 0, "C": 1, "N": 2}

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


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Geological UPS join + age/width tests.")
    p.add_argument("--intervals", type=Path, default=INTERVALS_CSV)
    p.add_argument("--kin-csv", type=Path, default=KIN_CSV)
    p.add_argument("--join-csv", type=Path, default=JOIN_CSV)
    p.add_argument("--out-summary", type=Path, default=OUT_SUMMARY)
    p.add_argument("--out-md", type=Path, default=OUT_MD)
    return p.parse_args(argv)


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


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p <= 0:
        return "<1e-16"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_n(x) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.1f}"


def in_bbox(lon, lat, lon_min, lon_max, lat_min, lat_max):
    lat_ok = (lat >= lat_min) & (lat <= lat_max)
    if lon_min <= lon_max:
        lon_ok = (lon >= lon_min) & (lon <= lon_max)
    else:
        lon_ok = (lon >= lon_min) | (lon <= lon_max)
    return lat_ok & lon_ok


def load_intervals(path: Path) -> pd.DataFrame:
    iv = pd.read_csv(path)
    required = {
        "zone_name",
        "segment",
        "lon_min",
        "lon_max",
        "lat_min",
        "lat_max",
        "ups_class",
        "t_start_Ma",
        "t_end_Ma",
    }
    missing = required - set(iv.columns)
    if missing:
        raise ValueError(f"intervals CSV missing columns: {sorted(missing)}")
    iv = iv.copy()
    iv["ups_class"] = iv["ups_class"].str.strip().str.upper()
    bad = ~iv["ups_class"].isin(CLASS_ORDER)
    if bad.any():
        raise ValueError(f"bad ups_class values: {iv.loc[bad, 'ups_class'].unique()}")
    if (iv["t_start_Ma"] < iv["t_end_Ma"]).any():
        raise ValueError("t_start_Ma must be older (larger) than t_end_Ma")
    iv["priority"] = iv["ups_class"].map(CLASS_PRIORITY)
    iv = iv.sort_values(["priority", "t_start_Ma"], ascending=[True, False]).reset_index(
        drop=True
    )
    return iv


def assign_geo_ups(kin: pd.DataFrame, iv: pd.DataFrame) -> pd.DataFrame:
    """Assign geological UPS / segment / zone_name to each trench sample."""
    n = len(kin)
    geo = np.array(["Unknown"] * n, dtype=object)
    segment = np.array([""] * n, dtype=object)
    zone = np.array([""] * n, dtype=object)
    filled = np.zeros(n, dtype=bool)

    lon = kin["lon"].to_numpy(dtype=float)
    lat = kin["lat"].to_numpy(dtype=float)
    time = kin["time_Ma"].to_numpy(dtype=float)

    for row in iv.itertuples(index=False):
        time_ok = (time >= float(row.t_end_Ma)) & (time <= float(row.t_start_Ma))
        box_ok = in_bbox(
            lon,
            lat,
            float(row.lon_min),
            float(row.lon_max),
            float(row.lat_min),
            float(row.lat_max),
        )
        hit = time_ok & box_ok & ~filled
        if not np.any(hit):
            continue
        geo[hit] = row.ups_class
        segment[hit] = row.segment
        zone[hit] = row.zone_name
        filled[hit] = True

    out = kin.copy()
    out["geo_ups"] = geo
    out["geo_segment"] = segment
    out["geo_zone"] = zone
    out["geo_classified"] = filled
    return out


def aggregate_segment_time(df: pd.DataFrame) -> pd.DataFrame:
    """Named-segment × time blocks. One row per classified segment-time."""
    sub = df.loc[df["geo_classified"]].copy()
    if sub.empty:
        return pd.DataFrame()

    def _mode(s: pd.Series):
        v = s.dropna()
        if v.empty:
            return ""
        return v.value_counts().index[0]

    g = sub.groupby(["time_Ma", "geo_segment"], sort=True)
    rows = []
    for (t, seg), part in g:
        ups = _mode(part["geo_ups"])
        zone = _mode(part["geo_zone"])
        age = part["seafloor_age_Ma"].to_numpy(dtype=float)
        rows.append(
            {
                "time_Ma": float(t),
                "segment": seg,
                "geo_ups": ups,
                "geo_zone": zone,
                "n_samples": int(len(part)),
                "n_age": int(np.isfinite(age).sum()),
                "median_age_Ma": float(np.nanmedian(age)) if np.isfinite(age).any() else np.nan,
                "median_width_km": float(part["trench_zone_width_km"].median()),
                "median_edge_km": float(part["dist_nearest_trench_edge_km"].median()),
                "median_lon": float(part["lon"].median()),
                "median_lat": float(part["lat"].median()),
            }
        )
    out = pd.DataFrame(rows)
    out["is_E"] = (out["geo_ups"] == "E").astype(int)
    out["is_C"] = (out["geo_ups"] == "C").astype(int)
    out["age_ge_55"] = np.where(
        np.isfinite(out["median_age_Ma"]),
        (out["median_age_Ma"] >= AGE_GATE_MA).astype(float),
        np.nan,
    )
    out["edge_near"] = (out["median_edge_km"] < EDGE_NEAR_KM).astype(int)
    return out


def mwu(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return {
            "n_a": int(len(a)),
            "n_b": int(len(b)),
            "median_a": float(np.median(a)) if len(a) else np.nan,
            "median_b": float(np.median(b)) if len(b) else np.nan,
            "U": np.nan,
            "p": np.nan,
        }
    U, p = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "U": float(U),
        "p": float(p),
    }


def gate_table(flag_E, age, classified_ok):
    """P(E | age>=55) vs P(E | age<55) among rows that are classified and have age."""
    flag_E = np.asarray(flag_E, dtype=bool)
    age = np.asarray(age, dtype=float)
    classified_ok = np.asarray(classified_ok, dtype=bool)
    ok = classified_ok & np.isfinite(age)
    old = ok & (age >= AGE_GATE_MA)
    young = ok & (age < AGE_GATE_MA)
    n_old = int(old.sum())
    n_young = int(young.sum())
    n_e_old = int((flag_E & old).sum())
    n_e_young = int((flag_E & young).sum())
    p_old = (n_e_old / n_old) if n_old else np.nan
    p_young = (n_e_young / n_young) if n_young else np.nan
    # 2x2: rows age>=55 / age<55; cols E / not-E
    a, b = n_e_old, n_old - n_e_old
    c, d = n_e_young, n_young - n_e_young
    if min(a, b, c, d) >= 0 and (a + b) > 0 and (c + d) > 0 and (a + c) > 0 and (b + d) > 0:
        table = np.array([[a, b], [c, d]], dtype=int)
        _, p_fish = fisher_exact(table, alternative="two-sided")
        # Haldane-Anscombe if a zero cell
        aa, bb, cc, dd = (x + 0.5 if min(a, b, c, d) == 0 else x for x in (a, b, c, d))
        odds = (aa * dd) / (bb * cc)
    else:
        p_fish = np.nan
        odds = np.nan
    return {
        "n_old": n_old,
        "n_young": n_young,
        "n_e_old": n_e_old,
        "n_e_young": n_e_young,
        "p_old": p_old,
        "p_young": p_young,
        "odds": odds,
        "p_fisher": p_fish,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
    }


def spearman_flag(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return np.nan, np.nan, int(m.sum())
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def submap_sanity(df0: pd.DataFrame, join_path: Path, iv: pd.DataFrame) -> dict:
    if not join_path.exists():
        return {"ok": False, "reason": f"missing {join_path}"}
    j = pd.read_csv(join_path)
    lon_col = "muller_lon" if "muller_lon" in j.columns else "Lon"
    lat_col = "muller_lat" if "muller_lat" in j.columns else "Lat"
    if "UPS_3" not in j.columns:
        j["UPS_3"] = j["Geodetic_UPS"].map(UPS3_MAP)
    # Assign from 0 Ma intervals using Müller coordinates on the join table.
    geo = np.array(["Unknown"] * len(j), dtype=object)
    seg = np.array([""] * len(j), dtype=object)
    filled = np.zeros(len(j), dtype=bool)
    lon = j[lon_col].to_numpy(dtype=float)
    lat = j[lat_col].to_numpy(dtype=float)
    for row in iv.itertuples(index=False):
        if float(row.t_end_Ma) > 0.0 or float(row.t_start_Ma) < 0.0:
            # interval must include 0 Ma
            if not (float(row.t_end_Ma) <= 0.0 <= float(row.t_start_Ma)):
                continue
        else:
            if not (float(row.t_end_Ma) <= 0.0 <= float(row.t_start_Ma)):
                continue
        box_ok = in_bbox(
            lon,
            lat,
            float(row.lon_min),
            float(row.lon_max),
            float(row.lat_min),
            float(row.lat_max),
        )
        hit = box_ok & ~filled
        if not np.any(hit):
            continue
        geo[hit] = row.ups_class
        seg[hit] = row.segment
        filled[hit] = True
    j = j.copy()
    j["geo_ups"] = geo
    j["geo_segment"] = seg
    j["geo_classified"] = filled
    ov = j.loc[j["geo_classified"]].copy()
    # Relaxed agreement: E↔E, C↔C, N↔{N,SS}
    def agree_strict(g, u):
        return g == u

    def agree_relaxed(g, u):
        if g == "E":
            return u == "E"
        if g == "C":
            return u == "C"
        if g == "N":
            return u in ("N", "SS")
        return False

    ov["agree_strict"] = [
        agree_strict(g, u) for g, u in zip(ov["geo_ups"], ov["UPS_3"])
    ]
    ov["agree_relaxed"] = [
        agree_relaxed(g, u) for g, u in zip(ov["geo_ups"], ov["UPS_3"])
    ]
    # drop Undefined from rates
    defined = ov.loc[ov["UPS_3"].isin(["C", "N", "E", "SS"])].copy()
    ct = pd.crosstab(defined["geo_ups"], defined["UPS_3"]).reindex(
        index=["E", "C", "N"], columns=["E", "C", "N", "SS"], fill_value=0
    )
    by_seg = (
        defined.groupby("geo_segment")
        .agg(
            n=("UPS_3", "size"),
            geo=("geo_ups", lambda s: s.mode().iloc[0] if len(s) else ""),
            submap=("UPS_3", lambda s: s.value_counts().to_dict()),
            n_agree_relaxed=("agree_relaxed", "sum"),
        )
        .reset_index()
    )
    return {
        "ok": True,
        "n_join": int(len(j)),
        "n_overlap": int(len(ov)),
        "n_defined": int(len(defined)),
        "n_agree_strict": int(defined["agree_strict"].sum()) if len(defined) else 0,
        "n_agree_relaxed": int(defined["agree_relaxed"].sum()) if len(defined) else 0,
        "ct": ct,
        "by_seg": by_seg,
        "overlap": defined,
    }


def fig_age(df: pd.DataFrame, blocks: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))

    cls = df.loc[df["geo_classified"] & np.isfinite(df["seafloor_age_Ma"])]
    ax = axes[0]
    bins = np.arange(0, 161, 5)
    for lab, color, alpha in (("C", CLASS_COLORS["C"], 0.55), ("E", CLASS_COLORS["E"], 0.55)):
        vals = cls.loc[cls["geo_ups"] == lab, "seafloor_age_Ma"].to_numpy()
        if len(vals) == 0:
            continue
        ax.hist(
            vals,
            bins=bins,
            density=True,
            color=color,
            alpha=alpha,
            label=f"{lab} n={len(vals)} med={np.median(vals):.0f} Ma",
            histtype="stepfilled",
            edgecolor=color,
        )
    ax.axvline(AGE_GATE_MA, color="#9B2226", ls="--", lw=1.1, label="55 Ma gate")
    ax.set_xlim(0, 160)
    ax.set_xlabel("Seafloor age at trench (Ma)")
    ax.set_ylabel("Density")
    ax.set_title("Classified samples, pooled 0–60 Ma")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    # One polyline per segment so E is not a zigzag across unrelated trenches.
    # Tethyan Hellenic/Calabria ages (~200–300 Ma) omitted (off-scale).
    tethyan = {"Hellenic", "Calabria"}
    vis = blocks.loc[
        blocks["geo_ups"].isin(["E", "C"]) & ~blocks["segment"].isin(tethyan)
    ].copy()
    # Distinct E segments get thin teal lines; C (Andes, Japan) in red.
    e_segs = sorted(vis.loc[vis["geo_ups"] == "E", "segment"].unique())
    cmap_e = plt.cm.Blues(np.linspace(0.45, 0.9, max(len(e_segs), 2)))
    for i, seg in enumerate(e_segs):
        b = vis.loc[(vis["segment"] == seg) & (vis["geo_ups"] == "E")].sort_values(
            "time_Ma"
        )
        ax.plot(
            b["time_Ma"],
            b["median_age_Ma"],
            color=cmap_e[i],
            marker="o",
            ms=2.8,
            lw=1.1,
            label=seg if seg in ("Mariana", "Tonga", "Scotia", "Andaman", "NewHebrides") else None,
        )
    for seg, ls in (("Andes", "-"), ("Japan", "--")):
        b = vis.loc[(vis["segment"] == seg) & (vis["geo_ups"] == "C")].sort_values(
            "time_Ma"
        )
        if b.empty:
            continue
        ax.plot(
            b["time_Ma"],
            b["median_age_Ma"],
            color=CLASS_COLORS["C"],
            marker="s",
            ms=3.0,
            lw=1.6,
            ls=ls,
            label=f"{seg} C",
        )
    ax.axhline(AGE_GATE_MA, color="#9B2226", ls="--", lw=1.0)
    ax.set_xlim(60, 0)
    ax.set_ylim(0, 160)
    ax.set_xlabel("Reconstruction time (Ma)")
    ax.set_ylabel("Median seafloor age (Ma)")
    ax.set_title("Segment × time (Hellenic/Calabria omitted)")
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper left")

    fig.suptitle(
        "Geological UPS — seafloor age at trench (not Müller back-arc MORs)",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def fig_width_edge(blocks: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.4))
    ce = blocks.loc[blocks["geo_ups"].isin(["E", "C"])].copy()

    def _box(ax, col, ylabel, title):
        data, labels, colors = [], [], []
        for lab in ("E", "C"):
            v = ce.loc[ce["geo_ups"] == lab, col].to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            if len(v) == 0:
                continue
            data.append(v)
            labels.append(f"{lab}\nn={len(v)}")
            colors.append(CLASS_COLORS[lab])
        if not data:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            return
        bp = ax.boxplot(
            data,
            tick_labels=labels,
            patch_artist=True,
            medianprops={"color": "0.1", "lw": 1.4},
            whiskerprops={"color": "0.3"},
            capprops={"color": "0.3"},
            flierprops={"marker": "o", "ms": 3, "mec": "0.4", "mfc": "none"},
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
            patch.set_edgecolor(c)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    _box(
        axes[0],
        "median_width_km",
        "Trench-zone width (km)",
        "Width (Schellart: E narrower?)",
    )
    _box(
        axes[1],
        "median_edge_km",
        "Distance to trench edge (km)",
        "Edge (Schellart: E nearer edge?)",
    )

    ax = axes[2]
    for lab in ("E", "C", "N"):
        b = blocks.loc[blocks["geo_ups"] == lab]
        if b.empty:
            continue
        ax.scatter(
            b["median_width_km"],
            b["median_edge_km"],
            c=CLASS_COLORS[lab],
            s=28,
            alpha=0.85,
            label=lab,
            edgecolors="0.2",
            linewidths=0.3,
        )
    ax.axhline(EDGE_NEAR_KM, color="0.5", ls=":", lw=1.0)
    ax.set_xlabel("Trench-zone width (km)")
    ax.set_ylabel("Median dist. to edge (km)")
    ax.set_title("Segment × time")
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle(
        "Geological UPS — Schellart width / edge (segment × time)",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def fig_timeline(iv: pd.DataFrame, path: Path) -> None:
    plot = iv.copy()
    # display order: compression / Andes first, then Pacific BABs north to south-ish
    seg_order = [
        "Andes",
        "Cascadia",
        "Japan",
        "Izu-Bonin",
        "Mariana",
        "Ryukyu",
        "Andaman",
        "Manus",
        "NewHebrides",
        "Tonga",
        "Kermadec",
        "Scotia",
        "Hellenic",
        "Calabria",
    ]
    plot["seg_rank"] = plot["segment"].map(
        {s: i for i, s in enumerate(seg_order)}
    ).fillna(99)
    plot = plot.sort_values(["seg_rank", "t_start_Ma"], ascending=[True, False])
    plot["y"] = np.arange(len(plot))[::-1]

    fig, ax = plt.subplots(figsize=(11.4, 7.6))
    for r in plot.itertuples(index=False):
        left = float(r.t_end_Ma)
        width = float(r.t_start_Ma) - float(r.t_end_Ma)
        ax.barh(
            r.y,
            width,
            left=left,
            height=0.72,
            color=CLASS_COLORS[r.ups_class],
            edgecolor="0.15",
            linewidth=0.4,
            alpha=0.9,
        )
    ax.set_yticks(plot["y"].to_numpy())
    labels = [
        f"{z.replace('_', ' ')}  ({s})"
        for z, s in zip(plot["zone_name"], plot["segment"])
    ]
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(62, -1)
    ax.set_xlabel("Age (Ma)")
    ax.set_title("Geological UPS intervals (published ages; not Müller topologies)")
    handles = [
        mpatches.Patch(
            color=CLASS_COLORS[k],
            label={"E": "E  extension / rifting", "C": "C  shortening", "N": "N  no back-arc / neutral"}[k],
        )
        for k in CLASS_ORDER
    ]
    ax.legend(handles=handles, loc="upper right", frameon=True, fontsize=9, fancybox=False, framealpha=0.92)
    ax.axvline(0, color="0.7", lw=0.6)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def md_table(rows, headers) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        lines.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(lines)


def write_report(
    path: Path,
    iv: pd.DataFrame,
    df: pd.DataFrame,
    blocks: pd.DataFrame,
    stats: dict,
    submap: dict,
) -> None:
    n_iv = len(iv)
    n_samp = int(len(df))
    n_cls = int(df["geo_classified"].sum())
    n_unk = n_samp - n_cls
    n_blocks = int(len(blocks))
    vc = df.loc[df["geo_classified"], "geo_ups"].value_counts()
    bvc = blocks["geo_ups"].value_counts() if len(blocks) else pd.Series(dtype=int)

    age_s = stats["age_samples"]
    age_b = stats["age_blocks"]
    gate_s = stats["gate_samples"]
    gate_b = stats["gate_blocks"]
    gate_s_ec = stats["gate_samples_ec"]
    gate_b_ec = stats["gate_blocks_ec"]
    w_b = stats["width_blocks"]
    e_b = stats["edge_blocks"]

    lines = []
    lines.append("# Geological UPS test — Cenozoic, independent of Müller back-arc MORs")
    lines.append("")
    lines.append("Built by `scripts/geological_ups_test.py`.")
    lines.append("")
    lines.append("**Outcome:** published geological / marine-geophysical UPS intervals")
    lines.append("(`data/geological_ups/cenozoic_ups_intervals.csv`).")
    lines.append("**Predictors:** Müller et al. (2019) trench kinematics already in")
    lines.append("`results/trench_kinematics_0-60Ma_Muller2019.csv` (seafloor age, trench-zone")
    lines.append("width, distance to nearest trench edge). Trench kinematics were **not**")
    lines.append("re-extracted.")
    lines.append("")
    lines.append("This is **not** `bab_on` and **not** deforming-network dilatation. Those")
    lines.append("proxies are reconstruction-internal (the topologies already know where they")
    lines.append("opened back-arc basins). This test asks whether an *independent* geological")
    lines.append("UPS series still shows a Sdrolias-style age gate or a Schellart-style")
    lines.append("narrow / near-edge preference.")
    lines.append("")
    lines.append("## Sparse by construction — read this first")
    lines.append("")
    lines.append(f"- **n intervals:** {n_iv} (E={int((iv.ups_class=='E').sum())}, "
                 f"C={int((iv.ups_class=='C').sum())}, N={int((iv.ups_class=='N').sum())})")
    lines.append(f"- **n named segments:** {iv['segment'].nunique()}")
    lines.append(f"- **n trench sample-times:** {n_samp} (61 times, 0–60 Ma, 1 Myr)")
    lines.append(f"- **n classified sample-times:** {n_cls} ({100*n_cls/n_samp:.1f}%)")
    lines.append(f"- **n Unknown:** {n_unk}")
    lines.append(f"- **n segment × time blocks (classified):** {n_blocks}")
    lines.append(f"- sample-level classified: E={int(vc.get('E',0))}, C={int(vc.get('C',0))}, N={int(vc.get('N',0))}")
    lines.append(f"- block-level: E={int(bvc.get('E',0))}, C={int(bvc.get('C',0))}, N={int(bvc.get('N',0))}")
    lines.append("")
    lines.append("Hand intervals cover well-dated back-arc basins plus Central Andean")
    lines.append("shortening and Cascadia (no back-arc). Most of the global trench catalog")
    lines.append("is `Unknown` and is **not** treated as C. Ages have published ranges;")
    lines.append("midpoints are used. Where literature conflicts, the CSV `notes` field")
    lines.append("records the conflict. **Do not read a p-value on 76k points as independent")
    lines.append("information** — neighbouring tessellated samples along a trench are not")
    lines.append("independent. The inference unit is **named segment × time**.")
    lines.append("")
    lines.append("Cenozoic C is dominated by the Central Andes box (one wide, long-lived")
    lines.append("shortening orogen). Cenozoic E is several short-lived, mostly western-Pacific")
    lines.append("and Mediterranean basins. That imbalance is geological, not a sampling trick,")
    lines.append("but it means E vs C is largely “Pacific BABs vs Andes”.")
    lines.append("")
    lines.append("## How samples are classified")
    lines.append("")
    lines.append("A trench sample is classified if its longitude/latitude at that reconstruction")
    lines.append("time falls inside an active interval bbox (`t_end_Ma ≤ time ≤ t_start_Ma`).")
    lines.append("Longitude wrap (`lon_min > lon_max`) is used for Tonga and Kermadec.")
    lines.append("If two intervals hit the same sample, **E then C then N**. Otherwise Unknown.")
    lines.append("Bboxes follow the trench in this catalog, not the basin floor; false-positive")
    lines.append("risk is documented in `data/geological_ups/README.md`.")
    lines.append("")
    lines.append("Aegean, Okinawa, Havre, and Tyrrhenian count as **E during rifting**, even")
    lines.append("when there is no MOR. Cascadia is **N** (no back-arc; margin-parallel")
    lines.append("shortening in Washington is not Andean C).")
    lines.append("")
    lines.append("## Age — Sdrolias gate")
    lines.append("")
    lines.append("Sdrolias & Müller (2006) associated back-arc opening with subducting")
    lines.append(f"lithosphere older than ~{AGE_GATE_MA:.0f} Ma. Test among classified rows")
    lines.append("with finite `seafloor_age_Ma`. Hellenic / Calabria ages in this catalog are")
    lines.append("often Mesozoic (Tethyan) or NaN; they still enter E when classified.")
    lines.append("Dropping those two segments (sensitivity below) makes E **younger**, not older.")
    lines.append("")
    lines.append("The Cenozoic geological catalog does **not** reproduce the present-day")
    lines.append("old-slab → back-arc association. West Philippine / early IBM extension")
    lines.append("(55–30 Ma) sits on ~17–40 Ma Pacific lithosphere in these age grids,")
    lines.append("while Central Andean shortening sits on older Nazca. Gate OR is < 1")
    lines.append("(younger slabs more often E among classified blocks).")
    lines.append("")
    lines.append("### Segment × time (preferred)")
    lines.append("")
    lines.append(
        f"- Median age **E vs C**: {fmt_n(age_b['median_a'])} Ma (n={age_b['n_a']}) vs "
        f"{fmt_n(age_b['median_b'])} Ma (n={age_b['n_b']}); Mann–Whitney p={fmt_p(age_b['p'])}"
    )
    lines.append(
        f"- Median age **E vs not-E** (C+N): {fmt_n(stats['age_blocks_e_not']['median_a'])} Ma "
        f"(n={stats['age_blocks_e_not']['n_a']}) vs "
        f"{fmt_n(stats['age_blocks_e_not']['median_b'])} Ma "
        f"(n={stats['age_blocks_e_not']['n_b']}); p={fmt_p(stats['age_blocks_e_not']['p'])}"
    )
    lines.append(
        f"- Gate among classified blocks: P(E | age≥55) = "
        f"{gate_b['p_old'] if np.isfinite(gate_b['p_old']) else np.nan:.3f} "
        f"({gate_b['n_e_old']}/{gate_b['n_old']}) vs P(E | age<55) = "
        f"{gate_b['p_young'] if np.isfinite(gate_b['p_young']) else np.nan:.3f} "
        f"({gate_b['n_e_young']}/{gate_b['n_young']}); odds ratio = "
        f"{fmt_n(gate_b['odds'])}; Fisher p={fmt_p(gate_b['p_fisher'])}"
    )
    lines.append(
        f"- Gate **E vs C only** (drop N): P(E | age≥55) = "
        f"{gate_b_ec['p_old'] if np.isfinite(gate_b_ec['p_old']) else np.nan:.3f} "
        f"({gate_b_ec['n_e_old']}/{gate_b_ec['n_old']}) vs P(E | age<55) = "
        f"{gate_b_ec['p_young'] if np.isfinite(gate_b_ec['p_young']) else np.nan:.3f} "
        f"({gate_b_ec['n_e_young']}/{gate_b_ec['n_young']}); OR = "
        f"{fmt_n(gate_b_ec['odds'])}; p={fmt_p(gate_b_ec['p_fisher'])}"
    )
    age_nt = stats["age_blocks_notethys"]
    gate_nt = stats["gate_blocks_notethys"]
    lines.append(
        f"- Sensitivity **drop Hellenic+Calabria** (Tethyan ages): median age E vs C "
        f"{fmt_n(age_nt['median_a'])} vs {fmt_n(age_nt['median_b'])} Ma "
        f"(n={age_nt['n_a']}/{age_nt['n_b']}); "
        f"P(E|≥55)={gate_nt['p_old']:.3f} ({gate_nt['n_e_old']}/{gate_nt['n_old']}) vs "
        f"P(E|<55)={gate_nt['p_young']:.3f} ({gate_nt['n_e_young']}/{gate_nt['n_young']}); "
        f"OR={fmt_n(gate_nt['odds'])}; p={fmt_p(gate_nt['p_fisher'])}"
    )
    lines.append("")
    lines.append("### Sample-level (autocorrelated; supplementary)")
    lines.append("")
    lines.append(
        f"- Median age **E vs C**: {fmt_n(age_s['median_a'])} Ma (n={age_s['n_a']}) vs "
        f"{fmt_n(age_s['median_b'])} Ma (n={age_s['n_b']}); p={fmt_p(age_s['p'])}"
    )
    lines.append(
        f"- Gate: P(E | age≥55) = "
        f"{gate_s['p_old'] if np.isfinite(gate_s['p_old']) else np.nan:.3f} "
        f"({gate_s['n_e_old']}/{gate_s['n_old']}) vs P(E | age<55) = "
        f"{gate_s['p_young'] if np.isfinite(gate_s['p_young']) else np.nan:.3f} "
        f"({gate_s['n_e_young']}/{gate_s['n_young']}); OR = {fmt_n(gate_s['odds'])}; "
        f"p={fmt_p(gate_s['p_fisher'])}"
    )
    lines.append(
        f"- Gate E vs C only: OR = {fmt_n(gate_s_ec['odds'])}; "
        f"p={fmt_p(gate_s_ec['p_fisher'])} "
        f"(E|old {gate_s_ec['n_e_old']}/{gate_s_ec['n_old']}; "
        f"E|young {gate_s_ec['n_e_young']}/{gate_s_ec['n_young']})"
    )
    lines.append("")
    lines.append("Sample-level p-values will look decisive because n is large and points along")
    lines.append("a trench are copies of each other. Believe the block-level n.")
    lines.append("")
    lines.append("## Width and edge — Schellart direction")
    lines.append("")
    lines.append("Schellart et al. (2007, *Nature*): rollback / back-arc opening concentrate")
    lines.append("on **narrow** slabs and **near a lateral slab edge**; interiors of wide slabs")
    lines.append("are more compressive. Operationally, at segment × time:")
    lines.append("")
    lines.append("- narrower `trench_zone_width_km` → more E → Spearman ρ(width, is_E) **negative**")
    lines.append("- smaller `dist_nearest_trench_edge_km` → more E → ρ(edge, is_E) **negative**")
    lines.append(f"- near-edge split at {EDGE_NEAR_KM:.0f} km")
    lines.append("")
    lines.append("`trench_zone_width_km` is the along-trench length of a connected polyline")
    lines.append("sharing a subducting-plate ID in Müller2019, **not** a tomographic slab width.")
    lines.append("Median edge distance is mechanically related to width (a uniform sample on a")
    lines.append("zone of length W has median distance-to-edge ≈ W/4).")
    lines.append("")
    lines.append(
        f"- Median width **E vs C**: {fmt_n(w_b['median_a'])} km (n={w_b['n_a']}) vs "
        f"{fmt_n(w_b['median_b'])} km (n={w_b['n_b']}); p={fmt_p(w_b['p'])}"
    )
    lines.append(
        f"- Median edge **E vs C**: {fmt_n(e_b['median_a'])} km (n={e_b['n_a']}) vs "
        f"{fmt_n(e_b['median_b'])} km (n={e_b['n_b']}); p={fmt_p(e_b['p'])}"
    )
    r_w, p_w, n_w = stats["spearman_width"]
    r_e, p_e, n_e = stats["spearman_edge"]
    lines.append(
        f"- Spearman ρ(width, is_E) = {r_w:.3f} (p={fmt_p(p_w)}, n={n_w}) "
        f"— Schellart sign is negative"
    )
    lines.append(
        f"- Spearman ρ(edge, is_E) = {r_e:.3f} (p={fmt_p(p_e)}, n={n_e}) "
        f"— Schellart sign is negative"
    )
    ne = stats["near_edge"]
    lines.append(
        f"- P(E | median edge < {EDGE_NEAR_KM:.0f} km) = {ne['p_near']:.3f} "
        f"({ne['n_e_near']}/{ne['n_near']}) vs P(E | ≥{EDGE_NEAR_KM:.0f}) = "
        f"{ne['p_far']:.3f} ({ne['n_e_far']}/{ne['n_far']})"
    )
    lines.append("")
    lines.append("Andes is a wide C orogen; Scotia / Manus / Andaman are short E segments.")
    lines.append("A width contrast in the E-vs-C medians is **consistent with** Schellart")
    lines.append("but is also the geometry of boxing those systems. Spearman ρ(width, is_E)")
    lines.append("across *all* classified blocks (including N and IBM) is near zero because")
    lines.append("Müller `trench_zone_width_km` often glues IBM–Japan–Kuril–Aleutian into a")
    lines.append("10,000–20,000 km polyline — that is **not** a tomographic slab width and")
    lines.append("is not the Schellart 2007 width of the Mariana or Tonga slab. Scotia")
    lines.append("(~800 km) and Manus (~1,100 km) are the honest narrow E end-member.")
    lines.append("Do not over-claim.")
    lines.append("")
    lines.append("## 0 Ma sanity vs SubMap geodetic UPS")
    lines.append("")
    if not submap.get("ok"):
        lines.append(f"SubMap join file not available ({submap.get('reason','')}).")
    else:
        n_def = submap["n_defined"]
        a_s = submap["n_agree_strict"]
        a_r = submap["n_agree_relaxed"]
        lines.append(
            f"Geological UPS at 0 Ma assigned to SubMap transects using Müller coordinates "
            f"in `{JOIN_CSV.name}`. Overlap classified: **{submap['n_overlap']}** of "
            f"{submap['n_join']} joined transects; **{n_def}** with defined `UPS_3` "
            f"(not Undefined)."
        )
        lines.append("")
        lines.append("Agreement:")
        lines.append(f"- **strict** (E↔E, C↔C, N↔N): {a_s}/{n_def} = {a_s/n_def:.2f}" if n_def else "- none")
        lines.append(
            f"- **relaxed** (E↔E, C↔C, N↔{{N,SS}}): {a_r}/{n_def} = {a_r/n_def:.2f}"
            if n_def
            else "- none"
        )
        lines.append("")
        lines.append("Relaxed N↔SS is the intended Cascadia rule (Wells et al. 1998;")
        lines.append("SubMap Cascadia is SS, we labelled N). Andaman is an expected")
        lines.append("**disagreement**: geological E (oblique Central Andaman spreading) vs")
        lines.append("SubMap SS.")
        lines.append("")
        lines.append("Crosstab geological UPS × SubMap `UPS_3` (defined only):")
        lines.append("")
        ct = submap["ct"]
        lines.append("| geo\\\\SubMap | E | C | N | SS |")
        lines.append("|---|---:|---:|---:|---:|")
        for idx in ct.index:
            vals = " | ".join(str(int(ct.loc[idx, c])) for c in ct.columns)
            lines.append(f"| {idx} | {vals} |")
        lines.append("")
        lines.append("| segment | n | geo | SubMap UPS_3 | n relaxed-agree |")
        lines.append("|---|---:|---|---|---:|")
        for r in submap["by_seg"].itertuples(index=False):
            lines.append(
                f"| {r.geo_segment} | {int(r.n)} | {r.geo} | {r.submap} | {int(r.n_agree_relaxed)} |"
            )
    lines.append("")
    lines.append("## Intervals used")
    lines.append("")
    lines.append("Full citations: `data/geological_ups/README.md`.")
    lines.append("")
    lines.append("| zone | segment | class | start Ma | end Ma | age basis (short) |")
    lines.append("|---|---|---|---:|---:|---|")
    show = iv.sort_values(["segment", "t_start_Ma"], ascending=[True, False])
    for r in show.itertuples(index=False):
        basis = str(r.age_basis).replace("|", "/")
        lines.append(
            f"| {r.zone_name} | {r.segment} | {r.ups_class} | {r.t_start_Ma:g} | {r.t_end_Ma:g} | {basis} |"
        )
    lines.append("")
    lines.append("## Segment × time coverage")
    lines.append("")
    cov = (
        blocks.groupby("segment")
        .agg(
            n_times=("time_Ma", "nunique"),
            tmin=("time_Ma", "min"),
            tmax=("time_Ma", "max"),
            classes=("geo_ups", lambda s: ",".join(sorted(s.unique()))),
            n_E=("geo_ups", lambda s: int((s == "E").sum())),
            n_C=("geo_ups", lambda s: int((s == "C").sum())),
            n_N=("geo_ups", lambda s: int((s == "N").sum())),
            med_age=("median_age_Ma", "median"),
            med_width=("median_width_km", "median"),
        )
        .reset_index()
        .sort_values("segment")
    )
    lines.append("| segment | n times | t range Ma | classes | nE | nC | nN | med age | med width km |")
    lines.append("|---|---:|---|---|---:|---:|---:|---:|---:|")
    for r in cov.itertuples(index=False):
        lines.append(
            f"| {r.segment} | {int(r.n_times)} | {r.tmax:g}–{r.tmin:g} | {r.classes} | "
            f"{int(r.n_E)} | {int(r.n_C)} | {int(r.n_N)} | {fmt_n(r.med_age)} | {fmt_n(r.med_width)} |"
        )
    lines.append("")
    lines.append("## Reading the result (no winner claim unless n is honest)")
    lines.append("")
    # Honest narrative filled from numbers, kept cautious.
    lines.append(stats["verdict"])
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append(f"- `{OUT_SUMMARY.relative_to(ROOT)}` — segment × time summary")
    lines.append(f"- `{FIG_AGE.relative_to(ROOT)}`")
    lines.append(f"- `{FIG_WIDTH.relative_to(ROOT)}`")
    lines.append(f"- `{FIG_TIMELINE.relative_to(ROOT)}`")
    lines.append("- `data/geological_ups/cenozoic_ups_intervals.csv`")
    lines.append("- `data/geological_ups/README.md`")
    lines.append("")
    lines.append("## Cite")
    lines.append("")
    lines.append("- Geological ages: `data/geological_ups/README.md` (primary papers per basin).")
    lines.append("- Müller, R. D., et al. (2019). *Tectonics*. https://doi.org/10.1029/2018TC005462")
    lines.append("- Sdrolias, M., & Müller, R. D. (2006). *G3*. https://doi.org/10.1029/2005GC001090")
    lines.append("- Schellart, W. P., et al. (2007). *Nature*. https://doi.org/10.1038/nature05737")
    lines.append("- Schellart, W. P. (2024). *Earth-Sci. Rev.* https://doi.org/10.1016/j.earscirev.2024.104755")
    lines.append("- SubMap: Cerpa, Lallemand & Heuret (2025). https://doi.org/10.31223/x56m8h")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def build_verdict(stats, submap, n_iv, n_cls, n_blocks) -> str:
    age_b = stats["age_blocks"]
    gate_b = stats["gate_blocks"]
    w_b = stats["width_blocks"]
    r_w, p_w, _ = stats["spearman_width"]
    bits = []
    bits.append(
        f"This is a **sparse hand catalog** ({n_iv} intervals; {n_cls} classified "
        f"sample-times; {n_blocks} segment×time blocks). Ages have published ranges."
    )
    e_older = (
        np.isfinite(age_b["median_a"])
        and np.isfinite(age_b["median_b"])
        and age_b["median_a"] > age_b["median_b"]
    )
    gate_dir = (
        np.isfinite(gate_b["odds"]) and gate_b["odds"] > 1.0 and gate_b["n_old"] >= 5
    )
    gate_flip = (
        np.isfinite(gate_b["odds"]) and gate_b["odds"] < 1.0 and gate_b["n_old"] >= 5
    )
    width_dir = (
        np.isfinite(w_b["median_a"])
        and np.isfinite(w_b["median_b"])
        and w_b["median_a"] < w_b["median_b"]
    )
    # Require honest n and a consistent direction before any "supports" language.
    strong_age = (
        e_older
        and gate_dir
        and age_b["n_a"] >= 8
        and age_b["n_b"] >= 8
        and np.isfinite(age_b["p"])
        and age_b["p"] < 0.05
        and np.isfinite(gate_b["p_fisher"])
        and gate_b["p_fisher"] < 0.05
    )
    strong_width = (
        width_dir
        and w_b["n_a"] >= 8
        and w_b["n_b"] >= 8
        and np.isfinite(r_w)
        and r_w < 0
        and np.isfinite(p_w)
        and p_w < 0.05
    )
    if strong_age:
        bits.append(
            "At segment×time, classified E sits on older slabs than C, and the 55 Ma "
            "gate has OR>1 with p<0.05. That is the Sdrolias direction, on independent "
            "UPS, but n is still the number of *segments×times* not 76k points, and C "
            "is mostly Andes."
        )
    elif gate_flip:
        bits.append(
            "At segment×time the 55 Ma gate runs the **opposite** way of Sdrolias & Müller "
            "(2006): P(E | age≥55) < P(E | age<55) (OR<1). Median age of E is younger than C. "
            "A large part of that is West Philippine / early IBM geological E sitting on "
            "young (≲40 Ma) Pacific lithosphere in this age grid, plus Andes C on older Nazca. "
            "Do not declare a Sdrolias-gate winner; if anything the independent catalog "
            "does not reproduce the present-day old-slab → back-arc association through "
            "the Cenozoic."
        )
    else:
        bits.append(
            "At segment×time the 55 Ma gate is **not** a clean, well-powered confirmation "
            "(see n and OR above). Do not declare a Sdrolias-gate winner from this table."
        )
    if strong_width:
        bits.append(
            "E blocks are narrower than C blocks in the Schellart direction (ρ(width, is_E)<0). "
            "That is consistent with Schellart et al. (2007 / 2024) but is also the geometry "
            "of boxing Scotia/Manus vs Andes. Not a knockout independent test of slab width."
        )
    else:
        bits.append(
            "Width / edge do **not** give a strong, well-powered Schellart confirmation "
            "at this n. No hypothesis winner."
        )
    if submap.get("ok") and submap["n_defined"]:
        frac = submap["n_agree_relaxed"] / submap["n_defined"]
        bits.append(
            f"0 Ma relaxed agreement with SubMap geodetic UPS is {frac:.2f} "
            f"({submap['n_agree_relaxed']}/{submap['n_defined']}). Cascadia N↔SS is counted "
            "as agree; Andaman E vs SS is a known mismatch (oblique spreading)."
        )
    bits.append(
        "**No hypothesis is declared the winner.** The point of this file is an "
        "independent UPS series plus honest n, not a knockout of Heuret–Lallemand vs Schellart."
    )
    return " ".join(bits)


def main(argv=None) -> int:
    args = parse_args(argv)
    _setup_style()
    print("Loading intervals ...", args.intervals)
    iv = load_intervals(args.intervals)
    print(f"  {len(iv)} intervals, segments={sorted(iv['segment'].unique())}")
    print("Loading kinematics ...", args.kin_csv)
    kin = pd.read_csv(args.kin_csv)
    print(f"  {len(kin)} rows, times={kin['time_Ma'].nunique()}")

    df = assign_geo_ups(kin, iv)
    n_cls = int(df["geo_classified"].sum())
    print(f"  classified sample-times: {n_cls} / {len(df)}")
    print("  class counts:", df.loc[df.geo_classified, "geo_ups"].value_counts().to_dict())

    blocks = aggregate_segment_time(df)
    print(f"  segment×time blocks: {len(blocks)}")
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    blocks.to_csv(args.out_summary, index=False)
    print("  wrote", args.out_summary)

    # ---- stats ----
    samp = df.loc[df["geo_classified"]].copy()
    samp_age = samp.loc[np.isfinite(samp["seafloor_age_Ma"])]
    age_s = mwu(
        samp_age.loc[samp_age["geo_ups"] == "E", "seafloor_age_Ma"],
        samp_age.loc[samp_age["geo_ups"] == "C", "seafloor_age_Ma"],
    )
    blk_age = blocks.loc[np.isfinite(blocks["median_age_Ma"])]
    age_b = mwu(
        blk_age.loc[blk_age["geo_ups"] == "E", "median_age_Ma"],
        blk_age.loc[blk_age["geo_ups"] == "C", "median_age_Ma"],
    )
    age_b_enot = mwu(
        blk_age.loc[blk_age["geo_ups"] == "E", "median_age_Ma"],
        blk_age.loc[blk_age["geo_ups"] != "E", "median_age_Ma"],
    )
    gate_s = gate_table(
        samp_age["geo_ups"] == "E",
        samp_age["seafloor_age_Ma"],
        np.ones(len(samp_age), dtype=bool),
    )
    samp_ec = samp_age.loc[samp_age["geo_ups"].isin(["E", "C"])]
    gate_s_ec = gate_table(
        samp_ec["geo_ups"] == "E",
        samp_ec["seafloor_age_Ma"],
        np.ones(len(samp_ec), dtype=bool),
    )
    gate_b = gate_table(
        blk_age["geo_ups"] == "E",
        blk_age["median_age_Ma"],
        np.ones(len(blk_age), dtype=bool),
    )
    blk_ec = blk_age.loc[blk_age["geo_ups"].isin(["E", "C"])]
    gate_b_ec = gate_table(
        blk_ec["geo_ups"] == "E",
        blk_ec["median_age_Ma"],
        np.ones(len(blk_ec), dtype=bool),
    )
    w_b = mwu(
        blocks.loc[blocks["geo_ups"] == "E", "median_width_km"],
        blocks.loc[blocks["geo_ups"] == "C", "median_width_km"],
    )
    e_b = mwu(
        blocks.loc[blocks["geo_ups"] == "E", "median_edge_km"],
        blocks.loc[blocks["geo_ups"] == "C", "median_edge_km"],
    )
    sp_w = spearman_flag(blocks["median_width_km"], blocks["is_E"])
    sp_e = spearman_flag(blocks["median_edge_km"], blocks["is_E"])
    near = blocks["median_edge_km"] < EDGE_NEAR_KM
    n_near = int(near.sum())
    n_far = int((~near).sum())
    n_e_near = int(((blocks["geo_ups"] == "E") & near).sum())
    n_e_far = int(((blocks["geo_ups"] == "E") & ~near).sum())
    near_edge = {
        "n_near": n_near,
        "n_far": n_far,
        "n_e_near": n_e_near,
        "n_e_far": n_e_far,
        "p_near": n_e_near / n_near if n_near else np.nan,
        "p_far": n_e_far / n_far if n_far else np.nan,
    }
    # Sensitivity: drop Tethyan Hellenic/Calabria (Mesozoic ages, not Pacific BABs).
    tethyan = {"Hellenic", "Calabria"}
    blk_nt = blk_age.loc[~blk_age["segment"].isin(tethyan)]
    age_b_nt = mwu(
        blk_nt.loc[blk_nt["geo_ups"] == "E", "median_age_Ma"],
        blk_nt.loc[blk_nt["geo_ups"] == "C", "median_age_Ma"],
    )
    gate_b_nt = gate_table(
        blk_nt["geo_ups"] == "E",
        blk_nt["median_age_Ma"],
        np.ones(len(blk_nt), dtype=bool),
    )
    stats = {
        "age_samples": age_s,
        "age_blocks": age_b,
        "age_blocks_e_not": age_b_enot,
        "age_blocks_notethys": age_b_nt,
        "gate_samples": gate_s,
        "gate_samples_ec": gate_s_ec,
        "gate_blocks": gate_b,
        "gate_blocks_ec": gate_b_ec,
        "gate_blocks_notethys": gate_b_nt,
        "width_blocks": w_b,
        "edge_blocks": e_b,
        "spearman_width": sp_w,
        "spearman_edge": sp_e,
        "near_edge": near_edge,
    }

    print("0 Ma SubMap sanity ...")
    df0 = df.loc[df["time_Ma"] == 0.0]
    submap = submap_sanity(df0, args.join_csv, iv)
    stats["verdict"] = build_verdict(stats, submap, len(iv), n_cls, len(blocks))

    print("Figures ...")
    fig_age(df, blocks, FIG_AGE)
    fig_width_edge(blocks, FIG_WIDTH)
    fig_timeline(iv, FIG_TIMELINE)
    print("  wrote", FIG_AGE)
    print("  wrote", FIG_WIDTH)
    print("  wrote", FIG_TIMELINE)

    write_report(args.out_md, iv, df, blocks, stats, submap)
    print("  wrote", args.out_md)

    # stdout summary for the user
    print("\n========== SHORT SUMMARY ==========")
    print(f"n intervals: {len(iv)}")
    print(f"n classified sample-times: {n_cls} / {len(df)}")
    print(
        f"median age E vs C (samples): {fmt_n(age_s['median_a'])} vs {fmt_n(age_s['median_b'])} Ma "
        f"(n={age_s['n_a']}/{age_s['n_b']})"
    )
    print(
        f"median age E vs C (segment×time): {fmt_n(age_b['median_a'])} vs {fmt_n(age_b['median_b'])} Ma "
        f"(n={age_b['n_a']}/{age_b['n_b']}; p={fmt_p(age_b['p'])})"
    )
    print(
        f"gate samples  P(E|≥55)={gate_s['p_old']:.3f} vs P(E|<55)={gate_s['p_young']:.3f}  "
        f"OR={fmt_n(gate_s['odds'])}"
    )
    print(
        f"gate blocks   P(E|≥55)={gate_b['p_old']:.3f} vs P(E|<55)={gate_b['p_young']:.3f}  "
        f"OR={fmt_n(gate_b['odds'])}  p={fmt_p(gate_b['p_fisher'])}"
    )
    print(
        f"width E vs C (blocks): {fmt_n(w_b['median_a'])} vs {fmt_n(w_b['median_b'])} km  "
        f"ρ(width,is_E)={sp_w[0]:.3f} p={fmt_p(sp_w[1])}"
    )
    if submap.get("ok") and submap["n_defined"]:
        print(
            f"0 Ma SubMap relaxed agreement: {submap['n_agree_relaxed']}/{submap['n_defined']} "
            f"= {submap['n_agree_relaxed']/submap['n_defined']:.2f}"
        )
    else:
        print("0 Ma SubMap: no overlap / missing join")
    print(stats["verdict"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
