#!/usr/bin/env python3
"""Cenozoic upper-plate strain PROXY from Müller et al. (2019) topologies.

Joins a back-arc spreading on/off flag (layer A) to the existing trench
catalog and, where deforming networks exist, samples dilatation ~200 km
into the overriding plate (layer B). Then tests the Sdrolias–Müller
(~55 Ma) age gate plus trench-zone width and distance-to-edge.

This proxy is built from the SAME reconstruction as the age grids. It is a
Sdrolias & Müller (2006)-style kinematic association, not independent
geology. Do **not** treat bab_on as an observed UPS time series.

Does **not** re-extract trench kinematics. Does **not** use v_T vs v_OP as UPS.

Run:
    /workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
        /workspace/subduction-test/scripts/cenozoic_ups_proxy.py

0 Ma only (sanity / prototype):
    .../cenozoic_ups_proxy.py --tmax 0
"""

from __future__ import annotations

import argparse
import sys
import time as time_mod
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pygplates
from scipy.stats import chi2, mannwhitneyu, spearmanr
from sklearn.linear_model import LogisticRegression

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import extract_trench_kinematics_0Ma as k0  # noqa: E402
from gplately import PlateModelManager, PlateReconstruction

# ---------------------------------------------------------------------------
# Paths / defaults
# ---------------------------------------------------------------------------
ROOT = Path("/workspace/subduction-test")
DATA_DIR = k0.DATA_DIR
RESULTS_DIR = k0.RESULTS_DIR
KIN_CSV = RESULTS_DIR / "trench_kinematics_0-60Ma_Muller2019.csv"
JOIN_CSV = RESULTS_DIR / "joined_submap_muller2019_0Ma.csv"
OUT_CSV = RESULTS_DIR / "trench_kinematics_0-60Ma_with_ups_proxy.csv"
REPORT_MD = RESULTS_DIR / "cenozoic_age_gate.md"
FIG_BAB_AGE = RESULTS_DIR / "fig_bab_vs_age.png"
FIG_GATE = RESULTS_DIR / "fig_age_gate_through_time.png"
FIG_MAP = RESULTS_DIR / "fig_bab_map_0Ma.png"

EARTH_RADIUS_KM = k0.EARTH_RADIUS_KM
TESSELLATION_DEG = k0.TESSELLATION_DEG
MAX_DIST_KM = 1500.0
MIN_ORTH_CM_YR = 1.0
# Practical along-trench cap so a major-ocean MOR that sits *along-strike*
# of the overriding plate (Gulf of California vs Middle America) is not
# counted as a back-arc. True BABs sit behind the local trench segment
# (nearest ridge typically |d_parallel| ≲ 150 km).
MAX_ALONG_TRENCH_KM = 400.0
DILAT_OFFSET_KM = 200.0
AGE_GATE_MA = 55.0
AGE_SWEEP_MA = (30, 40, 45, 50, 55, 60, 70, 80, 100)
SECONDS_PER_MYR = 365.25 * 24.0 * 3600.0 * 1.0e6

# 0 Ma spot-check boxes (lon_min, lon_max, lat_min, lat_max).
SANITY_REGIONS = {
    "Mariana": (140.0, 150.0, 10.0, 25.0),
    "Tonga": (-180.0, -170.0, -30.0, -15.0),
    "Scotia": (-35.0, -20.0, -62.0, -54.0),
    "Hellenic": (18.0, 30.0, 33.0, 40.0),
    "Andes": (-80.0, -68.0, -40.0, -5.0),
    "Cascadia": (-130.0, -120.0, 40.0, 51.0),
    "Andaman": (90.0, 98.0, 5.0, 15.0),
    "MAT_Mexico": (-106.0, -99.0, 15.0, 21.0),
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Cenozoic UPS proxy (back-arc spreading) from Müller2019."
    )
    p.add_argument("--kin-csv", type=Path, default=KIN_CSV)
    p.add_argument("--out-csv", type=Path, default=OUT_CSV)
    p.add_argument("--tmin", type=int, default=None, help="youngest time (Ma)")
    p.add_argument("--tmax", type=int, default=None, help="oldest time (Ma)")
    p.add_argument(
        "--skip-dilatation",
        action="store_true",
        help="skip layer B (deforming-network dilatation)",
    )
    p.add_argument(
        "--skip-sanity-abort",
        action="store_true",
        help="do not abort the 0–60 run if the 0 Ma spot check fails",
    )
    p.add_argument(
        "--from-csv",
        type=Path,
        default=None,
        help="skip ridge tessellation; recompute stats/figures from this proxy CSV",
    )
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


def load_reconstruction():
    print("Loading Müller2019 via PlateModelManager ...")
    model = PlateModelManager().get_model("Muller2019", data_dir=str(DATA_DIR))
    recon = PlateReconstruction(
        model.get_rotation_model(),
        model.get_topologies(),
        model.get_static_polygons(),
    )
    return model, recon


def tessellate_ridges(recon, time):
    """Tessellate mid-ocean ridges. Raises if the call fails or returns nothing."""
    try:
        gdf = recon.tessellate_mid_ocean_ridges(
            float(time),
            tessellation_threshold_radians=np.radians(TESSELLATION_DEG),
            ignore_warnings=True,
            return_geodataframe=True,
            output_obliquity_and_normal_and_left_right_plates=True,
        )
    except Exception as err:
        raise RuntimeError(
            f"Ridge tessellation failed at {time} Ma: {err}"
        ) from err
    if gdf is None or len(gdf) == 0:
        raise RuntimeError(
            f"Ridge tessellation returned no segments at {time} Ma. Stopping."
        )
    lon = gdf.geometry.x.to_numpy(dtype=float)
    lat = gdf.geometry.y.to_numpy(dtype=float)
    vel = gdf["velocity (cm/yr)"].to_numpy(dtype=float)
    obliq = gdf["obliquity (degrees)"].to_numpy(dtype=float)
    left = gdf["left plate ID"].to_numpy(dtype=int)
    right = gdf["right plate ID"].to_numpy(dtype=int)
    orth = vel * np.cos(np.radians(obliq))
    return lon, lat, left, right, orth


def ridge_containing_plate_ids(snapshot, rlon, rlat):
    """Reconstruction plate ID of the topology containing each ridge point."""
    pts = list(zip(rlat.tolist(), rlon.tolist()))
    locs = snapshot.get_point_locations(pts)
    pids = np.full(len(pts), np.nan)
    for i, loc in enumerate(locs):
        feat = None
        b = loc.located_in_resolved_boundary()
        n = loc.located_in_resolved_network()
        if b:
            feat = b.get_feature()
        elif n:
            feat = n.get_feature()
        if feat is None:
            continue
        try:
            pids[i] = int(feat.get_reconstruction_plate_id())
        except Exception:
            continue
    return pids


def _equirect_en(tlon, tlat, rlon, rlat):
    """East/north km of ridges relative to trenches. Shapes (N,) and (M,) → (N, M)."""
    dlat = np.radians(rlat)[None, :] - np.radians(tlat)[:, None]
    dlon = np.radians(rlon)[None, :] - np.radians(tlon)[:, None]
    dlon = (dlon + np.pi) % (2.0 * np.pi) - np.pi
    mean_lat = 0.5 * (np.radians(rlat)[None, :] + np.radians(tlat)[:, None])
    east = EARTH_RADIUS_KM * dlon * np.cos(mean_lat)
    north = EARTH_RADIUS_KM * dlat
    return east, north


def match_backarc(
    tlon,
    tlat,
    taz_deg,
    trench_id,
    sub_id,
    rlon,
    rlat,
    rleft,
    rright,
    rorth,
    rtopo_pid,
):
    """Vectorised trench↔ridge match for one time.

    A ridge counts as a back-arc of a trench sample if:
      * orthogonal spreading ≥ 1 cm/yr
      * great-circle distance ≤ 1500 km
      * landward of the trench (d_landward > 0; trench normal toward OP)
      * |along-trench offset| ≤ 400 km (kills along-strike major MORs)
      * on the overriding plate: left or right plate ID == trench_plate_id,
        or the ridge point lies in that plate
      * not on the subducting plate (left/right/containing ID)
    """
    n = len(tlon)
    if len(rlon) == 0:
        return (
            np.zeros(n, dtype=int),
            np.full(n, np.nan),
            np.full(n, np.nan),
        )

    active = np.isfinite(rorth) & (rorth >= MIN_ORTH_CM_YR)
    if not np.any(active):
        return (
            np.zeros(n, dtype=int),
            np.full(n, np.nan),
            np.full(n, np.nan),
        )
    rlon = rlon[active]
    rlat = rlat[active]
    rleft = rleft[active]
    rright = rright[active]
    rorth = rorth[active]
    if rtopo_pid is not None:
        rtopo_pid = rtopo_pid[active]

    east, north = _equirect_en(tlon, tlat, rlon, rlat)
    dist = np.hypot(east, north)
    az = np.radians(taz_deg)
    caz = np.cos(az)[:, None]
    saz = np.sin(az)[:, None]
    d_land = north * caz + east * saz
    d_par = north * (-saz) + east * caz

    plate = (rleft[None, :] == trench_id[:, None]) | (
        rright[None, :] == trench_id[:, None]
    )
    if rtopo_pid is not None:
        finite = np.isfinite(rtopo_pid)
        if np.any(finite):
            plate[:, finite] |= rtopo_pid[finite][None, :] == trench_id[:, None]

    on_sub = (rleft[None, :] == sub_id[:, None]) | (
        rright[None, :] == sub_id[:, None]
    )
    if rtopo_pid is not None:
        finite = np.isfinite(rtopo_pid)
        if np.any(finite):
            on_sub[:, finite] |= rtopo_pid[finite][None, :] == sub_id[:, None]

    valid = (
        (dist <= MAX_DIST_KM)
        & (d_land > 0.0)
        & (np.abs(d_par) <= MAX_ALONG_TRENCH_KM)
        & plate
        & ~on_sub
    )
    bab_on = valid.any(axis=1).astype(int)
    dist_valid = np.where(valid, dist, np.inf)
    min_j = dist_valid.argmin(axis=1)
    min_dist = dist_valid[np.arange(n), min_j]
    min_orth = rorth[min_j].astype(float)
    none = ~np.isfinite(min_dist)
    min_dist[none] = np.nan
    min_orth[none] = np.nan
    return bab_on, min_dist, min_orth


def sample_dilatation(snapshot, tlon, tlat, taz_deg):
    """Dilatation rate (1e-15 / s) ~200 km into the overriding plate."""
    olon, olat = k0._destination_lonlat(tlon, tlat, taz_deg, DILAT_OFFSET_KM)
    pts = list(zip(olat.tolist(), olon.tolist()))
    strain_rates, locations = snapshot.get_point_strain_rates(
        pts, return_point_locations=True
    )
    n = len(pts)
    dilat = np.full(n, np.nan)
    in_net = np.full(n, np.nan)
    n_net = 0
    n_bound = 0
    n_out = 0
    for i, (sr, loc) in enumerate(zip(strain_rates, locations)):
        if sr is None:
            n_out += 1
            continue
        dilat[i] = sr.get_dilatation_rate() * 1.0e15
        if loc.located_in_resolved_network():
            in_net[i] = 1.0
            n_net += 1
        else:
            in_net[i] = 0.0
            n_bound += 1
    stats = {"n_network": n_net, "n_boundary": n_bound, "n_outside": n_out}
    return dilat, in_net, stats


def proxy_one_time(recon, time, trench, do_dilatation=True):
    """Compute layer A (and optionally B) for one reconstruction time."""
    rlon, rlat, rleft, rright, rorth = tessellate_ridges(recon, time)
    n_ridge_raw = len(rlon)
    n_ridge_active = int(np.sum(rorth >= MIN_ORTH_CM_YR))

    snapshot = recon.topological_snapshot(float(time))
    rtopo = ridge_containing_plate_ids(snapshot, rlon, rlat)

    tlon = trench["lon"].to_numpy(dtype=float)
    tlat = trench["lat"].to_numpy(dtype=float)
    taz = trench["trench_normal_azimuth_deg"].to_numpy(dtype=float)
    tid = trench["trench_plate_id"].to_numpy(dtype=int)
    sid = trench["subducting_plate_id"].to_numpy(dtype=int)

    bab_on, dist_ba, orth_ba = match_backarc(
        tlon, tlat, taz, tid, sid, rlon, rlat, rleft, rright, rorth, rtopo
    )

    dilat = np.full(len(trench), np.nan)
    in_net = np.full(len(trench), np.nan)
    dil_stats = {"n_network": 0, "n_boundary": 0, "n_outside": len(trench)}
    if do_dilatation:
        dilat, in_net, dil_stats = sample_dilatation(snapshot, tlon, tlat, taz)

    out = trench.copy()
    out["bab_on"] = bab_on
    out["dist_to_backarc_ridge_km"] = dist_ba
    out["backarc_spreading_cm_yr"] = orth_ba
    out["dilat_rate_1e-15_s"] = dilat
    out["dilat_in_network"] = in_net
    meta = {
        "time": float(time),
        "n_trench": int(len(trench)),
        "n_ridge_raw": int(n_ridge_raw),
        "n_ridge_active": int(n_ridge_active),
        "n_bab_on": int(bab_on.sum()),
        **dil_stats,
    }
    return out, meta


def region_mask(df, box):
    lo, hi, la, lb = box
    return (df["lon"] >= lo) & (df["lon"] <= hi) & (df["lat"] >= la) & (df["lat"] <= lb)


def sanity_table_0ma(df0):
    rows = []
    for name, box in SANITY_REGIONS.items():
        m = region_mask(df0, box)
        s = df0.loc[m]
        n = int(m.sum())
        n_on = int(s["bab_on"].sum()) if n else 0
        frac = (n_on / n) if n else np.nan
        dmed = (
            float(s.loc[s["bab_on"] == 1, "dist_to_backarc_ridge_km"].median())
            if n_on
            else np.nan
        )
        dil_med = (
            float(s["dilat_rate_1e-15_s"].median()) if n else np.nan
        )
        dil_pos = (
            float((s["dilat_rate_1e-15_s"] > 0).mean()) if n else np.nan
        )
        rows.append(
            {
                "region": name,
                "n": n,
                "n_bab_on": n_on,
                "frac_bab_on": frac,
                "med_dist_ba_km": dmed,
                "med_dilat_1e-15_s": dil_med,
                "frac_dilat_pos": dil_pos,
            }
        )
    return pd.DataFrame(rows)


def check_sanity_0ma(table, abort=True):
    """Mariana/Tonga/Scotia on; Andes/Cascadia off. Hellenic is deforming, not MOR."""
    want_on = {"Mariana": 0.5, "Tonga": 0.5, "Scotia": 0.5}
    want_off = {"Andes": 0.15, "Cascadia": 0.15, "MAT_Mexico": 0.15}
    problems = []
    by = table.set_index("region")
    for name, thresh in want_on.items():
        frac = float(by.loc[name, "frac_bab_on"])
        if not np.isfinite(frac) or frac < thresh:
            problems.append(f"{name} frac_bab_on={frac:.2f} (want ≥ {thresh})")
    for name, thresh in want_off.items():
        frac = float(by.loc[name, "frac_bab_on"])
        if np.isfinite(frac) and frac > thresh:
            problems.append(f"{name} frac_bab_on={frac:.2f} (want ≤ {thresh})")
    hell = float(by.loc["Hellenic", "frac_bab_on"])
    print("\n=== 0 Ma spot check (layer A = MOR back-arc) ===")
    print(table.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    if problems:
        msg = "0 Ma sanity FAILED: " + "; ".join(problems)
        print(msg)
        if abort:
            raise RuntimeError(msg)
    else:
        print(
            "0 Ma sanity passed: Mariana/Tonga/Scotia ON; "
            "Andes/Cascadia/MAT OFF."
        )
        if hell < 0.15:
            print(
                "Note: Hellenic bab_on=0 — Müller2019 encodes Aegean extension "
                "as a deforming network, not a mid-ocean ridge. Layer B "
                "(dilatation) is the relevant proxy there."
            )
    return problems


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
def _finite_pair(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


def logit_univariate(x, y):
    """Unpenalized logistic y ~ x. Returns dict with coef, OR per unit, p (LRT)."""
    x, y = _finite_pair(x, y)
    if x.size < 20 or np.unique(y).size < 2 or np.nanstd(x) == 0:
        return None
    X = x.reshape(-1, 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = LogisticRegression(
            penalty=None, solver="lbfgs", max_iter=2000
        ).fit(X, y)
        null = LogisticRegression(
            penalty=None, solver="lbfgs", max_iter=2000
        ).fit(np.ones((len(y), 1)), y)
    coef = float(full.coef_[0, 0])
    intercept = float(full.intercept_[0])
    p_full = full.predict_proba(X)[:, 1]
    p_null = null.predict_proba(np.ones((len(y), 1)))[:, 1]
    p_full = np.clip(p_full, 1e-12, 1 - 1e-12)
    p_null = np.clip(p_null, 1e-12, 1 - 1e-12)
    ll_full = float(np.sum(y * np.log(p_full) + (1 - y) * np.log(1 - p_full)))
    ll_null = float(np.sum(y * np.log(p_null) + (1 - y) * np.log(1 - p_null)))
    lrt = 2.0 * (ll_full - ll_null)
    p_lrt = float(chi2.sf(max(lrt, 0.0), 1))
    pseudo_r2 = float(1.0 - (ll_full / ll_null)) if ll_null != 0 else np.nan
    return {
        "n": int(len(y)),
        "coef": coef,
        "intercept": intercept,
        "odds_ratio_per_unit": float(np.exp(coef)),
        "lrt": lrt,
        "p_lrt": p_lrt,
        "mcfadden_r2": pseudo_r2,
        "n_on": int(y.sum()),
    }


def logit_age_with_time(df, age_col, y_col, time_col="time_Ma"):
    """Logistic y ~ age + time dummies. Returns age coefficient and LRT vs time-only."""
    d = df[[age_col, y_col, time_col]].dropna()
    if len(d) < 30 or d[y_col].nunique() < 2:
        return None
    y = d[y_col].to_numpy(dtype=float)
    age = d[age_col].to_numpy(dtype=float)
    times = d[time_col].to_numpy(dtype=float)
    uniq = np.sort(np.unique(times))
    if uniq.size < 2:
        # single time slice: reduce to univariate logistic
        uni = logit_univariate(age, y)
        if uni is None:
            return None
        return {
            "n": uni["n"],
            "n_on": uni["n_on"],
            "coef_age": uni["coef"],
            "odds_ratio_per_Ma": uni["odds_ratio_per_unit"],
            "lrt_vs_time_only": uni["lrt"],
            "p_lrt": uni["p_lrt"],
        }
    # drop first time as reference
    dummies = np.column_stack([(times == t).astype(float) for t in uniq[1:]])
    X_full = np.column_stack([age, dummies])
    X_red = dummies
    if X_red.ndim == 1:
        X_red = X_red.reshape(-1, 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = LogisticRegression(
            penalty=None, solver="lbfgs", max_iter=4000
        ).fit(X_full, y)
        red = LogisticRegression(
            penalty=None, solver="lbfgs", max_iter=4000
        ).fit(X_red, y)
    coef_age = float(full.coef_[0, 0])
    p_full = np.clip(full.predict_proba(X_full)[:, 1], 1e-12, 1 - 1e-12)
    p_red = np.clip(red.predict_proba(X_red)[:, 1], 1e-12, 1 - 1e-12)
    ll_full = float(np.sum(y * np.log(p_full) + (1 - y) * np.log(1 - p_full)))
    ll_red = float(np.sum(y * np.log(p_red) + (1 - y) * np.log(1 - p_red)))
    lrt = 2.0 * (ll_full - ll_red)
    p_lrt = float(chi2.sf(max(lrt, 0.0), 1))
    return {
        "n": int(len(y)),
        "n_on": int(y.sum()),
        "coef_age": coef_age,
        "odds_ratio_per_Ma": float(np.exp(coef_age)),
        "lrt_vs_time_only": lrt,
        "p_lrt": p_lrt,
    }


def mantel_haenszel(tables):
    """tables: list of (a,b,c,d) = [[on|old, off|old],[on|young, off|young]] counts.

    Returns MH odds ratio and chi-square p (without continuity correction).
    """
    num = 0.0
    den = 0.0
    chi_num = 0.0
    chi_den = 0.0
    k = 0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n <= 1:
            continue
        k += 1
        num += a * d / n
        den += b * c / n
        n1 = a + b
        m1 = a + c
        ea = n1 * m1 / n
        va = n1 * (c + d) * m1 * (b + d) / (n * n * (n - 1)) if n > 1 else 0.0
        chi_num += a - ea
        chi_den += va
    if den == 0 or k == 0:
        return None
    or_mh = num / den
    chi = (chi_num ** 2) / chi_den if chi_den > 0 else np.nan
    p = float(chi2.sf(chi, 1)) if np.isfinite(chi) else np.nan
    return {"or_mh": float(or_mh), "chi2": float(chi), "p": p, "n_strata": k}


def gate_counts(age, bab, threshold):
    m = np.isfinite(age) & np.isfinite(bab)
    age = age[m]
    bab = bab[m].astype(int)
    old = age >= threshold
    young = ~old
    on_old = int(np.sum(old & (bab == 1)))
    off_old = int(np.sum(old & (bab == 0)))
    on_young = int(np.sum(young & (bab == 1)))
    off_young = int(np.sum(young & (bab == 0)))
    n_old = on_old + off_old
    n_young = on_young + off_young
    p_old = on_old / n_old if n_old else np.nan
    p_young = on_young / n_young if n_young else np.nan
    return {
        "threshold": threshold,
        "n_old": n_old,
        "n_young": n_young,
        "on_old": on_old,
        "off_old": off_old,
        "on_young": on_young,
        "off_young": off_young,
        "p_old": p_old,
        "p_young": p_young,
        "delta": (p_old - p_young) if np.isfinite(p_old) and np.isfinite(p_young) else np.nan,
    }


def zone_frame(df):
    """One row per (time_Ma, zone_id). bab_zone = 1 if any sample is bab_on."""
    g = df.groupby(["time_Ma", "zone_id"], sort=False)
    z = g.agg(
        n_points=("bab_on", "size"),
        n_bab_on=("bab_on", "sum"),
        bab_zone=("bab_on", "max"),
        seafloor_age_Ma=("seafloor_age_Ma", "median"),
        trench_zone_width_km=("trench_zone_width_km", "median"),
        dist_nearest_trench_edge_km=("dist_nearest_trench_edge_km", "median"),
        dilat_rate=("dilat_rate_1e-15_s", "median"),
        lon=("lon", "median"),
        lat=("lat", "median"),
    ).reset_index()
    z = z.rename(columns={"dilat_rate": "dilat_rate_1e-15_s"})
    z["bab_zone"] = z["bab_zone"].astype(int)
    z["bab_frac"] = z["n_bab_on"] / z["n_points"]
    return z


def per_time_stats(df):
    rows = []
    for t, s in df.groupby("time_Ma"):
        age = s["seafloor_age_Ma"].to_numpy(dtype=float)
        bab = s["bab_on"].to_numpy(dtype=float)
        m = np.isfinite(age)
        age_f, bab_f = age[m], bab[m]
        med_on = float(np.median(age_f[bab_f == 1])) if np.any(bab_f == 1) else np.nan
        med_off = float(np.median(age_f[bab_f == 0])) if np.any(bab_f == 0) else np.nan
        gate = gate_counts(age_f, bab_f, AGE_GATE_MA)
        sp = None
        if age_f.size >= 20 and np.unique(bab_f).size == 2:
            rho, p = spearmanr(age_f, bab_f)
            sp = (float(rho), float(p))
        rows.append(
            {
                "time_Ma": float(t),
                "n": int(len(s)),
                "n_age": int(m.sum()),
                "n_bab_on": int(s["bab_on"].sum()),
                "frac_bab_on": float(s["bab_on"].mean()),
                "med_age_on": med_on,
                "med_age_off": med_off,
                "p_bab_old": gate["p_old"],
                "p_bab_young": gate["p_young"],
                "n_old": gate["n_old"],
                "n_young": gate["n_young"],
                "spearman_age_rho": sp[0] if sp else np.nan,
                "spearman_age_p": sp[1] if sp else np.nan,
            }
        )
    return pd.DataFrame(rows)


def univariate_block(df, ycol, label):
    """Spearman + logistic for age, width, dist_to_edge vs binary y."""
    rows = []
    for pred, pretty in [
        ("seafloor_age_Ma", "seafloor age (Ma)"),
        ("trench_zone_width_km", "trench-zone width (km)"),
        ("dist_nearest_trench_edge_km", "distance to trench edge (km)"),
    ]:
        x, y = _finite_pair(df[pred], df[ycol])
        if x.size < 20 or np.unique(y).size < 2:
            rows.append({"level": label, "predictor": pretty, "n": int(x.size)})
            continue
        rho, p_s = spearmanr(x, y)
        logit = logit_univariate(x, y)
        rows.append(
            {
                "level": label,
                "predictor": pretty,
                "n": int(x.size),
                "n_on": int(y.sum()),
                "spearman_rho": float(rho),
                "spearman_p": float(p_s),
                "logit_coef": logit["coef"] if logit else np.nan,
                "logit_OR": logit["odds_ratio_per_unit"] if logit else np.nan,
                "logit_p_lrt": logit["p_lrt"] if logit else np.nan,
                "logit_r2": logit["mcfadden_r2"] if logit else np.nan,
            }
        )
    return pd.DataFrame(rows)


def submap_sanity(df0):
    """Join 0 Ma bab_on onto SubMap rows; crosstab vs UPS_3 (E vs C)."""
    if not JOIN_CSV.is_file():
        return None, "joined SubMap CSV not found"
    sub = pd.read_csv(JOIN_CSV)

    def hav(lon1, lat1, lon2, lat2):
        lon1 = np.radians(np.asarray(lon1, dtype=float))[:, None]
        lat1 = np.radians(np.asarray(lat1, dtype=float))[:, None]
        lon2 = np.radians(np.asarray(lon2, dtype=float))[None, :]
        lat2 = np.radians(np.asarray(lat2, dtype=float))[None, :]
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
        a = np.clip(a, 0.0, 1.0)
        return EARTH_RADIUS_KM * 2.0 * np.arcsin(np.sqrt(a))

    d = hav(sub["muller_lon"], sub["muller_lat"], df0["lon"], df0["lat"])
    idx = d.argmin(axis=1)
    dmin = d[np.arange(len(sub)), idx]
    sub = sub.copy()
    sub["bab_on"] = df0["bab_on"].to_numpy()[idx]
    sub["proxy_join_km"] = dmin
    # keep already-matched SubMap rows (join_distance_km exists)
    ct_full = pd.crosstab(sub["UPS_3"], sub["bab_on"], dropna=False)
    ce = sub[sub["UPS_3"].isin(["C", "E"])].copy()
    ct_ce = pd.crosstab(ce["UPS_3"], ce["bab_on"])
    # rates
    rates = (
        ce.groupby("UPS_3")["bab_on"]
        .agg(n="size", n_bab_on="sum", frac_bab_on="mean")
        .reindex(["C", "E"])
    )
    geo = (
        sub.groupby("Geodetic_UPS")["bab_on"]
        .agg(n="size", n_bab_on="sum", frac_bab_on="mean")
    )
    return {
        "n_sub": int(len(sub)),
        "med_proxy_join_km": float(np.median(dmin)),
        "ct_full": ct_full,
        "ct_ce": ct_ce,
        "rates_ce": rates,
        "geo": geo,
        "n_ce": int(len(ce)),
    }, None


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def fig_bab_vs_age(df, per_time, path):
    _setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))

    ax = axes[0]
    on = df.loc[df["bab_on"] == 1, "seafloor_age_Ma"].to_numpy(dtype=float)
    off = df.loc[df["bab_on"] == 0, "seafloor_age_Ma"].to_numpy(dtype=float)
    on = on[np.isfinite(on)]
    off = off[np.isfinite(off)]
    bins = np.linspace(0, min(200, np.nanpercentile(np.concatenate([on, off]), 98)), 40)
    ax.hist(
        off,
        bins=bins,
        density=True,
        histtype="stepfilled",
        color="0.75",
        edgecolor="0.4",
        alpha=0.9,
        label=f"bab_off  n={len(off)}  med={np.median(off):.0f} Ma",
    )
    ax.hist(
        on,
        bins=bins,
        density=True,
        histtype="stepfilled",
        color="#0A9396",
        edgecolor="#005F73",
        alpha=0.55,
        label=f"bab_on   n={len(on)}  med={np.median(on):.0f} Ma",
    )
    ax.axvline(AGE_GATE_MA, color="#9B2226", ls="--", lw=1.2, label="55 Ma gate")
    ax.set_xlabel("Seafloor age at trench (Ma)")
    ax.set_ylabel("Density")
    ax.set_xlim(0, bins[-1])
    ax.set_title("Pooled 0–60 Ma (point samples)")
    ax.legend(fontsize=8, frameon=True, loc="upper right")
    ax.grid(True, linestyle=":", linewidth=0.5, color="0.8")

    ax = axes[1]
    t = per_time["time_Ma"].to_numpy(dtype=float)
    ax.plot(
        t,
        per_time["med_age_on"],
        color="#0A9396",
        lw=2.0,
        marker="o",
        ms=3.0,
        label="median age | bab_on",
    )
    ax.plot(
        t,
        per_time["med_age_off"],
        color="0.45",
        lw=2.0,
        marker="o",
        ms=3.0,
        label="median age | bab_off",
    )
    ax.axhline(AGE_GATE_MA, color="#9B2226", ls="--", lw=1.0)
    ax.set_xlabel("Reconstruction time (Ma)")
    ax.set_ylabel("Median seafloor age (Ma)")
    t0, t1 = float(np.nanmin(t)), float(np.nanmax(t))
    ax.set_xlim(t1 + (0.5 if t1 == t0 else 0.0), t0 - (0.5 if t1 == t0 else 0.0))
    ax.set_title("Median age vs time")
    ax.legend(fontsize=8, frameon=True, loc="upper left")
    ax.grid(True, linestyle=":", linewidth=0.5, color="0.8")

    fig.suptitle(
        "Müller et al. (2019) — seafloor age at trench vs back-arc spreading proxy",
        y=1.02,
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def fig_age_gate(per_time, path):
    _setup_style()
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    t = per_time["time_Ma"].to_numpy(dtype=float)
    ax.plot(
        t,
        per_time["p_bab_old"],
        color="#005F73",
        lw=2.1,
        marker="o",
        ms=3.2,
        label=r"$P(\mathrm{bab\_on}\mid \mathrm{age}\geq 55\,\mathrm{Ma})$",
    )
    ax.plot(
        t,
        per_time["p_bab_young"],
        color="#BB3E03",
        lw=2.1,
        marker="o",
        ms=3.2,
        label=r"$P(\mathrm{bab\_on}\mid \mathrm{age}< 55\,\mathrm{Ma})$",
    )
    ax.plot(
        t,
        per_time["frac_bab_on"],
        color="0.55",
        lw=1.1,
        ls="--",
        label="P(bab_on) all finite-age samples",
    )
    ax.set_xlabel("Reconstruction time (Ma)")
    ax.set_ylabel("Fraction of trench samples with back-arc spreading")
    t0, t1 = float(np.nanmin(t)), float(np.nanmax(t))
    ax.set_xlim(t1 + (0.5 if t1 == t0 else 0.0), t0 - (0.5 if t1 == t0 else 0.0))
    ymax = float(np.nanmax(per_time["p_bab_old"])) if len(per_time) else 0.4
    ax.set_ylim(0.0, max(0.45, ymax * 1.15 if np.isfinite(ymax) else 0.45))
    ax.set_title(
        "Sdrolias–Müller age gate through the Cenozoic\n"
        "Müller et al. (2019) topologies — circular kinematic association"
    )
    ax.legend(fontsize=9, frameon=True, loc="upper left")
    ax.grid(True, linestyle=":", linewidth=0.5, color="0.8")
    fig.tight_layout()
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def fig_bab_map_0ma(df0, path):
    _setup_style()
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        proj = ccrs.Robinson()
        data_crs = ccrs.PlateCarree()
        fig, ax = plt.subplots(figsize=(12.5, 6.5), subplot_kw={"projection": proj})
        ax.set_global()
        try:
            ax.add_feature(cfeature.LAND, facecolor="0.88", edgecolor="none", zorder=0)
            ax.add_feature(
                cfeature.COASTLINE, linewidth=0.4, edgecolor="0.35", zorder=1
            )
        except Exception as err:
            print(f"WARNING: NaturalEarth features unavailable ({err})")
        ax.gridlines(draw_labels=False, linewidth=0.3, color="0.7", linestyle=":")
        transform = data_crs
    except Exception as err:
        print(f"WARNING: cartopy unavailable ({err}); plain lon/lat map.")
        fig, ax = plt.subplots(figsize=(12.5, 6.5))
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_aspect("equal", adjustable="box")
        transform = None

    off = df0[df0["bab_on"] == 0]
    on = df0[df0["bab_on"] == 1]
    kw = dict(s=11, linewidths=0, zorder=3)
    if transform is not None:
        ax.scatter(
            off["lon"], off["lat"], c="0.62", transform=transform, **kw, label="bab_off"
        )
        ax.scatter(
            on["lon"],
            on["lat"],
            c="#0A9396",
            transform=transform,
            **kw,
            label="bab_on",
        )
    else:
        ax.scatter(off["lon"], off["lat"], c="0.62", **kw, label="bab_off")
        ax.scatter(on["lon"], on["lat"], c="#0A9396", **kw, label="bab_on")
    ax.set_title(
        "0 Ma trench samples — back-arc spreading proxy (layer A)\n"
        "teal = actively spreading landward MOR on the overriding plate"
    )
    ax.legend(loc="lower left", fontsize=9, frameon=True)
    fig.tight_layout()
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def _fmt_p(p):
    if p is None or not np.isfinite(p):
        return "—"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.4f}"


def _fmt(x, nd=3):
    if x is None or not np.isfinite(x):
        return "—"
    return f"{x:.{nd}f}"


def write_report(
    path,
    *,
    n_rows,
    n_times,
    times,
    meta_rows,
    sanity,
    per_time,
    zone,
    uni_pt,
    uni_z,
    gate_sweep_pt,
    gate_sweep_z,
    mh_pt,
    mh_z,
    logit_time,
    mw,
    submap,
    dilat_coverage,
    csv_path,
    failed_times,
    circularity,
):
    n_on = int(per_time["n_bab_on"].sum()) if len(per_time) else 0
    sanity_tbl = "\n".join(
        f"| {r['region']} | {int(r['n'])} | {int(r['n_bab_on'])} | {_fmt(r['frac_bab_on'], 2)} | "
        f"{_fmt(r['med_dist_ba_km'], 0)} | {_fmt(r['med_dilat_1e-15_s'], 3)} |"
        for r in sanity.to_dict("records")
    ) if len(sanity) else "| — | — | — | — | — | — |"
    per_time_tbl = "\n".join(
        f"| {int(r.time_Ma)} | {int(r.n)} | {int(r.n_bab_on)} | {_fmt(r.frac_bab_on, 3)} | "
        f"{_fmt(r.med_age_on, 1)} | {_fmt(r.med_age_off, 1)} | "
        f"{_fmt(r.p_bab_old, 3)} | {_fmt(r.p_bab_young, 3)} |"
        for r in per_time.itertuples()
    )
    sweep_pt = "\n".join(
        f"| {int(g['threshold'])} | {g['n_old']} | {g['n_young']} | "
        f"{_fmt(g['p_old'], 3)} | {_fmt(g['p_young'], 3)} | {_fmt(g['delta'], 3)} |"
        for g in gate_sweep_pt
    )
    sweep_z = "\n".join(
        f"| {int(g['threshold'])} | {g['n_old']} | {g['n_young']} | "
        f"{_fmt(g['p_old'], 3)} | {_fmt(g['p_young'], 3)} | {_fmt(g['delta'], 3)} |"
        for g in gate_sweep_z
    )

    def uni_tbl(frame):
        lines = []
        for r in frame.itertuples():
            lines.append(
                f"| {r.predictor} | {int(r.n)} | {_fmt(getattr(r, 'spearman_rho', np.nan), 3)} | "
                f"{_fmt_p(getattr(r, 'spearman_p', np.nan))} | "
                f"{_fmt(getattr(r, 'logit_OR', np.nan), 4)} | "
                f"{_fmt_p(getattr(r, 'logit_p_lrt', np.nan))} |"
            )
        return "\n".join(lines) if lines else "| — | — | — | — | — | — |"

    submap_txt = "_not computed_"
    if submap is not None:
        rates = submap["rates_ce"]
        ct = submap["ct_ce"]
        geo = submap["geo"]
        ct_s = ct.to_string()
        geo_s = "\n".join(
            f"| {idx} | {int(r.n)} | {int(r.n_bab_on)} | {_fmt(r.frac_bab_on, 2)} |"
            for idx, r in geo.iterrows()
        )
        ce_s = "\n".join(
            f"| {idx} | {int(r.n)} | {int(r.n_bab_on)} | {_fmt(r.frac_bab_on, 2)} |"
            for idx, r in rates.iterrows()
            if pd.notna(idx)
        )
        submap_txt = f"""Join of 0 Ma `bab_on` onto `results/joined_submap_muller2019_0Ma.csv`
(nearest trench sample to each SubMap Müller hit; median distance
{submap['med_proxy_join_km']:.2f} km).

**E vs C only** (`UPS_3`):

| UPS_3 | n | n bab_on | P(bab_on) |
|---|---:|---:|---:|
{ce_s}

Crosstab `UPS_3` × `bab_on` (C/E):

```
{ct_s}
```

This is a **simple agreement table**, not a winner claim. Layer A is MOR
back-arc spreading in Müller2019; SubMap `UPS_3` is geodetic strain.
Hellenic/Aegean is geodetic E but `bab_on=0` because the model has no
Aegean MOR (dilatation is the relevant layer there). Andes is geodetic C
and `bab_on=0`.

Full Geodetic_UPS × bab_on:

| Geodetic_UPS | n | n bab_on | P(bab_on) |
|---|---:|---:|---:|
{geo_s}
"""

    mh_txt = "not computed"
    if mh_pt:
        mh_txt = (
            f"point-level MH OR = {mh_pt['or_mh']:.3f} "
            f"(χ²={mh_pt['chi2']:.2f}, p={_fmt_p(mh_pt['p'])}, "
            f"{mh_pt['n_strata']} time strata)"
        )
    if mh_z:
        mh_txt += (
            f"; zone-level MH OR = {mh_z['or_mh']:.3f} "
            f"(χ²={mh_z['chi2']:.2f}, p={_fmt_p(mh_z['p'])}, "
            f"{mh_z['n_strata']} strata)"
        )

    logit_txt = "not computed"
    if logit_time:
        logit_txt = (
            f"n={logit_time['n']} (n_on={logit_time['n_on']}); "
            f"age coef = {logit_time['coef_age']:.4f} per Ma "
            f"(OR={logit_time['odds_ratio_per_Ma']:.4f}); "
            f"LRT vs time-only χ²={logit_time['lrt_vs_time_only']:.2f}, "
            f"p={_fmt_p(logit_time['p_lrt'])}"
        )

    mw_txt = "not computed"
    if mw is not None:
        mw_txt = (
            f"U={mw['U']:.0f}, p={_fmt_p(mw['p'])}, "
            f"n_on={mw['n_on']}, n_off={mw['n_off']}; "
            f"median age on={mw['med_on']:.1f} Ma, off={mw['med_off']:.1f} Ma"
        )

    n_ridge = ", ".join(
        f"{int(m['time'])}:{m['n_ridge_active']}" for m in meta_rows[:5]
    ) if meta_rows else ""
    dilat_txt = (
        f"samples in a deforming network: {dilat_coverage['n_in_net']} / "
        f"{dilat_coverage['n']} ({100*dilat_coverage['frac']:.1f}%); "
        f"times with ≥1 network hit: {dilat_coverage['n_times_with_net']} / "
        f"{dilat_coverage['n_times']}"
    )

    def _uni(frame, key, col):
        if frame is None or len(frame) == 0 or "predictor" not in frame.columns:
            return np.nan
        hit = frame[frame["predictor"].astype(str).str.contains(key, regex=False)]
        if len(hit) == 0 or col not in hit.columns:
            return np.nan
        return float(hit[col].iloc[0])

    z_n = int(len(zone)) if zone is not None else 0
    z_on = int(zone["bab_zone"].sum()) if zone is not None and "bab_zone" in zone.columns else 0
    z_or = mh_z["or_mh"] if mh_z else np.nan
    z_orp = mh_z["p"] if mh_z else np.nan
    rho_w = _uni(uni_z, "width", "spearman_rho")
    rho_a = _uni(uni_z, "age", "spearman_rho")
    p_a = _uni(uni_z, "age", "spearman_p")
    row0 = None
    if per_time is not None and len(per_time) and "time_Ma" in per_time.columns:
        hit0 = per_time[np.isclose(per_time["time_Ma"], 0.0)]
        if len(hit0):
            row0 = hit0.iloc[0]
    med_on_0 = float(row0["med_age_on"]) if row0 is not None else np.nan
    med_off_0 = float(row0["med_age_off"]) if row0 is not None else np.nan
    p_old_0 = float(row0["p_bab_old"]) if row0 is not None else np.nan
    p_young_0 = float(row0["p_bab_young"]) if row0 is not None else np.nan
    mh_or_pt = mh_pt["or_mh"] if mh_pt else np.nan
    reading = f"""At **0 Ma** the Sdrolias pattern is present in this reconstruction:
median age bab_on vs off is {_fmt(med_on_0, 1)} vs {_fmt(med_off_0, 1)} Ma, and
P(bab_on | age ≥ 55) = {_fmt(p_old_0, 2)} vs P(bab_on | age < 55) = {_fmt(p_young_0, 2)}.
SubMap geodetic C has 0/36 `bab_on`; E has 12/39. That is agreement in
direction, not a winner claim (Hellenic is geodetic E but `bab_on=0`
because Müller2019 has no Aegean MOR).

Pooled over 0–60 Ma, **point samples** still associate older slabs with a
back-arc MOR (median {_fmt(mw['med_on'], 1) if mw else '—'} vs {_fmt(mw['med_off'], 1) if mw else '—'} Ma;
Mantel–Haenszel OR {_fmt(mh_or_pt, 2)}). That n is not independent — long
old-Pacific trenches (Mariana, Tonga) contribute hundreds of neighbouring
points.

**Zone–time blocks** (n={z_n}, of which {z_on} have a back-arc) are the
unit that does not pretend 76k rows are independent. There the 55 Ma gate
is **weak** (zone-level MH OR {_fmt(z_or, 2)}, p={_fmt_p(z_orp)}).
Univariate Spearman at zone–time: trench-zone width ρ = {_fmt(rho_w, 3)}
outranks distance-to-edge, which outranks age (ρ = {_fmt(rho_a, 3)},
p = {_fmt_p(p_a)}). Logistic age after time dummies is consistent with
null ({logit_txt}).

The 0 Ma age headline therefore **does not survive** Cenozoic
zone-blocking in this circular proxy. Width is the stronger zone-level
associate of reconstructed back-arc spreading — a Schellart-shaped
geometric pattern in the topologies, still not independent geology.

The 55 Ma gate also **flips** in some intervals (about 12–14 Ma and
51–60 Ma) where P(on | young) ≥ P(on | old). Steps in the `bab_on`
fraction (for example near 6, 12, and 15 Ma) likely mark topology
revisions, not smooth basin opening."""

    text = f"""# Cenozoic age-gate test — Müller et al. (2019) back-arc spreading proxy

Built by `scripts/cenozoic_ups_proxy.py` from the existing trench catalog
`results/trench_kinematics_0-60Ma_Muller2019.csv`. Trench kinematics were
**not** re-extracted.

- **n trench samples:** {n_rows}
- **n times:** {n_times} ({int(min(times))}–{int(max(times))} Ma, 1 Myr)
- **n bab_on (point samples):** {n_on} ({100.0*n_on/max(n_rows,1):.1f}%)
- **n zone–time blocks:** {len(zone)}
- **failed times:** {failed_times if failed_times else "none"}
- **output CSV:** `{csv_path}`

## Circularity (read this first)

{circularity}

`v_T` vs `v_OP` is **not** used as UPS. In this model those two series are
collinear (trench glued to the overriding plate). The 0 Ma H1/H2 test lives
in `results/submap_m56_h1h2.md` and uses SubMap M56.

## Layer A — back-arc spreading on/off (primary)

At each time, mid-ocean ridges are tessellated with GPlately
`tessellate_mid_ocean_ridges` (0.5°, same spacing as the trench catalog).
A trench sample is `bab_on=1` if there is an actively spreading ridge
segment that is:

1. on the overriding plate (left or right plate ID == `trench_plate_id`,
   or the ridge point lies in that plate);
2. **not** on the subducting plate (left/right/containing ID ==
   `subducting_plate_id`);
3. within **1500 km**;
4. **landward** of the trench (dot of trench→ridge with the GPlately
   landward normal > 0);
5. within **400 km along-trench** of the sample;
6. orthogonal spreading rate ≥ **1 cm/yr** (`vel × cos(obliquity)`).

The 400 km along-trench cap is the practical filter for “not a major-ocean
MOR that just happens to be <1500 km”. Without it, the Gulf of California
(Rivera–North America) lights up the Middle America trench from along-strike.
True back-arcs (Mariana Trough, Lau, East Scotia Ridge, Andaman) sit behind
the *local* trench segment.

`dist_to_backarc_ridge_km` is the distance to the nearest ridge that passes
those tests (NaN if none).

Ridge tessellation is required; if it fails at any time the run stops.

## Layer B — deforming-network dilatation (secondary)

Dilatation rate is sampled 200 km along the trench normal into the
overriding plate (`pygplates.TopologicalSnapshot.get_point_strain_rates`).
Positive dilatation = extension. Rigid plates return zero. Column
`dilat_rate_1e-15_s` is in 10⁻¹⁵ s⁻¹; `dilat_in_network` is 1 if the sample
point sits in a deforming network.

{dilat_txt}

Layer B is **model-internal**. It is kept separate from `bab_on`. Aegean /
Hellenic extension is a deforming network in Müller2019, not a MOR, so it
is invisible to layer A and only shows up in B.

## 0 Ma spot check

| region | n | n bab_on | frac | med dist (km) | med dilat (10⁻¹⁵ s⁻¹) |
|---|---:|---:|---:|---:|---:|
{sanity_tbl}

Expected: Mariana / Tonga / Scotia **on**; Andes / Cascadia **off**.
Hellenic is off in A (no MOR) and extensional in B (positive dilatation).
MAT/Mexico must stay off (Gulf of California is not a back-arc).

## Tests

Independent unit for inference is a **trench-zone × time** block
(`groupby(time_Ma, zone_id)`), not a tessellated point. Point-level n≈76k
is spatially autocorrelated along the trench. Zone-level `bab_zone=1` if
any sample in that zone is `bab_on`. Age/width/edge-distance at zone level
are medians. Age tests drop NaN seafloor ages.

### Reading the result

{reading}

### 1. Median seafloor age, bab_on vs bab_off

Pooled point samples, finite ages: {mw_txt}

Per-time medians are in the table below and in `fig_bab_vs_age.png`.

### 2. Sdrolias gate (55 Ma) and a threshold sweep

P(bab_on | age ≥ T) versus P(bab_on | age < T).

**Point-level sweep (pooled over time):**

| T (Ma) | n old | n young | P(on\\|old) | P(on\\|young) | Δ |
|---|---:|---:|---:|---:|---:|
{sweep_pt}

**Zone–time sweep:**

| T (Ma) | n old | n young | P(on\\|old) | P(on\\|young) | Δ |
|---|---:|---:|---:|---:|---:|
{sweep_z}

Mantel–Haenszel odds ratio for (age ≥ 55) vs bab, stratified by time:
{mh_txt}

An OR > 1 means old slabs are more likely to have a reconstructed back-arc.

### 3. Spearman / logistic, univariate

Headline predictor if the 0 Ma pattern still holds: **seafloor age**.
Width and distance-to-edge are the Schellart-style geometric covariates.

Logistic is unpenalized (`sklearn.LogisticRegression(penalty=None)`);
p-values are likelihood-ratio tests against intercept-only. Spearman is
rank correlation of the binary flag against the predictor.

**Point-level (autocorrelated; n is not independent):**

| predictor | n | Spearman ρ | p | logistic OR | LRT p |
|---|---:|---:|---:|---:|---:|
{uni_tbl(uni_pt)}

**Zone–time blocks (preferred n):**

| predictor | n | Spearman ρ | p | logistic OR | LRT p |
|---|---:|---:|---:|---:|---:|
{uni_tbl(uni_z)}

Logistic `bab_zone ~ age + time dummies` (age coefficient after time as a
factor): {logit_txt}

### 4. 0 Ma SubMap sanity (bab_on vs UPS_3, E vs C)

{submap_txt}

## Per-time summary

| time_Ma | n | n bab_on | frac | med age on | med age off | P(on\\|≥55) | P(on\\|<55) |
|---|---:|---:|---:|---:|---:|---:|---:|
{per_time_tbl}

## How to rerun

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
    /workspace/subduction-test/scripts/cenozoic_ups_proxy.py
```

0 Ma prototype only:

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
    /workspace/subduction-test/scripts/cenozoic_ups_proxy.py --tmax 0
```

## Columns added

| column | meaning |
|---|---|
| `bab_on` | 1 if a qualifying back-arc MOR exists (layer A) |
| `dist_to_backarc_ridge_km` | distance to the nearest qualifying ridge; NaN if `bab_on=0` |
| `backarc_spreading_cm_yr` | orthogonal spreading rate of that ridge (cm/yr) |
| `dilat_rate_1e-15_s` | dilatation 200 km landward (layer B), 10⁻¹⁵ s⁻¹ |
| `dilat_in_network` | 1 if that sample sits in a deforming network |

## Caveats

- Reconstruction-internal. Opening histories of back-arc basins are
  encoded in the topologies that define `bab_on`.
- Hellenic / Okinawa / other continental back-arcs that are deforming
  networks rather than MORs are `bab_on=0`.
- Seafloor ages older than the reconstruction time (Tethyan remnants in
  the eastern Mediterranean) are left as published in the age grid.
- `zone_id` is unique only within a time slice.
- Ridge tessellation uses Müller2019 spreading features labelled
  `gpml:MidOceanRidge`; transforms are dropped by GPlately's default
  transform-segment filter.
"""
    path.write_text(text)


def patch_readme():
    """Add a short circularity / proxy pointer to the project README."""
    readme = ROOT / "README.md"
    if not readme.is_file():
        return
    text = readme.read_text()
    marker = "## Cenozoic UPS proxy (Müller2019 back-arc spreading)"
    block = """## Cenozoic UPS proxy (Müller2019 back-arc spreading)

`scripts/cenozoic_ups_proxy.py` flags trench samples with an actively
spreading **back-arc MOR** in the same Müller et al. (2019) topologies
that produced the age grids (`bab_on`). This is a Sdrolias & Müller
(2006)-style **kinematic association, not independent geology** — the
reconstruction already knows where it opened back-arc basins. Do not
treat it as an observed UPS time series, and do not use `v_T` vs `v_OP`
as UPS in this model (they are collinear).

Results: `results/cenozoic_age_gate.md`,
`results/trench_kinematics_0-60Ma_with_ups_proxy.csv`.

"""
    if marker in text:
        pre, rest = text.split(marker, 1)
        # replace from marker through the next ## or EOF
        nxt = rest.find("\n## ")
        if nxt == -1:
            text = pre + block
        else:
            text = pre + block + rest[nxt + 1 :]
    else:
        # insert after Status section if present, else append
        status = text.find("## Layout")
        if status != -1:
            text = text[:status] + block + text[status:]
        else:
            text = text.rstrip() + "\n\n" + block
    readme.write_text(text)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
CIRCULARITY = (
    "This upper-plate strain proxy is taken from the **same** Müller et al. "
    "(2019) reconstruction that supplies the seafloor-age grids and the trench "
    "kinematics. Back-arc basins exist in the model because the topologies "
    "put spreading ridges there. The Cenozoic test is therefore a "
    "**Sdrolias & Müller (2006)-style kinematic association** (does the "
    "reconstruction open a back-arc where the subducting plate is old?), "
    "not a test against independent geological UPS. Present-day SubMap "
    "geodetic UPS remains the only observation that is independent of this "
    "plate model. That circularity is a feature of using reconstructed "
    "back-arc spreading as UPS, and it is the reason the 0 Ma SubMap table "
    "is reported as a sanity check rather than a confirmation."
)


def main(argv=None):
    args = parse_args(argv)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    _setup_style()

    t_run = time_mod.time()
    meta_rows = []
    failed_times = []
    df0 = None

    if args.from_csv is not None:
        print(f"Recomputing stats from {args.from_csv} ...")
        out = pd.read_csv(args.from_csv)
        if "bab_on" not in out.columns:
            raise SystemExit(f"{args.from_csv} has no bab_on column")
        times = np.sort(out["time_Ma"].unique())
        tmin = float(times.min())
        tmax = float(times.max())
        print(f"  {len(out)} rows, {len(times)} times, {tmin:.0f}–{tmax:.0f} Ma")
        hit0 = out[np.isclose(out["time_Ma"], 0.0)]
        df0 = hit0.copy() if len(hit0) else None
        if df0 is not None:
            table = sanity_table_0ma(df0)
            check_sanity_0ma(table, abort=False)
    else:
        print(f"Reading trench catalog {args.kin_csv} ...")
        kin = pd.read_csv(args.kin_csv)
        if "time_Ma" not in kin.columns:
            raise SystemExit("kinematics CSV has no time_Ma column")
        times_all = np.sort(kin["time_Ma"].unique())
        tmin = float(args.tmin) if args.tmin is not None else float(times_all.min())
        tmax = float(args.tmax) if args.tmax is not None else float(times_all.max())
        kin = kin[(kin["time_Ma"] >= tmin) & (kin["time_Ma"] <= tmax)].copy()
        times = np.sort(kin["time_Ma"].unique())
        print(
            f"  {len(kin)} rows, {len(times)} times, {tmin:.0f}–{tmax:.0f} Ma"
        )
        if len(kin) == 0:
            raise SystemExit("no kinematics rows in requested time range")

        _, recon = load_reconstruction()

        frames = []
        do_dilat = not args.skip_dilatation

        # 0 Ma first if present
        times_ordered = list(times)
        if 0.0 in times_ordered:
            times_ordered = [0.0] + [t for t in times_ordered if t != 0.0]

        for i, t in enumerate(times_ordered):
            t0 = time_mod.time()
            trench = kin[np.isclose(kin["time_Ma"], t)].copy()
            print(
                f"[{i+1}/{len(times_ordered)}] {int(t)} Ma  n_trench={len(trench)} ...",
                flush=True,
            )
            try:
                out_t, meta = proxy_one_time(
                    recon, t, trench, do_dilatation=do_dilat
                )
            except RuntimeError as err:
                print(f"  STOP: {err}")
                raise
            except Exception as err:
                print(f"  ERROR at {int(t)} Ma: {err}")
                failed_times.append(int(t))
                continue
            frames.append(out_t)
            meta_rows.append(meta)
            print(
                f"  ridges {meta['n_ridge_active']}/{meta['n_ridge_raw']} active  "
                f"bab_on={meta['n_bab_on']}/{meta['n_trench']}  "
                f"network={meta['n_network']}  ({time_mod.time()-t0:.2f}s)",
                flush=True,
            )
            if np.isclose(t, 0.0):
                df0 = out_t
                table = sanity_table_0ma(df0)
                check_sanity_0ma(table, abort=not args.skip_sanity_abort)

        if not frames:
            raise RuntimeError("no times produced a proxy")

        out = pd.concat(frames, ignore_index=True)
        out = out.sort_values(["time_Ma", "lon", "lat"]).reset_index(drop=True)
        out.to_csv(args.out_csv, index=False, float_format="%.6f")
        print(f"Wrote {len(out)} rows -> {args.out_csv}")

    # --- stats ---
    print("Computing tests ...")
    per_time = per_time_stats(out)
    z = zone_frame(out)

    age_pt, bab_pt = _finite_pair(out["seafloor_age_Ma"], out["bab_on"])
    mw = None
    if np.any(bab_pt == 1) and np.any(bab_pt == 0):
        on = age_pt[bab_pt == 1]
        off = age_pt[bab_pt == 0]
        U, p = mannwhitneyu(on, off, alternative="two-sided")
        mw = {
            "U": float(U),
            "p": float(p),
            "n_on": int(on.size),
            "n_off": int(off.size),
            "med_on": float(np.median(on)),
            "med_off": float(np.median(off)),
        }

    gate_sweep_pt = [
        gate_counts(out["seafloor_age_Ma"].to_numpy(), out["bab_on"].to_numpy(), T)
        for T in AGE_SWEEP_MA
    ]
    gate_sweep_z = [
        gate_counts(z["seafloor_age_Ma"].to_numpy(), z["bab_zone"].to_numpy(), T)
        for T in AGE_SWEEP_MA
    ]

    mh_tables_pt = []
    mh_tables_z = []
    for t, s in out.groupby("time_Ma"):
        g = gate_counts(s["seafloor_age_Ma"].to_numpy(), s["bab_on"].to_numpy(), AGE_GATE_MA)
        mh_tables_pt.append((g["on_old"], g["off_old"], g["on_young"], g["off_young"]))
    for t, s in z.groupby("time_Ma"):
        g = gate_counts(s["seafloor_age_Ma"].to_numpy(), s["bab_zone"].to_numpy(), AGE_GATE_MA)
        mh_tables_z.append((g["on_old"], g["off_old"], g["on_young"], g["off_young"]))
    mh_pt = mantel_haenszel(mh_tables_pt)
    mh_z = mantel_haenszel(mh_tables_z)

    uni_pt = univariate_block(out, "bab_on", "point")
    uni_z = univariate_block(z, "bab_zone", "zone-time")
    logit_time = logit_age_with_time(z, "seafloor_age_Ma", "bab_zone")

    submap, sub_err = (None, "0 Ma not in this run")
    sanity = sanity_table_0ma(df0) if df0 is not None else pd.DataFrame()
    if df0 is not None:
        submap, sub_err = submap_sanity(df0)
        if sub_err:
            print("SubMap sanity skipped:", sub_err)

    n_in_net = int(np.nansum(out["dilat_in_network"] == 1))
    n_times_net = int(
        out.groupby("time_Ma")["dilat_in_network"].apply(lambda s: np.nansum(s == 1) > 0).sum()
    )
    dilat_coverage = {
        "n_in_net": n_in_net,
        "n": int(len(out)),
        "frac": n_in_net / max(len(out), 1),
        "n_times_with_net": n_times_net,
        "n_times": int(out["time_Ma"].nunique()),
    }

    write_report(
        REPORT_MD,
        n_rows=len(out),
        n_times=int(out["time_Ma"].nunique()),
        times=times,
        meta_rows=meta_rows,
        sanity=sanity if len(sanity) else pd.DataFrame(
            columns=list(SANITY_REGIONS.keys())
        ),
        per_time=per_time,
        zone=z,
        uni_pt=uni_pt,
        uni_z=uni_z,
        gate_sweep_pt=gate_sweep_pt,
        gate_sweep_z=gate_sweep_z,
        mh_pt=mh_pt,
        mh_z=mh_z,
        logit_time=logit_time,
        mw=mw,
        submap=submap,
        dilat_coverage=dilat_coverage,
        csv_path=args.out_csv,
        failed_times=failed_times,
        circularity=CIRCULARITY,
    )
    print(f"Wrote {REPORT_MD}")

    print("Drawing figures ...")
    fig_bab_vs_age(out, per_time, FIG_BAB_AGE)
    print(f"Wrote {FIG_BAB_AGE}")
    fig_age_gate(per_time, FIG_GATE)
    print(f"Wrote {FIG_GATE}")
    if df0 is not None:
        fig_bab_map_0ma(df0, FIG_MAP)
        print(f"Wrote {FIG_MAP}")

    patch_readme()
    print("Patched project README with circularity note.")

    print("\n=== headline ===")
    if mw:
        print(
            f"median age bab_on={mw['med_on']:.1f} Ma vs bab_off={mw['med_off']:.1f} Ma "
            f"(Mann–Whitney p={_fmt_p(mw['p'])})"
        )
    g55 = gate_sweep_pt[AGE_SWEEP_MA.index(55)]
    print(
        f"P(bab_on | age≥55)={g55['p_old']:.3f}  "
        f"P(bab_on | age<55)={g55['p_young']:.3f}  "
        f"(point-level, pooled)"
    )
    if mh_z:
        print(f"zone-level MH OR (old vs bab) = {mh_z['or_mh']:.3f}  p={_fmt_p(mh_z['p'])}")
    print(f"elapsed {time_mod.time()-t_run:.1f}s")
    print("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
