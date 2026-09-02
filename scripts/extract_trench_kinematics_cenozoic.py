#!/usr/bin/env python3
"""Extract Cenozoic (0–60 Ma) trench kinematics from Müller et al. (2019).

Reuses helpers from ``extract_trench_kinematics_0Ma.py``. This product is
trench kinematics only: it does **not** join SubMap / UPS (present-day).

Predictors extracted because they are independent of the Müller2019
``v_T_perp ≈ −v_OP_perp`` collinearity (trench glued to the overriding plate):

- seafloor age at trench (time-dependent Müller2019 age grids)
- trench-zone width
- distance to nearest trench edge
- convergence rate / obliquity
- lon, lat, time, plate IDs

``v_T_*`` and ``v_OP_*`` are retained for continuity with the 0 Ma product and
to document the collinearity. They are **not** independent rollback predictors
in this model and must not be used to test H1/H2 (that was SubMap M56 at 0 Ma).

Run:
    /workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
        /workspace/subduction-test/scripts/extract_trench_kinematics_cenozoic.py

Optional:
    .../extract_trench_kinematics_cenozoic.py --dt 5
"""

from __future__ import annotations

import argparse
import sys
import time as time_mod
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gplately import PlateModelManager, PlateReconstruction
from gplately.grids import Raster

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import extract_trench_kinematics_0Ma as k0  # noqa: E402

# ---------------------------------------------------------------------------
# Paths / defaults
# ---------------------------------------------------------------------------
DATA_DIR = k0.DATA_DIR
RESULTS_DIR = k0.RESULTS_DIR
CSV_PATH = RESULTS_DIR / "trench_kinematics_0-60Ma_Muller2019.csv"
README_PATH = RESULTS_DIR / "README_kinematics_cenozoic.md"
FIG_PATH = RESULTS_DIR / "fig_age_at_trench_vs_time.png"

T_MIN_MA = 0
T_MAX_MA = 60
DEFAULT_DT_MYR = 1  # 1 Myr steps; pass --dt 5 if downloads are too slow
AGE_DOWNLOAD_BATCH = 8


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Cenozoic trench kinematics from Müller et al. (2019)."
    )
    p.add_argument("--tmin", type=int, default=T_MIN_MA, help="youngest time (Ma)")
    p.add_argument("--tmax", type=int, default=T_MAX_MA, help="oldest time (Ma)")
    p.add_argument(
        "--dt",
        type=int,
        default=DEFAULT_DT_MYR,
        help="time step in Myr (1 is default; 5 if age-grid download is too large)",
    )
    p.add_argument(
        "--out-csv",
        type=Path,
        default=CSV_PATH,
        help="output CSV path",
    )
    return p.parse_args(argv)


def times_ma(tmin, tmax, dt):
    if dt <= 0:
        raise ValueError("--dt must be a positive integer")
    if tmax < tmin:
        raise ValueError("--tmax must be >= --tmin")
    t = np.arange(int(tmin), int(tmax) + 1, int(dt), dtype=int)
    if t.size == 0 or t[-1] != int(tmax):
        # Always include tmax even if it is not on the dt grid.
        t = np.unique(np.append(t, int(tmax)))
    return t


def local_age_grid_path(model, time):
    """Return the expected local path for a Müller2019 age grid, or None."""
    cfg = model.get_cfg()
    template = cfg["TimeDepRasters"]["AgeGrids"]
    fname = template.format(time).split("/")[-1]
    path = Path(model.get_model_dir()) / "Rasters" / "AgeGrids" / fname
    if path.is_file() and path.stat().st_size > 1_000_000:
        return path
    return None


def ensure_age_grids(model, times, batch_size=AGE_DOWNLOAD_BATCH):
    """Download missing Müller2019 age grids; return {time: path}."""
    paths = {}
    missing = []
    for t in times:
        local = local_age_grid_path(model, int(t))
        if local is not None:
            paths[int(t)] = str(local)
        else:
            missing.append(int(t))
    print(
        f"Age grids on disk: {len(paths)}/{len(times)}; "
        f"need to download {len(missing)}"
    )
    if not missing:
        return paths

    for i in range(0, len(missing), batch_size):
        batch = missing[i : i + batch_size]
        t0 = time_mod.time()
        print(f"  downloading AgeGrids Ma={list(batch)} ...", flush=True)
        got = model.get_age_grids(batch)
        for t, p in zip(batch, got):
            if not Path(p).is_file():
                raise FileNotFoundError(f"Age grid download failed for {t} Ma: {p}")
            paths[int(t)] = p
        print(
            f"    got {len(batch)} grids in {time_mod.time() - t0:.1f}s "
            f"({Path(got[0]).name} ...)",
            flush=True,
        )
    return paths


def _col(gdf, *names):
    for name in names:
        if name in gdf.columns:
            return gdf[name].to_numpy()
    return None


def extract_one_time(recon, time, tess_rad, extra_names, age_path):
    """Tessellate and derive kinematics at one reconstruction time.

    Column mapping and derived-field logic follow
    ``extract_trench_kinematics_0Ma.py``. Upper-plate nature is omitted:
    the 0 Ma classifier uses present-day continental polygons.
    """
    gdf, extras_used, extra_err = k0.tessellate_with_fallback(
        recon, float(time), tess_rad, extra_names
    )
    if gdf is None or len(gdf) == 0:
        return (
            pd.DataFrame(),
            extras_used,
            extra_err,
            ["empty tessellation"],
        )

    lon = (
        gdf.geometry.x.to_numpy()
        if "geometry" in gdf.columns
        else gdf["lon"].to_numpy()
    )
    lat = (
        gdf.geometry.y.to_numpy()
        if "geometry" in gdf.columns
        else gdf["lat"].to_numpy()
    )

    sub_id = _col(gdf, "subducting plate ID")
    trench_id = _col(gdf, "trench plate ID")
    conv_mag = _col(gdf, "convergence velocity (cm/yr)")
    conv_obliq = _col(gdf, "convergence obliquity angle (degrees)")
    trench_az = _col(gdf, "trench normal angle (degrees)")
    arc_deg = _col(gdf, "length (degrees)")
    conv_orth = _col(gdf, "convergence velocity orthogonal component (cm/yr)")
    conv_par = _col(gdf, "convergence velocity parallel component (cm/yr)")
    trench_orth = _col(gdf, "trench absolute velocity orthogonal component (cm/yr)")
    trench_par = _col(gdf, "trench absolute velocity parallel component (cm/yr)")
    dist_edge_deg = _col(gdf, "distance to nearest trench edge (degrees)")

    missing = []

    if conv_orth is None or conv_par is None:
        if conv_mag is not None and conv_obliq is not None:
            rad = np.radians(conv_obliq)
            conv_orth = np.cos(rad) * np.abs(conv_mag)
            conv_par = np.sin(rad) * np.abs(conv_mag)
        else:
            missing.append("convergence trench-normal / trench-parallel")

    v_T_perp = None
    v_T_par = None
    if trench_orth is None:
        trench_abs_mag = _col(gdf, "trench velocity (cm/yr)")
        trench_abs_obliq = _col(gdf, "trench obliquity angle (degrees)")
        if trench_abs_mag is not None and trench_abs_obliq is not None:
            rad = np.radians(trench_abs_obliq)
            trench_orth = np.cos(rad) * np.abs(trench_abs_mag)
            trench_par = np.sin(rad) * np.abs(trench_abs_mag)
        else:
            missing.append("v_T_perp / v_T_parallel")
    if trench_orth is not None:
        # Flip GPlately trench-absolute orthogonal so oceanward retreat > 0.
        v_T_perp = -np.asarray(trench_orth, dtype=float)
        v_T_par = np.asarray(trench_par, dtype=float)

    v_east, v_north = k0.plate_abs_velocity_east_north(
        recon.rotation_model,
        lon,
        lat,
        trench_id,
        float(time),
        k0.VELOCITY_DELTA_TIME_MYR,
    )
    v_OP_perp, v_OP_par = k0.trench_normal_components(v_east, v_north, trench_az)

    # Seafloor age: prefer 50 km oceanward of the trench (toward the
    # subducting plate), fall back to the trench node. NaN where the age
    # grid has no seafloor.
    age_raster = Raster(str(age_path))
    age_trench = np.asarray(
        age_raster.interpolate(lon, lat, method="linear"), dtype=float
    )
    lon_sp, lat_sp = k0._destination_lonlat(
        lon, lat, trench_az + 180.0, k0.OFFSET_INTO_PLATE_KM
    )
    age_sp = np.asarray(
        age_raster.interpolate(lon_sp, lat_sp, method="linear"), dtype=float
    )
    seafloor_age = np.where(np.isfinite(age_sp), age_sp, age_trench)

    if dist_edge_deg is None:
        missing.append("distance to nearest trench/slab edge")
        dist_edge_km = np.full(len(lon), np.nan)
    else:
        dist_edge_km = np.asarray(dist_edge_deg, dtype=float) * k0.DEG_TO_KM

    arc_km = np.asarray(arc_deg, dtype=float) * k0.DEG_TO_KM
    tess_km = k0.TESSELLATION_DEG * k0.DEG_TO_KM
    zone_id, zone_width_km = k0.connected_zone_width_km(
        lon,
        lat,
        sub_id.astype(int),
        trench_id.astype(int),
        arc_km,
        tess_km,
    )

    n = len(lon)
    out = pd.DataFrame(
        {
            "time_Ma": np.full(n, float(time)),
            "lon": lon,
            "lat": lat,
            "subducting_plate_id": sub_id.astype(int),
            "trench_plate_id": trench_id.astype(int),
            "trench_normal_azimuth_deg": trench_az,
            "arc_segment_length_km": arc_km,
            "conv_mag_cm_yr": conv_mag,
            "conv_obliquity_deg": conv_obliq,
            "conv_perp_cm_yr": conv_orth,
            "conv_parallel_cm_yr": conv_par,
            "v_T_perp": v_T_perp if v_T_perp is not None else np.full(n, np.nan),
            "v_T_parallel": v_T_par if v_T_par is not None else np.full(n, np.nan),
            "v_OP_perp": v_OP_perp,
            "v_OP_parallel": v_OP_par,
            "seafloor_age_Ma": seafloor_age,
            "dist_nearest_trench_edge_km": dist_edge_km,
            "zone_id": zone_id,
            "trench_zone_width_km": zone_width_km,
        }
    )
    if dist_edge_deg is None:
        out = out.drop(columns=["dist_nearest_trench_edge_km"])
    return out, extras_used, extra_err, missing


def make_age_vs_time_figure(df, path):
    """Median seafloor age at trench vs reconstruction time, with IQR band."""
    g = df.groupby("time_Ma")["seafloor_age_Ma"]
    times = g.median().index.to_numpy(dtype=float)
    med = g.median().to_numpy(dtype=float)
    q10 = g.quantile(0.10).to_numpy(dtype=float)
    q25 = g.quantile(0.25).to_numpy(dtype=float)
    q75 = g.quantile(0.75).to_numpy(dtype=float)
    q90 = g.quantile(0.90).to_numpy(dtype=float)
    n_valid = g.count().reindex(times).to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    ax.fill_between(
        times, q10, q90, color="0.82", linewidth=0, zorder=1, label="10–90 percentile"
    )
    ax.fill_between(
        times, q25, q75, color="0.58", linewidth=0, zorder=2, label="IQR (25–75)"
    )
    ax.plot(times, med, color="0.08", lw=2.1, zorder=3, label="median")
    ax.set_xlabel("Reconstruction time (Ma)")
    ax.set_ylabel("Seafloor age at trench (Ma)")
    tmin = float(np.nanmin(times))
    tmax = float(np.nanmax(times))
    ax.set_xlim(tmax, tmin)  # geological: present (0 Ma) on the right
    ymax = np.nanmax(q90)
    ax.set_ylim(0.0, max(50.0, float(ymax) * 1.08) if np.isfinite(ymax) else 200.0)
    ax.grid(True, linestyle=":", linewidth=0.6, color="0.75", zorder=0)
    ax.set_title(
        "Müller et al. (2019) — seafloor age at trench vs time\n"
        "median of tessellated trench samples (NaNs on continents excluded)"
    )
    ax.legend(loc="upper left", frameon=True, fontsize=9)

    ax2 = ax.twinx()
    ax2.plot(
        times,
        n_valid,
        color="C0",
        lw=1.0,
        ls="--",
        alpha=0.7,
        zorder=4,
        label="n finite ages",
    )
    ax2.set_ylabel("n samples with finite seafloor age", color="C0")
    ax2.tick_params(axis="y", labelcolor="C0")
    ax2.set_ylim(0.0, max(float(np.nanmax(n_valid)) * 1.25, 10.0))
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", frameon=True, fontsize=9)

    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def _fmt_stats(stats):
    if stats is None:
        return "not computed"
    vmin, vmed, vmax, nvalid = stats
    return f"min={vmin:.3f}, median={vmed:.3f}, max={vmax:.3f} (n_valid={nvalid})"


def write_readme(
    path,
    *,
    n_rows,
    n_times,
    tmin,
    tmax,
    dt,
    times,
    columns,
    extras_used,
    extras_rejected,
    n_per_time,
    age_stats,
    conv_stats,
    vt_stats,
    vop_stats,
    missing_fields,
    age_grid_count,
    csv_path,
    fig_path,
    failed_times,
    collinearity_note,
):
    time_list = ", ".join(str(int(t)) for t in times)
    n_table = "\n".join(
        f"| {int(t)} | {n_per_time.get(float(t), n_per_time.get(int(t), 0))} |"
        for t in times
    )
    failed_txt = (
        ", ".join(str(t) for t in failed_times) if failed_times else "none"
    )
    text = f"""# Cenozoic trench kinematics — Müller et al. (2019)

Extracted by `scripts/extract_trench_kinematics_cenozoic.py`.
Reuses helpers from `scripts/extract_trench_kinematics_0Ma.py`.

This product is **trench kinematics only**. It is **not** joined to SubMap or
UPS. SubMap is present-day only; the 0 Ma H1 vs H2 comparison lives in
`results/submap_m56_h1h2.md` and must not be repeated from Müller2019
`v_T` vs `v_OP`.

- **n rows:** {n_rows}
- **n times:** {n_times}
- **time range:** {tmin}–{tmax} Ma
- **time step:** {dt} Myr
- **times (Ma):** {time_list}
- **tessellation:** {k0.TESSELLATION_DEG} deg (~{k0.TESSELLATION_DEG * k0.DEG_TO_KM:.1f} km)
- **velocity delta time:** {k0.VELOCITY_DELTA_TIME_MYR} Myr
- **Earth radius:** {k0.EARTH_RADIUS_KM} km (`pygplates.Earth.mean_radius_in_kms`)
- **age grids on disk for this run:** {age_grid_count}
- **GPlately extras accepted:** {', '.join(extras_used) if extras_used else '(none)'}
- **GPlately extras rejected:** {', '.join(extras_rejected) if extras_rejected else '(none)'}
- **failed times:** {failed_txt}

## Why these predictors (and not rollback from v_T)

In Müller et al. (2019) the trench polyline is typically attached to
`trench_plate_id`, so **`v_T_perp ≈ −v_OP_perp`** (trench glued to the
overriding plate). That collinearity is a model property, not a geodynamic
measurement of independent hinge rollback. Do **not** invent independent
trench rollback from this model, and do **not** claim H1/H2 from `v_T` vs
`v_OP` here.

Predictors that **are** independent of that collinearity, and that this CSV
is meant to carry through the Cenozoic:

| predictor | column |
|---|---|
| seafloor age at trench | `seafloor_age_Ma` |
| trench-zone width | `trench_zone_width_km` |
| distance to nearest trench edge | `dist_nearest_trench_edge_km` |
| convergence rate / obliquity | `conv_mag_cm_yr`, `conv_obliquity_deg`, `conv_perp_cm_yr`, `conv_parallel_cm_yr` |
| location, time, plate IDs | `lon`, `lat`, `time_Ma`, `subducting_plate_id`, `trench_plate_id` |

`v_T_*` and `v_OP_*` are written for continuity with the 0 Ma product and so
the collinearity can be checked. They are **not** independent predictors in
this model. {collinearity_note}

`upper_plate_nature` is **not** computed: the 0 Ma classifier uses present-day
continental polygons, which are not valid through time without reconstruction.

## Model citation

Müller, R. D., Zahirovic, S., Williams, S. E., Cannon, J., Seton, M., Bower, D. J.,
Tetley, M. G., Heine, C., Le Breton, E., Liu, S., Russell, S. H. J., Yang, T.,
Leonard, J., & Gurnis, M. (2019). A global plate model including lithospheric
deformation along major rifts and orogens since the Triassic. *Tectonics*,
38(6), 1884–1907. https://doi.org/10.1029/2018TC005462

Digital plate-model package (GPlately / GPlates):
https://doi.org/10.5281/zenodo.10525286
(version recorded in local metadata: `10.5281/zenodo.11601026`).

Seafloor age grids: Müller et al. (2019) Tectonics v2.0 netCDF age grids,
1 Myr snapshots, downloaded via PlateModelManager `AgeGrids` from
`https://repo.gplates.org/webdav/PlateModel_Age_SR_Grids/Muller_etal_2019_Tectonics/v2.0/`.

Local data directory: `/workspace/subduction-test/data/Muller2019`
(PlateModelManager also resolves the lowercase symlink `data/muller2019`).
Age grids live in `data/Muller2019/Rasters/AgeGrids/`.

## How to rerun

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
    /workspace/subduction-test/scripts/extract_trench_kinematics_cenozoic.py
```

Equivalent micromamba one-shot:

```bash
export MAMBA_ROOT_PREFIX=/workspace/subduction-test/.micromamba
/workspace/subduction-test/.micromamba/bin/micromamba run -n pygplates \\
    --root-prefix "$MAMBA_ROOT_PREFIX" python \\
    /workspace/subduction-test/scripts/extract_trench_kinematics_cenozoic.py
```

Coarser time step (downloads 13 age grids instead of 61):

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
    /workspace/subduction-test/scripts/extract_trench_kinematics_cenozoic.py --dt 5
```

Already-downloaded age grids are reused; missing grids are fetched in batches.

Outputs:

- `{csv_path}`
- `{fig_path}`
- this README

## Sign conventions (same as the 0 Ma product)

GPlately's trench **normal always points toward the overriding plate**.
Obliquity is the signed angle between that normal and the velocity vector,
range (−180, 180) deg: positive is clockwise from the trench normal when
viewed from above. Parallel components are positive 90 deg clockwise from
the trench normal.

- **convergence orthogonal > 0** = subducting plate moving toward the
  overriding plate (convergence). Convergence magnitude is negative if the
  plates are diverging (|obliquity| > 90).
- **`v_T_perp` > 0** = oceanward trench retreat. This is GPlately's trench-
  absolute orthogonal **flipped** (GPlately orth > 0 is landward advance).
- **`v_OP_perp` > 0** = overriding plate moving away from the trench
  (absolute velocity of `trench_plate_id` at the sample).

## Column dictionary

{chr(10).join(f'- `{c}`' for c in columns)}

| column | units | meaning |
|---|---|---|
| time_Ma | Ma | reconstruction time of this sample |
| lon, lat | deg | tessellated trench sample |
| subducting_plate_id | int | GPlately subducting plate ID |
| trench_plate_id | int | GPlately trench plate ID (used as overriding plate) |
| trench_normal_azimuth_deg | deg | clockwise from North, toward overriding plate |
| arc_segment_length_km | km | length of the tessellated arc segment this sample sits on |
| conv_mag_cm_yr | cm/yr | GPlately signed convergence magnitude (negative if diverging) |
| conv_obliquity_deg | deg | angle trench-normal → convergence velocity |
| conv_perp_cm_yr | cm/yr | trench-normal convergence; + = converging |
| conv_parallel_cm_yr | cm/yr | trench-parallel convergence component |
| v_T_perp | cm/yr | trench-normal trench abs. vel.; + = oceanward retreat. **Collinear with −v_OP_perp in this model — not an independent rollback predictor.** |
| v_T_parallel | cm/yr | trench-parallel trench abs. vel. |
| v_OP_perp | cm/yr | trench-normal OP abs. vel.; + = OP moving away from trench |
| v_OP_parallel | cm/yr | trench-parallel overriding-plate abs. vel. |
| seafloor_age_Ma | Myr | age-grid sample at this reconstruction time; NaN on continents / missing |
| dist_nearest_trench_edge_km | km | distance along trench to nearest slab/trench edge |
| zone_id | int | spatially connected polyline id **within this time slice** (same subducting plate ID). Not unique across times; group by `(time_Ma, zone_id)`. |
| trench_zone_width_km | km | along-trench length of that connected polyline |

## Sampling notes

- **Seafloor age** is sampled from the Müller2019 age grid at the same
  reconstruction time as the trench, preferring a point 50 km oceanward of
  the trench (toward the subducting plate) so Andean-type trenches that sit
  on the continental side of the raster still receive oceanic lithosphere
  age. The trench node is the fallback. NaN remains where the raster has no
  seafloor. Source-grid values older than the reconstruction time (e.g.
  Tethyan remnants) are left as published; they are not clipped or filled.
- **`trench_zone_width_km`** is the along-trench length of each spatially
  connected set of samples that share a subducting plate ID (neighbour
  gap ≤ 2.5 × tessellation). Disconnected trenches of the same plate stay
  separate. Zone IDs are recomputed independently at each time.
- Kinematic values are taken from GPlately / the rotation model; nothing is
  invented. Fields that could not be computed are omitted rather than filled.

## Rows per reconstruction time

| time_Ma | n |
|---|---|
{n_table}

## Missing / not computed

{missing_fields if missing_fields else "None — all requested derived fields were computed."}

## Ranges (all times pooled)

- seafloor age (Myr): {_fmt_stats(age_stats)}
- conv_perp_cm_yr (cm/yr): {_fmt_stats(conv_stats)}
- v_T_perp (cm/yr): {_fmt_stats(vt_stats)}
- v_OP_perp (cm/yr): {_fmt_stats(vop_stats)}
"""
    path.write_text(text)


def main(argv=None):
    args = parse_args(argv)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    times = times_ma(args.tmin, args.tmax, args.dt)
    out_csv = Path(args.out_csv)
    print(
        f"Cenozoic extraction: {int(args.tmin)}–{int(args.tmax)} Ma, "
        f"dt={int(args.dt)} Myr, n_times={len(times)}"
    )
    print("times:", list(map(int, times)))

    print("Loading Müller2019 via PlateModelManager ...")
    model = PlateModelManager().get_model("Muller2019", data_dir=str(DATA_DIR))
    recon = PlateReconstruction(
        model.get_rotation_model(),
        model.get_topologies(),
        model.get_static_polygons(),
    )

    accepted, rejected_by_sig = k0._enable_extras(
        recon.tessellate_subduction_zones, k0.REQUESTED_EXTRAS
    )
    if rejected_by_sig:
        print("Extras not in installed signature:", rejected_by_sig)

    print("Ensuring Müller2019 age grids for requested times ...")
    age_paths = ensure_age_grids(model, times)

    tess_rad = np.radians(k0.TESSELLATION_DEG)
    frames = []
    extras_used = list(accepted)
    extra_err = None
    missing_union = []
    failed_times = []
    n_per_time = {}

    t_run = time_mod.time()
    for i, t in enumerate(times):
        t0 = time_mod.time()
        age_path = age_paths[int(t)]
        print(
            f"[{i + 1}/{len(times)}] tessellating {int(t)} Ma  "
            f"(age grid: {Path(age_path).name}) ...",
            flush=True,
        )
        try:
            df_t, extras_used, extra_err, missing = extract_one_time(
                recon, t, tess_rad, accepted, age_path
            )
        except Exception as err:
            print(f"  ERROR at {int(t)} Ma: {err}")
            failed_times.append(int(t))
            continue
        for m in missing:
            if m not in missing_union:
                missing_union.append(m)
        n_per_time[int(t)] = int(len(df_t))
        age_med = (
            float(np.nanmedian(df_t["seafloor_age_Ma"])) if len(df_t) else np.nan
        )
        print(
            f"  n={len(df_t)}  median seafloor age={age_med:.2f} Ma  "
            f"({time_mod.time() - t0:.2f}s)",
            flush=True,
        )
        if extra_err is not None and i == 0:
            print("  Last extra-kwarg error (recovered):", extra_err)
        if len(df_t):
            frames.append(df_t)

    if not frames:
        raise RuntimeError("No trench samples were extracted at any time.")

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"Wrote {len(out)} rows -> {out_csv}")
    print("CSV columns:", list(out.columns))

    age_stats = k0._stats(out["seafloor_age_Ma"])
    conv_stats = k0._stats(out["conv_perp_cm_yr"])
    vt_stats = k0._stats(out["v_T_perp"]) if "v_T_perp" in out.columns else None
    vop_stats = k0._stats(out["v_OP_perp"]) if "v_OP_perp" in out.columns else None

    # Document observed collinearity; do not interpret it as rollback.
    collinearity_note = ""
    if vt_stats is not None and vop_stats is not None:
        pair = out[["v_T_perp", "v_OP_perp"]].dropna()
        if len(pair) > 10:
            residual = pair["v_T_perp"] + pair["v_OP_perp"]
            r = float(np.corrcoef(pair["v_T_perp"], pair["v_OP_perp"])[0, 1])
            collinearity_note = (
                f"Observed over this run: corr(v_T_perp, v_OP_perp) = {r:.4f}; "
                f"median(v_T_perp + v_OP_perp) = {float(np.median(residual)):.4f} cm/yr "
                f"(near zero is the expected glued-trench relation)."
            )

    write_readme(
        README_PATH,
        n_rows=len(out),
        n_times=int(out["time_Ma"].nunique()),
        tmin=int(args.tmin),
        tmax=int(args.tmax),
        dt=int(args.dt),
        times=times,
        columns=list(out.columns),
        extras_used=extras_used,
        extras_rejected=rejected_by_sig,
        n_per_time=n_per_time,
        age_stats=age_stats,
        conv_stats=conv_stats,
        vt_stats=vt_stats,
        vop_stats=vop_stats,
        missing_fields="; ".join(missing_union) if missing_union else "",
        age_grid_count=len(age_paths),
        csv_path=out_csv,
        fig_path=FIG_PATH,
        failed_times=failed_times,
        collinearity_note=collinearity_note,
    )
    print(f"Wrote {README_PATH}")

    print("Drawing median seafloor age vs time ...")
    make_age_vs_time_figure(out, FIG_PATH)
    print(f"Wrote {FIG_PATH}")

    def _p(label, stats):
        if stats is None:
            print(f"  {label}: no finite values")
            return
        vmin, vmed, vmax, nvalid = stats
        print(
            f"  {label}: min={vmin:.4g}  median={vmed:.4g}  max={vmax:.4g}  n={nvalid}"
        )

    print("\n=== pooled ranges ===")
    _p("seafloor_age_Ma", age_stats)
    _p("conv_perp_cm_yr", conv_stats)
    _p("v_T_perp (cm/yr)", vt_stats)
    _p("v_OP_perp (cm/yr)", vop_stats)
    if collinearity_note:
        print(" ", collinearity_note)
    print(f"failed times: {failed_times or 'none'}")
    print(f"elapsed {time_mod.time() - t_run:.1f}s")
    print("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
