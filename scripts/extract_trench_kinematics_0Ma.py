#!/usr/bin/env python3
"""Extract 0 Ma trench kinematics from the Müller et al. (2019) plate model.

Reproducible, local-data-only workflow for the subduction-dynamics project.
Does not join SubMap / UPS.

Run:
    /workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
        /workspace/subduction-test/scripts/extract_trench_kinematics_0Ma.py
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pygplates
from shapely.geometry import Point
from shapely.ops import unary_union
from shapely.prepared import prep

import gplately
from gplately import PlateModelManager, PlateReconstruction
from gplately.grids import Raster

# ---------------------------------------------------------------------------
# Paths (no hidden locations)
# ---------------------------------------------------------------------------
DATA_DIR = Path("/workspace/subduction-test/data")
RESULTS_DIR = Path("/workspace/subduction-test/results")
CSV_PATH = RESULTS_DIR / "trench_kinematics_0Ma_Muller2019.csv"
README_PATH = RESULTS_DIR / "README_kinematics_0Ma.md"
FIG_PATH = RESULTS_DIR / "fig_trench_vTperp_0Ma.png"

TIME_MA = 0.0
# ~0.5 deg ≈ 55.6 km along a great circle (within the requested 50–100 km).
TESSELLATION_DEG = 0.5
VELOCITY_DELTA_TIME_MYR = 1.0
EARTH_RADIUS_KM = float(pygplates.Earth.mean_radius_in_kms)
DEG_TO_KM = np.pi / 180.0 * EARTH_RADIUS_KM
# Offset used to classify upper-plate nature / rescue seafloor age.
OFFSET_INTO_PLATE_KM = 50.0

# Extra tessellate kwargs requested by the project. Order MUST match GPlately's
# append order so that GeoDataFrame column names stay aligned if we fall back
# to the ndarray interface.
REQUESTED_EXTRAS = [
    "output_distance_to_nearest_edge_of_trench",
    "output_convergence_velocity_components",
    "output_trench_absolute_velocity_components",
    "output_subducting_absolute_velocity",
    "output_subducting_absolute_velocity_components",
]


def _destination_lonlat(lon, lat, azimuth_deg, distance_km):
    """Spherical destination: azimuth clockwise from North, distance in km.

    Returns (lon, lat) in degrees, lon wrapped to [-180, 180].
    """
    lat1 = np.radians(lat)
    lon1 = np.radians(lon)
    az = np.radians(azimuth_deg)
    d = np.asarray(distance_km, dtype=float) / EARTH_RADIUS_KM
    sin_lat1 = np.sin(lat1)
    cos_lat1 = np.cos(lat1)
    sin_d = np.sin(d)
    cos_d = np.cos(d)
    lat2 = np.arcsin(sin_lat1 * cos_d + cos_lat1 * sin_d * np.cos(az))
    lon2 = lon1 + np.arctan2(
        np.sin(az) * sin_d * cos_lat1, cos_d - sin_lat1 * np.sin(lat2)
    )
    lon2_deg = (np.degrees(lon2) + 180.0) % 360.0 - 180.0
    return lon2_deg, np.degrees(lat2)


def _great_circle_km(lon1, lat1, lon2, lat2):
    """Vectorised great-circle distance in km (haversine)."""
    p1 = np.radians(lat1)
    p2 = np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlmb = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2.0) ** 2
    return 2.0 * EARTH_RADIUS_KM * np.arcsin(np.clip(np.sqrt(a), 0.0, 1.0))


def _enable_extras(func, requested):
    """Keep only extra kwargs that exist on the installed tessellate signature."""
    params = inspect.signature(func).parameters
    accepted, rejected = [], []
    for name in requested:
        if name in params:
            accepted.append(name)
        else:
            rejected.append(name)
    return accepted, rejected


def tessellate_with_fallback(recon, time, tess_rad, extra_names):
    """Call tessellate_subduction_zones; drop extras one-by-one if a kwarg errors."""
    extras = list(extra_names)
    last_err = None
    while True:
        kwargs = {name: True for name in extras}
        try:
            gdf = recon.tessellate_subduction_zones(
                time,
                tessellation_threshold_radians=tess_rad,
                ignore_warnings=True,
                return_geodataframe=True,
                **kwargs,
            )
            return gdf, extras, last_err
        except TypeError as err:
            last_err = err
            if not extras:
                raise
            dropped = extras.pop()
            print(f"WARNING: tessellate rejected extras, dropping '{dropped}': {err}")


def plate_abs_velocity_east_north(rotation_model, lons, lats, plate_ids, time, delta_time):
    """Absolute velocity of each sample's plate ID, in cm/yr, east and north.

    Uses the rotation model directly (not topology-located velocities), so the
    trench plate ID is honoured even if the sample sits slightly on the
    subducting side of the boundary.
    Stage rotation is [time + delta_time -> time], matching GPlately's default
    VelocityDeltaTimeType.t_plus_delta_t_to_t.
    """
    n = len(lons)
    v_east = np.full(n, np.nan, dtype=float)
    v_north = np.full(n, np.nan, dtype=float)
    pids = np.asarray(plate_ids, dtype=int)
    for pid in np.unique(pids):
        mask = pids == pid
        idxs = np.flatnonzero(mask)
        points = [
            pygplates.PointOnSphere(float(lats[i]), float(lons[i])) for i in idxs
        ]
        stage = rotation_model.get_rotation(time, int(pid), time + delta_time)
        if stage is None:
            continue
        vels = pygplates.calculate_velocities(
            points,
            stage,
            delta_time,
            velocity_units=pygplates.VelocityUnits.cms_per_yr,
        )
        ned = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(points, vels)
        for local_i, vec in enumerate(ned):
            if vec is None:
                continue
            # NED: x=north, y=east, z=down
            v_north[idxs[local_i]] = vec.get_x()
            v_east[idxs[local_i]] = vec.get_y()
    return v_east, v_north


def trench_normal_components(v_east, v_north, azimuth_deg):
    """Project east/north velocity onto trench-normal and trench-parallel.

    Trench normal azimuth is clockwise from North and points toward the
    overriding plate (GPlately convention).

    Parallel is 90 deg clockwise from the trench normal (same as GPlately
    extra-output components).
    """
    az = np.radians(np.asarray(azimuth_deg, dtype=float))
    n_north = np.cos(az)
    n_east = np.sin(az)
    # 90 deg clockwise: azimuth + 90
    p_north = -np.sin(az)
    p_east = np.cos(az)
    v_perp = v_north * n_north + v_east * n_east
    v_par = v_north * p_north + v_east * p_east
    return v_perp, v_par


def connected_zone_width_km(lons, lats, sub_ids, trench_ids, segment_len_km, tess_km):
    """Along-trench length of each connected subduction polyline.

    A zone is a spatially connected set of tessellated samples that share
    the same subducting plate ID. Neighbours are linked when their
    great-circle separation is <= 2.5 * tessellation spacing (~140 km at
    0.5 deg). Matching trench_plate_id is NOT required: Müller2019
    deforming networks assign many distinct trench IDs along one continuous
    subduction polyline (Hellenic, Andes, etc.), which would otherwise
    fragment a single slab into tens of tiny zones.

    Disconnected trenches of the same subducting plate (e.g. Pacific under
    Tonga vs Pacific under Japan) stay separate because they are not
    spatially linked at this threshold.

    trench_ids is accepted for API stability but is not used in grouping.

    Returns (zone_id, zone_width_km) aligned with the input samples.
    """
    n = len(lons)
    if n == 0:
        return np.array([], dtype=int), np.array([], dtype=float)

    xyz = np.column_stack(
        (
            np.cos(np.radians(lats)) * np.cos(np.radians(lons)),
            np.cos(np.radians(lats)) * np.sin(np.radians(lons)),
            np.sin(np.radians(lats)),
        )
    )
    parent = np.arange(n, dtype=int)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    neighbor_km = 2.5 * tess_km
    chord = 2.0 * np.sin((neighbor_km / EARTH_RADIUS_KM) / 2.0)
    chord2 = chord * chord

    from collections import defaultdict

    buckets = defaultdict(list)
    for i, sid in enumerate(np.asarray(sub_ids, dtype=int)):
        buckets[int(sid)].append(i)

    for idxs in buckets.values():
        if len(idxs) < 2:
            continue
        pts = xyz[idxs]
        for a in range(len(idxs)):
            d2 = np.sum((pts[a + 1 :] - pts[a]) ** 2, axis=1)
            for h in np.flatnonzero(d2 <= chord2):
                union(idxs[a], idxs[a + 1 + int(h)])

    roots = np.array([find(i) for i in range(n)])
    _, zone_id = np.unique(roots, return_inverse=True)
    widths = np.empty(n, dtype=float)
    for z in np.unique(zone_id):
        m = zone_id == z
        widths[m] = float(np.nansum(segment_len_km[m]))
    return zone_id, widths


def points_in_polygons(lons, lats, geom_union):
    """Boolean mask: sample in the (prepared) polygonal union."""
    prepared = prep(geom_union)
    out = np.zeros(len(lons), dtype=bool)
    for i, (x, y) in enumerate(zip(lons, lats)):
        out[i] = prepared.contains(Point(float(x), float(y)))
    return out


def write_readme(
    path,
    n_rows,
    columns,
    extras_used,
    extras_rejected,
    tess_deg,
    age_stats,
    vt_stats,
    vop_stats,
    conv_stats,
    missing_fields,
):
    def _fmt(stats):
        if stats is None:
            return "not computed"
        vmin, vmed, vmax, nvalid = stats
        return (
            f"min={vmin:.3f}, median={vmed:.3f}, max={vmax:.3f} "
            f"(n_valid={nvalid})"
        )

    text = f"""# 0 Ma trench kinematics — Müller et al. (2019)

Extracted by `scripts/extract_trench_kinematics_0Ma.py`.
This product is **trench kinematics only**. It is not joined to SubMap or UPS.

- **n rows:** {n_rows}
- **reconstruction time:** 0 Ma
- **tessellation:** {tess_deg} deg (~{tess_deg * DEG_TO_KM:.1f} km)
- **velocity delta time:** {VELOCITY_DELTA_TIME_MYR} Myr
- **Earth radius:** {EARTH_RADIUS_KM} km (`pygplates.Earth.mean_radius_in_kms`)
- **GPlately extras accepted:** {', '.join(extras_used) if extras_used else '(none)'}
- **GPlately extras rejected:** {', '.join(extras_rejected) if extras_rejected else '(none)'}

## Model citation

Müller, R. D., Zahirovic, S., Williams, S. E., Cannon, J., Seton, M., Bower, D. J.,
Tetley, M. G., Heine, C., Le Breton, E., Liu, S., Russell, S. H. J., Yang, T.,
Leonard, J., & Gurnis, M. (2019). A global plate model including lithospheric
deformation along major rifts and orogens since the Triassic. *Tectonics*,
38(6), 1884–1907. https://doi.org/10.1029/2018TC005462

Digital plate-model package (GPlately / GPlates):
https://doi.org/10.5281/zenodo.10525286
(version recorded in local metadata: `10.5281/zenodo.11601026`).

Local data directory: `/workspace/subduction-test/data/Muller2019`
(PlateModelManager also resolves the lowercase symlink `data/muller2019`).

## How to rerun

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
    /workspace/subduction-test/scripts/extract_trench_kinematics_0Ma.py
```

Equivalent micromamba one-shot:

```bash
export MAMBA_ROOT_PREFIX=/workspace/subduction-test/.micromamba
/workspace/subduction-test/.micromamba/bin/micromamba run -n pygplates \\
    --root-prefix "$MAMBA_ROOT_PREFIX" python \\
    /workspace/subduction-test/scripts/extract_trench_kinematics_0Ma.py
```

Outputs:

- `{CSV_PATH}`
- `{FIG_PATH}`
- this README

## Sign conventions (read this before using the CSV)

GPlately's trench **normal always points toward the overriding plate**.
Obliquity is the signed angle between that normal and the velocity vector,
range (−180, 180) deg: positive is clockwise from the trench normal when
viewed from above. Parallel components are positive 90 deg clockwise from
the trench normal.

GPlately also encodes a sign in some *magnitudes*:

- convergence magnitude is **negative if the plates are diverging**
  (|obliquity| > 90).
- trench absolute-velocity magnitude is **negative if the trench is moving
  toward the overriding plate** (|obliquity| < 90). That is landward trench
  *advance*, not rollback.

Orthogonal extra outputs use `cos(obliquity) * |magnitude|`, so:

- **convergence orthogonal > 0** = subducting plate moving toward the
  overriding plate (convergence).
- **trench-absolute orthogonal > 0** = trench moving toward the overriding
  plate (landward advance).

This CSV **flips the trench-absolute orthogonal** so that rollback is
positive, and computes overriding-plate motion from the rotation model.

| field | positive means |
|---|---|
| `conv_perp_cm_yr` | trench-normal convergence: subducting plate approaching the overriding plate |
| `conv_parallel_cm_yr` | trench-parallel relative motion, 90 deg clockwise from trench normal |
| `v_T_perp` | trench-normal **trench** absolute velocity: **oceanward trench retreat** (hinge rolling back) |
| `v_T_parallel` | trench-parallel trench absolute velocity, 90 deg clockwise from trench normal |
| `v_OP_perp` | trench-normal **overriding-plate** absolute velocity: upper plate moving **away from the trench** |
| `v_OP_parallel` | trench-parallel overriding-plate absolute velocity, 90 deg clockwise from trench normal |

Derivation:

```
v_T_perp      = - trench_absolute_velocity_orthogonal
               # GPlately orth > 0 is landward advance; we want retreat > 0
v_OP_perp     = v_OP · n_hat
               # n_hat = trench normal toward OP
               # so v_OP_perp > 0 is OP moving away from the trench
```

`v_OP_*` is the absolute velocity of `trench_plate_id` (the GPlately trench
plate ID; this can differ from a purely topological overriding plate in
deforming models) evaluated at the sample with the rotation model.

## Column dictionary

{chr(10).join(f'- `{c}`' for c in columns)}

Detailed meanings:

| column | units | meaning |
|---|---|---|
| lon, lat | deg | tessellated trench sample |
| subducting_plate_id | int | GPlately subducting plate ID |
| trench_plate_id | int | GPlately trench plate ID (used as overriding plate) |
| trench_normal_azimuth_deg | deg | clockwise from North, toward overriding plate |
| arc_segment_length_km | km | length of the tessellated arc segment this sample sits on |
| conv_mag_cm_yr | cm/yr | GPlately signed convergence magnitude (negative if diverging) |
| conv_obliquity_deg | deg | angle trench-normal → convergence velocity |
| conv_perp_cm_yr | cm/yr | trench-normal convergence; + = converging |
| conv_parallel_cm_yr | cm/yr | trench-parallel convergence component |
| v_T_perp | cm/yr | trench-normal trench abs. vel.; + = oceanward retreat |
| v_T_parallel | cm/yr | trench-parallel trench abs. vel. |
| v_OP_perp | cm/yr | trench-normal OP abs. vel.; + = OP moving away from trench |
| v_OP_parallel | cm/yr | trench-parallel OP abs. vel. |
| seafloor_age_Ma | Myr | 0 Ma age-grid sample; NaN on continents / missing |
| dist_nearest_trench_edge_km | km | distance along trench to nearest slab/trench edge |
| zone_id | int | spatially connected polyline id (same subducting plate ID) |
| trench_zone_width_km | km | along-trench length of the spatially connected polyline sharing the same subducting plate ID |
| upper_plate_nature | str | `continental` or `oceanic` |

## Sampling notes

- **Seafloor age** is sampled from the 0 Ma Müller age grid, preferring a
  point 50 km oceanward of the trench (toward the subducting plate) so
  Andean-type trenches that sit on the continental side of the raster still
  receive oceanic lithosphere age. The trench node is the fallback. NaN
  remains where the raster has no seafloor (continents, some collision
  zones). The raster itself contains values up to ~339 Ma in the eastern
  Mediterranean; those are source-grid values, not filled by this script.
- **`trench_zone_width_km`** is the along-trench length of each spatially
  connected set of samples that share a subducting plate ID (neighbour
  gap ≤ 2.5 × tessellation). Disconnected trenches of the same plate stay
  separate (Tonga vs Japan). Continuous chains stay together: the
  Aleutian–Kamchatka–Kuril–Japan–Izu–Bonin–Mariana Pacific polyline is one
  zone (~10,000 km); the Nazca–Andes trench is one zone (~6,900 km).
- In this rigid/deforming-network model the trench line is typically
  attached to `trench_plate_id`, so `v_T_perp ≈ −v_OP_perp`. Both columns
  are retained because they are defined independently.

## Missing / not computed

{missing_fields if missing_fields else "None — all requested derived fields were computed."}

Kinematic values are taken from GPlately / the rotation model; nothing is
invented. Fields that could not be computed are omitted rather than filled.

## Velocity / age ranges (this run)

- seafloor age (Myr): {_fmt(age_stats)}
- v_T_perp (cm/yr): {_fmt(vt_stats)}
- v_OP_perp (cm/yr): {_fmt(vop_stats)}
- conv_perp_cm_yr (cm/yr): {_fmt(conv_stats)}
"""
    path.write_text(text)


def _stats(arr):
    a = np.asarray(arr, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return None
    return float(np.min(a)), float(np.median(a)), float(np.max(a)), int(a.size)


def make_sanity_map(df, path):
    """Global scatter of trench samples coloured by v_T_perp."""
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        proj = ccrs.Robinson()
        data_crs = ccrs.PlateCarree()
        fig, ax = plt.subplots(figsize=(12.5, 6.5), subplot_kw={"projection": proj})
        ax.set_global()
        try:
            ax.add_feature(cfeature.LAND, facecolor="0.88", edgecolor="none", zorder=0)
            ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor="0.35", zorder=1)
        except Exception as err:
            print(f"WARNING: cartopy NaturalEarth features unavailable ({err}); plotting without coasts.")
        ax.gridlines(draw_labels=False, linewidth=0.3, color="0.7", linestyle=":")
        scatter_transform = data_crs
    except Exception as err:
        print(f"WARNING: cartopy unavailable ({err}); falling back to a plain lon/lat plot.")
        fig, ax = plt.subplots(figsize=(12.5, 6.5))
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_aspect("equal", adjustable="box")
        scatter_transform = None

    v = df["v_T_perp"].to_numpy()
    vmax = np.nanpercentile(np.abs(v), 95)
    vmax = max(float(vmax), 1.0)
    kw = dict(
        c=v,
        s=12,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        linewidths=0,
        zorder=3,
    )
    if scatter_transform is not None:
        sc = ax.scatter(df["lon"], df["lat"], transform=scatter_transform, **kw)
    else:
        sc = ax.scatter(df["lon"], df["lat"], **kw)
    cbar = fig.colorbar(sc, ax=ax, shrink=0.72, pad=0.02)
    cbar.set_label(r"$v_{T\perp}$ (cm/yr)  + oceanward retreat")
    ax.set_title(
        "Müller et al. (2019) trench samples at 0 Ma — coloured by $v_{T\\perp}$"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading Müller2019 via PlateModelManager ...")
    model = PlateModelManager().get_model("Muller2019", data_dir=str(DATA_DIR))
    recon = PlateReconstruction(
        model.get_rotation_model(),
        model.get_topologies(),
        model.get_static_polygons(),
    )

    accepted, rejected_by_sig = _enable_extras(
        recon.tessellate_subduction_zones, REQUESTED_EXTRAS
    )
    if rejected_by_sig:
        print("Extras not in installed signature:", rejected_by_sig)

    tess_rad = np.radians(TESSELLATION_DEG)
    print(
        f"Tessellating subduction zones at {TIME_MA} Ma, "
        f"threshold={TESSELLATION_DEG} deg ({tess_rad:.5f} rad) ..."
    )
    gdf, extras_used, extra_err = tessellate_with_fallback(
        recon, TIME_MA, tess_rad, accepted
    )
    if extra_err is not None:
        print("Last extra-kwarg error (recovered):", extra_err)

    print("\n=== Actual tessellate column layout ===")
    print("n samples:", len(gdf))
    print("columns:")
    for i, col in enumerate(gdf.columns):
        print(f"  [{i:02d}] {col}")
    print("dtypes:\n", gdf.dtypes)
    print("head:\n", gdf.head())
    print("======================================\n")

    lon = gdf.geometry.x.to_numpy() if "geometry" in gdf.columns else gdf["lon"].to_numpy()
    lat = gdf.geometry.y.to_numpy() if "geometry" in gdf.columns else gdf["lat"].to_numpy()
    # GeoDataFrame from GPlately stores lon/lat only in geometry, not as columns.

    def col(*names):
        for name in names:
            if name in gdf.columns:
                return gdf[name].to_numpy()
        return None

    sub_id = col("subducting plate ID")
    trench_id = col("trench plate ID")
    conv_mag = col("convergence velocity (cm/yr)")
    conv_obliq = col("convergence obliquity angle (degrees)")
    trench_az = col("trench normal angle (degrees)")
    arc_deg = col("length (degrees)")
    conv_orth = col("convergence velocity orthogonal component (cm/yr)")
    conv_par = col("convergence velocity parallel component (cm/yr)")
    trench_orth = col("trench absolute velocity orthogonal component (cm/yr)")
    trench_par = col("trench absolute velocity parallel component (cm/yr)")
    dist_edge_deg = col("distance to nearest trench edge (degrees)")

    missing = []

    # Convergence components. Prefer GPlately extras; otherwise reconstruct
    # from signed magnitude and obliquity using the same formula GPlately uses:
    #   orth = cos(obliq) * |mag|,  par = sin(obliq) * |mag|
    if conv_orth is None or conv_par is None:
        if conv_mag is not None and conv_obliq is not None:
            rad = np.radians(conv_obliq)
            conv_orth = np.cos(rad) * np.abs(conv_mag)
            conv_par = np.sin(rad) * np.abs(conv_mag)
            print("Computed conv components from magnitude/obliquity (extra kwarg missing).")
        else:
            missing.append("convergence trench-normal / trench-parallel")

    # v_T_perp: flip GPlately trench-absolute orthogonal so retreat is positive.
    # GPlately: trench_orth > 0  => trench moving toward OP (landward ADVANCE)
    # Ours:     v_T_perp > 0     => trench moving oceanward (RETREAT / rollback)
    if trench_orth is None:
        trench_abs_mag = col("trench velocity (cm/yr)")
        trench_abs_obliq = col("trench obliquity angle (degrees)")
        if trench_abs_mag is not None and trench_abs_obliq is not None:
            rad = np.radians(trench_abs_obliq)
            trench_orth = np.cos(rad) * np.abs(trench_abs_mag)
            trench_par = np.sin(rad) * np.abs(trench_abs_mag)
            print("Computed trench-abs components from magnitude/obliquity (extra kwarg missing).")
        else:
            missing.append("v_T_perp / v_T_parallel")
            v_T_perp = None
            v_T_par = None
    if trench_orth is not None:
        v_T_perp = -np.asarray(trench_orth, dtype=float)
        v_T_par = np.asarray(trench_par, dtype=float)

    # Overriding-plate velocity from the rotation model at trench_plate_id.
    print("Computing overriding-plate velocities from the rotation model ...")
    v_east, v_north = plate_abs_velocity_east_north(
        recon.rotation_model,
        lon,
        lat,
        trench_id,
        TIME_MA,
        VELOCITY_DELTA_TIME_MYR,
    )
    v_OP_perp, v_OP_par = trench_normal_components(v_east, v_north, trench_az)

    # Seafloor age at the trench. Prefer a sample 50 km toward the
    # subducting plate (oceanward, opposite the trench normal) so Andean-
    # type trenches that sit on the continental side of the age grid still
    # get oceanic lithosphere age. Fall back to the trench node. Leave NaN
    # where the 0 Ma age grid is missing (continents / no seafloor).
    # The raster itself contains values up to ~339 Ma (eastern Mediterranean);
    # those are source-grid values, not filled in by this script.
    print("Sampling 0 Ma seafloor age grid ...")
    age_path = model.get_age_grid(0)
    print("age grid:", age_path)
    age_raster = Raster(str(age_path))
    age_trench = np.asarray(age_raster.interpolate(lon, lat, method="linear"), dtype=float)
    lon_sp, lat_sp = _destination_lonlat(lon, lat, trench_az + 180.0, OFFSET_INTO_PLATE_KM)
    age_sp = np.asarray(age_raster.interpolate(lon_sp, lat_sp, method="linear"), dtype=float)
    seafloor_age = np.where(np.isfinite(age_sp), age_sp, age_trench)

    # Distance to nearest trench/slab edge.
    if dist_edge_deg is None:
        missing.append("distance to nearest trench/slab edge")
        dist_edge_km = None
    else:
        dist_edge_km = np.asarray(dist_edge_deg, dtype=float) * DEG_TO_KM

    arc_km = np.asarray(arc_deg, dtype=float) * DEG_TO_KM
    tess_km = TESSELLATION_DEG * DEG_TO_KM
    zone_id, zone_width_km = connected_zone_width_km(
        lon, lat, sub_id.astype(int), trench_id.astype(int), arc_km, tess_km
    )

    # Upper-plate nature: point ~50 km into the overriding plate vs continental polygons.
    print("Classifying upper-plate nature against continental polygons ...")
    cont_files = model.get_continental_polygons()
    if isinstance(cont_files, (list, tuple)):
        cont_path = cont_files[0]
    else:
        cont_path = cont_files
    print("continental polygons:", cont_path)
    cont_gdf = gpd.read_file(cont_path)
    cont_union = unary_union(cont_gdf.geometry)
    lon_op, lat_op = _destination_lonlat(lon, lat, trench_az, OFFSET_INTO_PLATE_KM)
    is_cont = points_in_polygons(lon_op, lat_op, cont_union)
    # If the offset point misses (coastal geometry), try the trench sample itself.
    is_cont_trench = points_in_polygons(lon, lat, cont_union)
    is_cont = np.logical_or(is_cont, is_cont_trench)
    upper_plate_nature = np.where(is_cont, "continental", "oceanic")

    out = pd.DataFrame(
        {
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
            "v_T_perp": v_T_perp,
            "v_T_parallel": v_T_par,
            "v_OP_perp": v_OP_perp,
            "v_OP_parallel": v_OP_par,
            "seafloor_age_Ma": seafloor_age,
            "dist_nearest_trench_edge_km": dist_edge_km
            if dist_edge_km is not None
            else np.full(len(lon), np.nan),
            "zone_id": zone_id,
            "trench_zone_width_km": zone_width_km,
            "upper_plate_nature": upper_plate_nature,
        }
    )
    if dist_edge_km is None:
        out = out.drop(columns=["dist_nearest_trench_edge_km"])

    out.to_csv(CSV_PATH, index=False, float_format="%.6f")
    print(f"Wrote {len(out)} rows -> {CSV_PATH}")
    print("CSV columns:", list(out.columns))

    age_stats = _stats(out["seafloor_age_Ma"])
    vt_stats = _stats(out["v_T_perp"])
    vop_stats = _stats(out["v_OP_perp"])
    conv_stats = _stats(out["conv_perp_cm_yr"])

    write_readme(
        README_PATH,
        n_rows=len(out),
        columns=list(out.columns),
        extras_used=extras_used,
        extras_rejected=rejected_by_sig,
        tess_deg=TESSELLATION_DEG,
        age_stats=age_stats,
        vt_stats=vt_stats,
        vop_stats=vop_stats,
        conv_stats=conv_stats,
        missing_fields="; ".join(missing) if missing else "",
    )
    print(f"Wrote {README_PATH}")

    print("Drawing sanity map ...")
    make_sanity_map(out, FIG_PATH)
    print(f"Wrote {FIG_PATH}")

    def _p(label, stats):
        if stats is None:
            print(f"  {label}: no finite values")
            return
        vmin, vmed, vmax, nvalid = stats
        print(
            f"  {label}: min={vmin:.4g}  median={vmed:.4g}  max={vmax:.4g}  n={nvalid}"
        )

    print("\n=== velocity / age ranges ===")
    _p("seafloor_age_Ma", age_stats)
    _p("v_T_perp (cm/yr)", vt_stats)
    _p("v_OP_perp (cm/yr)", vop_stats)
    _p("conv_perp_cm_yr", conv_stats)
    print(
        "upper_plate_nature counts:",
        out["upper_plate_nature"].value_counts().to_dict(),
    )
    print("n zones:", int(out["zone_id"].nunique()))
    print("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
