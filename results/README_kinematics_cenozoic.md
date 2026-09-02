# Cenozoic trench kinematics — Müller et al. (2019)

Extracted by `scripts/extract_trench_kinematics_cenozoic.py`.
Reuses helpers from `scripts/extract_trench_kinematics_0Ma.py`.

This product is **trench kinematics only**. It is **not** joined to SubMap or
UPS. SubMap is present-day only; the 0 Ma H1 vs H2 comparison lives in
`results/submap_m56_h1h2.md` and must not be repeated from Müller2019
`v_T` vs `v_OP`.

- **n rows:** 76529
- **n times:** 61
- **time range:** 0–60 Ma
- **time step:** 1 Myr
- **times (Ma):** 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60
- **tessellation:** 0.5 deg (~55.6 km)
- **velocity delta time:** 1.0 Myr
- **Earth radius:** 6371.009 km (`pygplates.Earth.mean_radius_in_kms`)
- **age grids on disk for this run:** 61
- **GPlately extras accepted:** output_distance_to_nearest_edge_of_trench, output_convergence_velocity_components, output_trench_absolute_velocity_components, output_subducting_absolute_velocity, output_subducting_absolute_velocity_components
- **GPlately extras rejected:** (none)
- **failed times:** none

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
this model. Observed over this run: corr(v_T_perp, v_OP_perp) = -0.9999; median(v_T_perp + v_OP_perp) = -0.0000 cm/yr (near zero is the expected glued-trench relation).

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
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \
    /workspace/subduction-test/scripts/extract_trench_kinematics_cenozoic.py
```

Equivalent micromamba one-shot:

```bash
export MAMBA_ROOT_PREFIX=/workspace/subduction-test/.micromamba
/workspace/subduction-test/.micromamba/bin/micromamba run -n pygplates \
    --root-prefix "$MAMBA_ROOT_PREFIX" python \
    /workspace/subduction-test/scripts/extract_trench_kinematics_cenozoic.py
```

Coarser time step (downloads 13 age grids instead of 61):

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \
    /workspace/subduction-test/scripts/extract_trench_kinematics_cenozoic.py --dt 5
```

Already-downloaded age grids are reused; missing grids are fetched in batches.

Outputs:

- `/workspace/subduction-test/results/trench_kinematics_0-60Ma_Muller2019.csv`
- `/workspace/subduction-test/results/fig_age_at_trench_vs_time.png`
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

- `time_Ma`
- `lon`
- `lat`
- `subducting_plate_id`
- `trench_plate_id`
- `trench_normal_azimuth_deg`
- `arc_segment_length_km`
- `conv_mag_cm_yr`
- `conv_obliquity_deg`
- `conv_perp_cm_yr`
- `conv_parallel_cm_yr`
- `v_T_perp`
- `v_T_parallel`
- `v_OP_perp`
- `v_OP_parallel`
- `seafloor_age_Ma`
- `dist_nearest_trench_edge_km`
- `zone_id`
- `trench_zone_width_km`

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
| 0 | 1238 |
| 1 | 1280 |
| 2 | 1281 |
| 3 | 1297 |
| 4 | 1336 |
| 5 | 1297 |
| 6 | 1310 |
| 7 | 1287 |
| 8 | 1282 |
| 9 | 1297 |
| 10 | 1271 |
| 11 | 1268 |
| 12 | 1260 |
| 13 | 1275 |
| 14 | 1271 |
| 15 | 1298 |
| 16 | 1300 |
| 17 | 1329 |
| 18 | 1350 |
| 19 | 1347 |
| 20 | 1342 |
| 21 | 1306 |
| 22 | 1300 |
| 23 | 1291 |
| 24 | 1234 |
| 25 | 1240 |
| 26 | 1176 |
| 27 | 1184 |
| 28 | 1189 |
| 29 | 1209 |
| 30 | 1249 |
| 31 | 1261 |
| 32 | 1235 |
| 33 | 1225 |
| 34 | 1215 |
| 35 | 1232 |
| 36 | 1178 |
| 37 | 1150 |
| 38 | 1144 |
| 39 | 1140 |
| 40 | 1176 |
| 41 | 1120 |
| 42 | 1116 |
| 43 | 1146 |
| 44 | 1140 |
| 45 | 1147 |
| 46 | 1214 |
| 47 | 1198 |
| 48 | 1285 |
| 49 | 1283 |
| 50 | 1295 |
| 51 | 1225 |
| 52 | 1327 |
| 53 | 1330 |
| 54 | 1335 |
| 55 | 1338 |
| 56 | 1289 |
| 57 | 1308 |
| 58 | 1295 |
| 59 | 1292 |
| 60 | 1296 |

## Missing / not computed

None — all requested derived fields were computed.

## Ranges (all times pooled)

- seafloor age (Myr): min=0.010, median=47.614, max=338.810 (n_valid=57186)
- conv_perp_cm_yr (cm/yr): min=-12.945, median=4.093, max=25.196 (n_valid=76529)
- v_T_perp (cm/yr): min=-13.626, median=1.015, max=21.897 (n_valid=76529)
- v_OP_perp (cm/yr): min=-21.897, median=-1.014, max=13.626 (n_valid=76529)
