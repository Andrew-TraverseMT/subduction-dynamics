# 0 Ma trench kinematics — Müller et al. (2019)

Extracted by `scripts/extract_trench_kinematics_0Ma.py`.
This product is **trench kinematics only**. It is not joined to SubMap or UPS.

- **n rows:** 1238
- **reconstruction time:** 0 Ma
- **tessellation:** 0.5 deg (~55.6 km)
- **velocity delta time:** 1.0 Myr
- **Earth radius:** 6371.009 km (`pygplates.Earth.mean_radius_in_kms`)
- **GPlately extras accepted:** output_distance_to_nearest_edge_of_trench, output_convergence_velocity_components, output_trench_absolute_velocity_components, output_subducting_absolute_velocity, output_subducting_absolute_velocity_components
- **GPlately extras rejected:** (none)

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
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \
    /workspace/subduction-test/scripts/extract_trench_kinematics_0Ma.py
```

Equivalent micromamba one-shot:

```bash
export MAMBA_ROOT_PREFIX=/workspace/subduction-test/.micromamba
/workspace/subduction-test/.micromamba/bin/micromamba run -n pygplates \
    --root-prefix "$MAMBA_ROOT_PREFIX" python \
    /workspace/subduction-test/scripts/extract_trench_kinematics_0Ma.py
```

Outputs:

- `/workspace/subduction-test/results/trench_kinematics_0Ma_Muller2019.csv`
- `/workspace/subduction-test/results/fig_trench_vTperp_0Ma.png`
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
- `upper_plate_nature`

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

None — all requested derived fields were computed.

Kinematic values are taken from GPlately / the rotation model; nothing is
invented. Fields that could not be computed are omitted rather than filled.

## Velocity / age ranges (this run)

- seafloor age (Myr): min=0.240, median=49.635, max=338.810 (n_valid=779)
- v_T_perp (cm/yr): min=-4.669, median=0.541, max=8.621 (n_valid=1238)
- v_OP_perp (cm/yr): min=-8.656, median=-0.557, max=4.790 (n_valid=1238)
- conv_perp_cm_yr (cm/yr): min=-8.123, median=3.968, max=14.932 (n_valid=1238)
