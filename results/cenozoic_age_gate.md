# Cenozoic age-gate test — Müller et al. (2019) back-arc spreading proxy

Built by `scripts/cenozoic_ups_proxy.py` from the existing trench catalog
`results/trench_kinematics_0-60Ma_Muller2019.csv`. Trench kinematics were
**not** re-extracted.

- **n trench samples:** 76529
- **n times:** 61 (0–60 Ma, 1 Myr)
- **n bab_on (point samples):** 3858 (5.0%)
- **n zone–time blocks:** 2077
- **failed times:** none
- **output CSV:** `/workspace/subduction-test/results/trench_kinematics_0-60Ma_with_ups_proxy.csv`

## Circularity (read this first)

This upper-plate strain proxy is taken from the **same** Müller et al. (2019) reconstruction that supplies the seafloor-age grids and the trench kinematics. Back-arc basins exist in the model because the topologies put spreading ridges there. The Cenozoic test is therefore a **Sdrolias & Müller (2006)-style kinematic association** (does the reconstruction open a back-arc where the subducting plate is old?), not a test against independent geological UPS. Present-day SubMap geodetic UPS remains the only observation that is independent of this plate model. That circularity is a feature of using reconstructed back-arc spreading as UPS, and it is the reason the 0 Ma SubMap table is reported as a sanity check rather than a confirmation.

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

samples in a deforming network: 26189 / 76529 (34.2%); times with ≥1 network hit: 61 / 61

Layer B is **model-internal**. It is kept separate from `bab_on`. Aegean /
Hellenic extension is a deforming network in Müller2019, not a MOR, so it
is invisible to layer A and only shows up in B.

## 0 Ma spot check

| region | n | n bab_on | frac | med dist (km) | med dilat (10⁻¹⁵ s⁻¹) |
|---|---:|---:|---:|---:|---:|
| Mariana | 40 | 39 | 0.97 | 225 | 0.000 |
| Tonga | 32 | 28 | 0.88 | 307 | 0.000 |
| Scotia | 17 | 17 | 1.00 | 227 | 0.000 |
| Hellenic | 34 | 0 | 0.00 | — | 0.790 |
| Andes | 66 | 0 | 0.00 | — | 0.060 |
| Cascadia | 22 | 0 | 0.00 | — | -0.143 |
| Andaman | 17 | 17 | 1.00 | 218 | 0.000 |
| MAT_Mexico | 17 | 0 | 0.00 | — | 0.001 |

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

At **0 Ma** the Sdrolias pattern is present in this reconstruction:
median age bab_on vs off is 97.5 vs 44.0 Ma, and
P(bab_on | age ≥ 55) = 0.25 vs P(bab_on | age < 55) = 0.02.
SubMap geodetic C has 0/36 `bab_on`; E has 12/39. That is agreement in
direction, not a winner claim (Hellenic is geodetic E but `bab_on=0`
because Müller2019 has no Aegean MOR).

Pooled over 0–60 Ma, **point samples** still associate older slabs with a
back-arc MOR (median 72.3 vs 46.4 Ma;
Mantel–Haenszel OR 2.49). That n is not independent — long
old-Pacific trenches (Mariana, Tonga) contribute hundreds of neighbouring
points.

**Zone–time blocks** (n=2077, of which 187 have a back-arc) are the
unit that does not pretend 76k rows are independent. There the 55 Ma gate
is **weak** (zone-level MH OR 1.21, p=0.2391).
Univariate Spearman at zone–time: trench-zone width ρ = 0.216
outranks distance-to-edge, which outranks age (ρ = 0.042,
p = 0.0717). Logistic age after time dummies is consistent with
null (n=1831 (n_on=177); age coef = 0.0016 per Ma (OR=1.0016); LRT vs time-only χ²=0.63, p=0.4257).

The 0 Ma age headline therefore **does not survive** Cenozoic
zone-blocking in this circular proxy. Width is the stronger zone-level
associate of reconstructed back-arc spreading — a Schellart-shaped
geometric pattern in the topologies, still not independent geology.

The 55 Ma gate also **flips** in some intervals (about 12–14 Ma and
51–60 Ma) where P(on | young) ≥ P(on | old). Steps in the `bab_on`
fraction (for example near 6, 12, and 15 Ma) likely mark topology
revisions, not smooth basin opening.

### 1. Median seafloor age, bab_on vs bab_off

Pooled point samples, finite ages: U=112862351, p=1.6e-193, n_on=3185, n_off=54001; median age on=72.3 Ma, off=46.4 Ma

Per-time medians are in the table below and in `fig_bab_vs_age.png`.

### 2. Sdrolias gate (55 Ma) and a threshold sweep

P(bab_on | age ≥ T) versus P(bab_on | age < T).

**Point-level sweep (pooled over time):**

| T (Ma) | n old | n young | P(on\|old) | P(on\|young) | Δ |
|---|---:|---:|---:|---:|---:|
| 30 | 39594 | 17592 | 0.070 | 0.023 | 0.047 |
| 40 | 33934 | 23252 | 0.076 | 0.026 | 0.050 |
| 45 | 30398 | 26788 | 0.081 | 0.027 | 0.053 |
| 50 | 27246 | 29940 | 0.083 | 0.031 | 0.052 |
| 55 | 24751 | 32435 | 0.083 | 0.035 | 0.048 |
| 60 | 22837 | 34349 | 0.083 | 0.037 | 0.046 |
| 70 | 18462 | 38724 | 0.091 | 0.039 | 0.052 |
| 80 | 13929 | 43257 | 0.100 | 0.041 | 0.059 |
| 100 | 7030 | 50156 | 0.107 | 0.048 | 0.059 |

**Zone–time sweep:**

| T (Ma) | n old | n young | P(on\|old) | P(on\|young) | Δ |
|---|---:|---:|---:|---:|---:|
| 30 | 1176 | 655 | 0.118 | 0.058 | 0.060 |
| 40 | 1021 | 810 | 0.127 | 0.058 | 0.069 |
| 45 | 941 | 890 | 0.120 | 0.072 | 0.048 |
| 50 | 836 | 995 | 0.104 | 0.090 | 0.014 |
| 55 | 771 | 1060 | 0.104 | 0.092 | 0.012 |
| 60 | 699 | 1132 | 0.094 | 0.098 | -0.004 |
| 70 | 572 | 1259 | 0.091 | 0.099 | -0.008 |
| 80 | 474 | 1357 | 0.101 | 0.095 | 0.006 |
| 100 | 269 | 1562 | 0.104 | 0.095 | 0.009 |

Mantel–Haenszel odds ratio for (age ≥ 55) vs bab, stratified by time:
point-level MH OR = 2.488 (χ²=643.01, p=7.4e-142, 61 time strata); zone-level MH OR = 1.205 (χ²=1.39, p=0.2391, 61 strata)

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
| seafloor age (Ma) | 57186 | 0.124 | 5.4e-195 | 1.0082 | 1.1e-104 |
| trench-zone width (km) | 76529 | 0.047 | 4.5e-39 | 1.0000 | 4.5e-50 |
| distance to trench edge (km) | 76529 | 0.021 | 3.2e-09 | 1.0002 | 3.2e-07 |

**Zone–time blocks (preferred n):**

| predictor | n | Spearman ρ | p | logistic OR | LRT p |
|---|---:|---:|---:|---:|---:|
| seafloor age (Ma) | 1831 | 0.042 | 0.0717 | 1.0014 | 0.4661 |
| trench-zone width (km) | 2077 | 0.216 | 2.6e-23 | 1.0003 | 3.9e-34 |
| distance to trench edge (km) | 2077 | 0.141 | 1.2e-10 | 1.0010 | 6.7e-06 |

Logistic `bab_zone ~ age + time dummies` (age coefficient after time as a
factor): n=1831 (n_on=177); age coef = 0.0016 per Ma (OR=1.0016); LRT vs time-only χ²=0.63, p=0.4257

### 4. 0 Ma SubMap sanity (bab_on vs UPS_3, E vs C)

Join of 0 Ma `bab_on` onto `results/joined_submap_muller2019_0Ma.csv`
(nearest trench sample to each SubMap Müller hit; median distance
0.00 km).

**E vs C only** (`UPS_3`):

| UPS_3 | n | n bab_on | P(bab_on) |
|---|---:|---:|---:|
| C | 36 | 0 | 0.00 |
| E | 39 | 12 | 0.31 |

Crosstab `UPS_3` × `bab_on` (C/E):

```
bab_on   0   1
UPS_3         
C       36   0
E       27  12
```

This is a **simple agreement table**, not a winner claim. Layer A is MOR
back-arc spreading in Müller2019; SubMap `UPS_3` is geodetic strain.
Hellenic/Aegean is geodetic E but `bab_on=0` because the model has no
Aegean MOR (dilatation is the relevant layer there). Andes is geodetic C
and `bab_on=0`.

Full Geodetic_UPS × bab_on:

| Geodetic_UPS | n | n bab_on | P(bab_on) |
|---|---:|---:|---:|
| C | 7 | 0 | 0.00 |
| C-SS | 29 | 0 | 0.00 |
| E | 7 | 2 | 0.29 |
| E-SS | 32 | 10 | 0.31 |
| Neutral | 16 | 1 | 0.06 |
| SS | 32 | 1 | 0.03 |
| SS-C | 50 | 0 | 0.00 |
| SS-E | 38 | 12 | 0.32 |
| Undefined | 30 | 2 | 0.07 |


## Per-time summary

| time_Ma | n | n bab_on | frac | med age on | med age off | P(on\|≥55) | P(on\|<55) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1238 | 120 | 0.097 | 97.5 | 44.0 | 0.247 | 0.023 |
| 1 | 1280 | 119 | 0.093 | 97.6 | 44.0 | 0.251 | 0.019 |
| 2 | 1281 | 114 | 0.089 | 96.5 | 45.4 | 0.219 | 0.019 |
| 3 | 1297 | 110 | 0.085 | 93.1 | 45.7 | 0.218 | 0.017 |
| 4 | 1336 | 75 | 0.056 | 87.2 | 46.8 | 0.129 | 0.018 |
| 5 | 1297 | 71 | 0.055 | 91.0 | 45.1 | 0.124 | 0.019 |
| 6 | 1310 | 37 | 0.028 | 69.6 | 53.3 | 0.041 | 0.024 |
| 7 | 1287 | 35 | 0.027 | 69.4 | 50.4 | 0.040 | 0.024 |
| 8 | 1282 | 34 | 0.027 | 69.2 | 57.1 | 0.036 | 0.022 |
| 9 | 1297 | 31 | 0.024 | 68.2 | 66.2 | 0.031 | 0.021 |
| 10 | 1271 | 32 | 0.025 | 67.1 | 64.3 | 0.033 | 0.021 |
| 11 | 1268 | 30 | 0.024 | 68.4 | 62.1 | 0.029 | 0.014 |
| 12 | 1260 | 17 | 0.013 | 48.6 | 63.6 | 0.011 | 0.023 |
| 13 | 1275 | 21 | 0.016 | 47.8 | 65.7 | 0.009 | 0.041 |
| 14 | 1271 | 21 | 0.017 | 46.4 | 64.4 | 0.011 | 0.037 |
| 15 | 1298 | 95 | 0.073 | 76.5 | 58.7 | 0.116 | 0.079 |
| 16 | 1300 | 74 | 0.057 | 79.5 | 57.9 | 0.100 | 0.043 |
| 17 | 1329 | 80 | 0.060 | 69.2 | 60.7 | 0.094 | 0.053 |
| 18 | 1350 | 79 | 0.059 | 62.3 | 57.0 | 0.096 | 0.050 |
| 19 | 1347 | 76 | 0.056 | 52.0 | 55.2 | 0.064 | 0.080 |
| 20 | 1342 | 73 | 0.054 | 56.8 | 54.5 | 0.070 | 0.065 |
| 21 | 1306 | 48 | 0.037 | 65.5 | 54.3 | 0.059 | 0.035 |
| 22 | 1300 | 48 | 0.037 | 85.4 | 52.1 | 0.074 | 0.021 |
| 23 | 1291 | 46 | 0.036 | 67.4 | 54.6 | 0.053 | 0.037 |
| 24 | 1234 | 44 | 0.036 | 69.5 | 58.1 | 0.079 | 0.007 |
| 25 | 1240 | 54 | 0.044 | 59.8 | 54.7 | 0.085 | 0.023 |
| 26 | 1176 | 78 | 0.066 | 51.9 | 54.2 | 0.059 | 0.106 |
| 27 | 1184 | 83 | 0.070 | 49.7 | 52.9 | 0.059 | 0.114 |
| 28 | 1189 | 62 | 0.052 | 53.1 | 49.7 | 0.058 | 0.068 |
| 29 | 1209 | 73 | 0.060 | 54.7 | 48.5 | 0.076 | 0.063 |
| 30 | 1249 | 83 | 0.066 | 55.3 | 46.8 | 0.095 | 0.064 |
| 31 | 1261 | 60 | 0.048 | 62.6 | 47.0 | 0.071 | 0.042 |
| 32 | 1235 | 70 | 0.057 | 51.6 | 48.3 | 0.065 | 0.065 |
| 33 | 1225 | 70 | 0.057 | 49.8 | 47.8 | 0.065 | 0.044 |
| 34 | 1215 | 81 | 0.067 | 60.2 | 46.8 | 0.076 | 0.046 |
| 35 | 1232 | 60 | 0.049 | 92.4 | 46.1 | 0.059 | 0.017 |
| 36 | 1178 | 46 | 0.039 | 119.5 | 45.6 | 0.069 | 0.000 |
| 37 | 1150 | 46 | 0.040 | 111.0 | 44.2 | 0.074 | 0.000 |
| 38 | 1144 | 55 | 0.048 | 113.5 | 44.3 | 0.089 | 0.000 |
| 39 | 1140 | 58 | 0.051 | 109.3 | 41.9 | 0.092 | 0.004 |
| 40 | 1176 | 53 | 0.045 | 107.8 | 39.4 | 0.103 | 0.000 |
| 41 | 1120 | 59 | 0.053 | 105.0 | 40.9 | 0.094 | 0.002 |
| 42 | 1116 | 62 | 0.056 | 103.0 | 41.9 | 0.109 | 0.002 |
| 43 | 1146 | 57 | 0.050 | 103.3 | 42.2 | 0.122 | 0.000 |
| 44 | 1140 | 58 | 0.051 | 100.8 | 40.7 | 0.137 | 0.000 |
| 45 | 1147 | 68 | 0.059 | 97.9 | 36.4 | 0.161 | 0.003 |
| 46 | 1214 | 89 | 0.073 | 90.2 | 35.9 | 0.211 | 0.003 |
| 47 | 1198 | 81 | 0.068 | 90.1 | 32.5 | 0.210 | 0.000 |
| 48 | 1285 | 63 | 0.049 | 97.7 | 28.3 | 0.111 | 0.003 |
| 49 | 1283 | 63 | 0.049 | 94.2 | 28.6 | 0.107 | 0.005 |
| 50 | 1295 | 62 | 0.048 | 94.6 | 26.9 | 0.099 | 0.002 |
| 51 | 1225 | 39 | 0.032 | 16.9 | 35.1 | 0.043 | 0.030 |
| 52 | 1327 | 38 | 0.029 | 15.6 | 39.8 | 0.042 | 0.026 |
| 53 | 1330 | 95 | 0.071 | 28.5 | 41.3 | 0.091 | 0.086 |
| 54 | 1335 | 97 | 0.073 | 38.0 | 40.6 | 0.086 | 0.088 |
| 55 | 1338 | 97 | 0.072 | 35.2 | 41.1 | 0.082 | 0.094 |
| 56 | 1289 | 69 | 0.054 | 53.3 | 44.0 | 0.082 | 0.063 |
| 57 | 1308 | 90 | 0.069 | 60.0 | 45.7 | 0.128 | 0.061 |
| 58 | 1295 | 36 | 0.028 | 42.4 | 47.5 | 0.000 | 0.060 |
| 59 | 1292 | 36 | 0.028 | 40.4 | 46.0 | 0.000 | 0.058 |
| 60 | 1296 | 37 | 0.029 | 38.5 | 42.7 | 0.000 | 0.058 |

## How to rerun

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \
    /workspace/subduction-test/scripts/cenozoic_ups_proxy.py
```

0 Ma prototype only:

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \
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
