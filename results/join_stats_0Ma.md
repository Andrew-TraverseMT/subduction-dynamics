# 0 Ma SubMap UPS ↔ Müller 2019 trench-kinematics join

First-milestone look only. Values are computed from the joined table; nothing is invented.

- **Kinematics:** `/workspace/subduction-test/results/trench_kinematics_0Ma_Muller2019.csv`
- **SubMap:** `/workspace/subduction-test/data/submap/submap_geodetic_UPS.csv` (n=260)
- **Joined table:** `/workspace/subduction-test/results/joined_submap_muller2019_0Ma.csv`
- **Join:** nearest Müller 2019 trench sample to each SubMap (Lon, Lat), great-circle distance with lon wrapping (SubMap lon is 0–360 in some zones; Müller is −180–180). Earth radius 6371.009 km.
- **Cutoff:** drop joins farther than 200 km.

## Match rates

- n SubMap transects: **260**
- n matched (≤ 200 km): **241**
- n dropped (> 200 km): **19**
- median join distance (matched): **25.46 km**
- mean join distance (matched): **43.14 km**
- max join distance (matched): **198.45 km**

Dropped transects (join > 200 km):

| Short_Name | Zone | Trench_Name | Lon | Lat | dist_km | Geodetic_UPS |
|---|---|---|---:|---:|---:|---|
| MED14 | MED | Makran | 57.50 | 24.70 | 276.4 | C-SS |
| MED15 | MED | Makran | 59.50 | 24.20 | 249.8 | C-SS |
| SAM18 | SAM | Venezuela | 287.00 | 13.40 | 263.6 | SS-C |
| SAM19 | SAM | Venezuela | 289.00 | 13.80 | 477.9 | SS |
| SAM20 | SAM | Venezuela | 291.00 | 13.50 | 675.9 | SS |
| SAM21 | SAM | Venezuela | 293.00 | 13.10 | 883.3 | SS |
| SAM22 | SAM | Muertos | 290.00 | 17.30 | 789.0 | SS-C |
| SAM23 | SAM | Muertos | 292.00 | 17.20 | 761.0 | SS-C |
| SAM24 | SAM | Muertos | 294.00 | 17.20 | 568.0 | SS-C |
| SAM25 | SAM | Hispaniola | 290.00 | 20.00 | 917.8 | SS-C |
| SAM26 | SAM | Puerto-Rico | 292.00 | 19.60 | 708.3 | SS-C |
| SAM27 | SAM | Puerto-Rico | 294.00 | 19.80 | 499.1 | SS |
| SAM28 | SAM | Antilles | 296.00 | 19.80 | 290.1 | SS |
| SEA19 | SEA | Java | 117.00 | -11.30 | 205.6 | SS-C |
| SEA21 | SEA | Timor | 121.00 | -12.10 | 311.6 | C-SS |
| SEA22 | SEA | Timor | 123.00 | -11.20 | 223.9 | C-SS |
| SEA33 | SEA | Seram | 130.00 | -2.40 | 336.5 | SS-C |
| SWP32 | SWP | Hikurangi | 175.00 | -42.30 | 343.5 | SS-C |
| SWP33 | SWP | Puysegur | 164.50 | -48.00 | 1378.0 | Undefined |

## UPS recode (`UPS_3`)

| Geodetic_UPS | UPS_3 | n matched |
|---|---|---:|
| C, C-SS | C | 36 |
| Neutral | N | 16 |
| E, E-SS | E | 39 |
| SS, SS-C, SS-E | SS | 120 |
| Undefined | Undefined (out of C/N/E test) | 30 |

Geodetic_UPS counts among matched rows:

| class | n |
|---|---:|
| C | 7 |
| C-SS | 29 |
| SS-C | 50 |
| SS | 32 |
| SS-E | 38 |
| E-SS | 32 |
| E | 7 |
| Neutral | 16 |
| Undefined | 30 |

M56 velocities were converted mm/yr → cm/yr (`*_cm_yr` columns) so they sit next to Müller cm/yr. SubMap sentinel −999 was set to NaN before conversion. `Seismic_UPS` is all −999 and is ignored.

## Spearman rank correlations (defined UPS only)

Defined UPS = `Geodetic_UPS != Undefined`, plus finite values of both series. Mixed strike-slip (UPS_3 = SS) **is included** here (it is a defined class); Undefined is not.

| pair | Spearman ρ | p | n |
|---|---:|---:|---:|
| `v_OP_perp` vs SubMap `M56_vdn` (cm/yr) | 0.193 | 0.0050 | 211 |
| `v_T_perp` vs SubMap `M56_vdn` (cm/yr) | -0.183 | 0.0076 | 211 |
| `seafloor_age_Ma` vs SubMap `M56_vdn` (cm/yr) | -0.295 | 0.0002 | 153 |

Sign conventions (do not mix them up):

- SubMap: shortening ⇒ 0 ≤ Ovd < 90°, `vdn > 0`; extension ⇒ 90° < Ovd ≤ 180°, `vdn < 0`.
- Müller `v_OP_perp > 0` = overriding plate moving **away from the trench** (upper-plate retreat).
- Müller `v_T_perp > 0` = **oceanward trench retreat**.

Because the signs differ, ρ(`v_OP_perp`, `M56_vdn`) need not be positive even if the two fields describe the same deformation sense (OP retreat ↔ extension ↔ `vdn < 0`). The observed sign of ρ is reported as computed; it is not a hypothesis verdict.

## Median kinematics by UPS_3

| UPS_3 | median `v_OP_perp` | median `v_T_perp` |
|---|---|---|
| C (C + C-SS) | -0.575 cm/yr (n=36) | 0.564 cm/yr (n=36) |
| N (Neutral) | -0.868 cm/yr (n=16) | 0.868 cm/yr (n=16) |
| E (E + E-SS) | -1.980 cm/yr (n=39) | 1.980 cm/yr (n=39) |
| SS (SS + SS-C + SS-E), shown separately | -0.733 cm/yr (n=120) | 0.733 cm/yr (n=120) |

C / N / E are the classes for a first H1/H2 look. SS is reported so it is not silently folded in.

## `v_T_perp` and `v_OP_perp` are not independent in Müller 2019

In this rigid/deforming-network model the trench line is typically attached to `trench_plate_id`, so **`v_T_perp ≈ −v_OP_perp`**. Among matched rows:

- Spearman ρ(`v_T_perp`, `v_OP_perp`) = -0.998 (n=241)
- median |`v_T_perp` + `v_OP_perp`| = 0.0001 cm/yr

Correlations of `v_T_perp` vs `M56_vdn` and of `v_OP_perp` vs `M56_vdn` are therefore **not independent tests of H1 (trench rollback) vs H2 (upper-plate motion)**. They are approximately sign-flipped versions of the same kinematic field. Do not treat a “winner” from this pair.

## First-look caveat (do not overclaim)

This is a nearest-point geographic join of two independently built 0 Ma products, not a common segmentation. UPS class is SubMap’s geodetic classification (from their `Ovd` / `vd`); the top panel of `fig_classic_replication.png` therefore largely recovers SubMap’s own definition of UPS, and is a sanity check rather than an independent hypothesis test. The Müller panels ask whether reconstruction kinematics line up with that classification. Sample sizes in pure C and pure E are small. No hypothesis is declared a winner here.

## Figures

- `/workspace/subduction-test/results/fig_vOP_vs_UPS.png`
- `/workspace/subduction-test/results/fig_vT_vs_UPS.png`
- `/workspace/subduction-test/results/fig_classic_replication.png`
- `/workspace/subduction-test/results/fig_age_vs_UPS.png`
- `/workspace/subduction-test/results/fig_vOP_vs_submap_vdn.png`
