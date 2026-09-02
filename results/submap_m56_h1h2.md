# SubMap M56 H1 vs H2 kinematic test (0 Ma)

Values are computed from the SubMap dump and, where used, the existing Müller join.
Nothing is invented. UPS ~ vdn is circular by construction and is reported only as a sanity check.

- **Compact UPS table:** `/workspace/subduction-test/data/submap/submap_geodetic_UPS.csv` (n=260)
- **Full Sub-DATA sheet (absolute velocities only):** `/workspace/subduction-test/data/submap/submap_transects_full.csv`
- **Script:** `/workspace/subduction-test/scripts/submap_m56_h1h2.py`
- **Units:** SubMap stores mm/yr; all velocities below are converted to **cm/yr** (÷ 10). Sentinel −999 → NaN.

## Sample

- n compact rows: **260**
- n dropped (Undefined and/or `M56_vd` missing/−999): **31**
- n kept: **229**
- C vs E test n: **79** (C n=40, E n=39)
- C/N/E ordinal n: **95** (N n=16)
- SS held out of C/N/E tests, shown separately in box plots and median tables (n=134)

### UPS_3 recode

C = C, C-SS; E = E, E-SS; N = Neutral; SS = SS, SS-C, SS-E (kept separate). Undefined dropped.

| Geodetic_UPS | UPS_3 | n |
|---|---|---:|
| C | C | 7 |
| C-SS | C | 33 |
| Neutral | N | 16 |
| E | E | 7 |
| E-SS | E | 32 |
| SS | SS | 37 |
| SS-C | SS | 59 |
| SS-E | SS | 38 |

`UPN` is a 1/2 upper-plate-nature code, **not** a velocity. It is not used as an H1/H2 proxy.
Among kept rows, UPN counts by UPS_3 (not used as a predictor):

| UPS_3 | UPN=1 | UPN=2 |
|---|---:|---:|
| C | 38 | 2 |
| N | 4 | 12 |
| E | 14 | 25 |
| SS | 104 | 30 |

## Kinematic identity (documented, not assumed)

Cerpa, Lallemand & Heuret (2025, §II.1) define **relative** velocities from MORVEL56 + an arc block:

- `vc` = convergence of undeformed SP and UP interiors
- `vd` = arc-block motion relative to UP (forearc-to-back-arc deformation)
- `vs` = subduction / SP consumption at the trench (`vs = V_SP + V_trench` as vectors)
- Vector statement in the paper: `vs = vc + vd`

SubMap **signed trench-normal scalars** (shortening ⇒ `vdn > 0`, extension ⇒ `vdn < 0`) obey

**`vcn ≈ vsn + vdn`**  (equivalently **`vsn ≈ vcn − vdn`**).

On the kept rows this is a data identity, not a fit:

- n finite: **229**
- median |vcn − (vsn + vdn)|: **0.000 mm/yr**
- n with exact 0 (integer mm/yr): **149**
- n with |residual| ≤ 1 mm/yr: **223**
- n with |residual| > 1 mm/yr: **6** (max 21.0 mm/yr)

Rows with |residual| > 1 mm/yr (left as in the database; not corrected):

| Short_Name | Trench_Name | Geodetic_UPS | vcn | vsn | vdn | vcn−(vsn+vdn) |
|---|---|---|---:|---:|---:|---:|
| NPA23 | Aleutians | SS-C | 60.0 | 59.0 | 5.0 | -4.0 |
| SAM15 | Panama | C-SS | 32.0 | 9.0 | 27.0 | -4.0 |
| SAM52 | Andean | SS-C | 68.0 | 63.0 | 3.0 | 2.0 |
| SAM55 | Andean | SS-C | 59.0 | 53.0 | 3.0 | 3.0 |
| SAM56 | Andean | C-SS | 55.0 | 49.0 | 4.0 | 2.0 |
| SEA61 | Ryukyus | SS-E | 32.0 | 66.0 | -55.0 | 21.0 |

Absolute trench-normal scalars in the full sheet (frame 1 = MORVEL56-NNR) satisfy

**`vdn ≈ Vupn1 − Vtn1`** (median |residual| = 0.000 mm/yr; 222/229 within 1 mm/yr)

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

Spearman ρ(`Vupn1`, `Vtn1`) = 0.655 (p = 1.79e-29, n = 229).
That is **not** Müller-style collinearity (`v_T ≈ −v_OP`). SubMap’s own absolute
velocities therefore *can* separate H1 from H2, but only after merging `Vupn1`/`Vtn1`
from the full sheet. The compact CSV alone cannot.

If those absolute columns are ignored, the compact table still supports a
**non-circular relative-velocity test**: does `vcn` or `vsn` predict UPS_3 (C vs E)
better? That is not a clean H1-vs-H2 split (`vcn` mixes both plates; `vsn` mixes SP
and trench), but it does not use the class-defining `vdn`.

## Median M56 velocities by UPS_3 (cm/yr)

| UPS_3 | n | median vcn | median vsn | median vdn |
|---|---:|---|---|---|
| C | 40 | 5.700 cm/yr (n=40) | 4.300 cm/yr (n=40) | 1.400 cm/yr (n=40) |
| N | 16 | 1.600 cm/yr (n=16) | 1.600 cm/yr (n=16) | 0.000 cm/yr (n=16) |
| E | 39 | 4.500 cm/yr (n=39) | 6.600 cm/yr (n=39) | -3.100 cm/yr (n=39) |
| SS | 134 | 4.900 cm/yr (n=134) | 4.450 cm/yr (n=134) | 0.100 cm/yr (n=134) |

Absolute-velocity medians (same rows; cm/yr):

| UPS_3 | median Vupn1 (H1) | median Vtn1 (H2) |
|---|---|---|
| C | 0.900 cm/yr (n=40) | -0.250 cm/yr (n=40) |
| N | -0.600 cm/yr (n=16) | -0.600 cm/yr (n=16) |
| E | -0.900 cm/yr (n=39) | 2.200 cm/yr (n=39) |
| SS | 0.650 cm/yr (n=134) | 0.800 cm/yr (n=134) |

## Spearman: C/N/E ordinal and C vs E

Score sign: C = +1, E = −1 (and N = 0), matching SubMap `vdn` (shortening positive).
`vdn` vs class is **circular**. `vcn` and `vsn` are the non-circular relative predictors.
`Vupn1` / `Vtn1` are the H1 / H2 absolute predictors.

| predictor | Spearman ρ vs C/N/E (C=1, N=0, E=−1) | p | n |
|---|---:|---:|---:|
| vcn | 0.301 | 0.0030 | 95 |
| vsn | -0.385 | 0.0001 | 95 |
| vdn | 0.925 | 6.82e-41 | 95 |
| Vupn1 | 0.460 | 2.78e-06 | 95 |
| Vtn1 | -0.299 | 0.0033 | 95 |

| predictor | Spearman ρ vs C vs E (C=1, E=−1) | p | n |
|---|---:|---:|---:|
| vcn | 0.314 | 0.0049 | 79 |
| vsn | -0.471 | 1.16e-05 | 79 |
| vdn | 0.866 | 6.43e-25 | 79 |
| Vupn1 | 0.498 | 2.99e-06 | 79 |
| Vtn1 | -0.341 | 0.0021 | 79 |

## OLS and logistic (C vs E only)

Binary y = 1 if C, 0 if E. Predictors in cm/yr. Nested model is `y ~ vcn + vsn`.
`vdn` is listed as a sanity check (quasi-complete separation is expected).

| predictor | OLS slope (y=1 if C, 0 if E) | R² | p (slope) | n |
|---|---:|---:|---:|---:|
| vcn | 0.066 | 0.125 | 0.0014 | 79 |
| vsn | -0.058 | 0.226 | 9.45e-06 | 79 |
| vdn | 0.098 | 0.566 | 1.34e-15 | 79 |
| Vupn1 | 0.070 | 0.236 | 5.81e-06 | 79 |
| Vtn1 | -0.048 | 0.118 | 0.0020 | 79 |
| vcn + vsn (nested) | — | 0.592 | — | 79 |

| predictor | McFadden R² | AIC | converged | quasi-separation | n |
|---|---:|---:|---|---|---:|
| vcn | 0.094 | 103.2 | yes | no | 79 |
| vsn | 0.213 | 90.1 | yes | no | 79 |
| vdn (circular sanity) | 1.000 | 4.0 | yes | yes | 79 |
| Vupn1 | 0.206 | 91.0 | yes | no | 79 |
| Vtn1 | 0.094 | 103.2 | yes | no | 79 |
| vcn + vsn (nested) | 1.000 | 6.0 | yes | yes | 79 |

Because `vdn ≈ vcn − vsn`, the two-predictor nested model reconstructs the
class-defining deformation rate. Nested logistic McFadden R² = 1 and nested OLS R²
near the univariate `vdn` R² are therefore **circular**, not evidence that `vcn` and
`vsn` jointly explain UPS independently of `vd`.

## Which predictor is stronger (plain language)

Non-circular relative velocities, C vs E (n=79):
Spearman ρ(vcn) = 0.314, ρ(vsn) = -0.471;
OLS R²(vcn) = 0.125, R²(vsn) = 0.226;
logistic McFadden R²(vcn) = 0.094, R²(vsn) = 0.213.
|ρ(vsn, C/E)| = 0.471 exceeds |ρ(vcn, C/E)| = 0.314 (gap 0.158). Univariate comparison only: nested `vcn+vsn` reconstructs `vdn` via the identity, so nested R² values are not an independent test (OLS R²=0.592; logistic McFadden R²=1.000).

`vdn` vs C/E is circular (ρ = 0.866) and is not a test.

H1 vs H2 from SubMap absolute velocities: SubMap Vupn1 and Vtn1 are only moderately rank-correlated, unlike Müller v_OP ≈ −v_T. On C vs E, |ρ| is 0.498 for Vupn1 (H1) and 0.341 for Vtn1 (H2). That is not an overwhelming H1-vs-H2 separation.

No hypothesis is declared a winner: the vcn-vs-vsn contrast is not overwhelming on this sample.

## Müller 2019 covariates (H2/H3; independent of SubMap class)

Joined table: `/workspace/subduction-test/results/joined_submap_muller2019_0Ma.csv` (nearest trench sample, ≤ 200 km; n matched = 241).
These fields are **not** derived from SubMap UPS or M56 velocities.

- n with defined UPS_3 after join: **211**
- n with finite seafloor age: **153**

| predictor | sample | Spearman ρ vs UPS score | p | n |
|---|---|---:|---:|---:|
| `seafloor_age_Ma` | C/N/E (C=1,N=0,E=−1) | -0.418 | 0.0006 | 64 |
| `seafloor_age_Ma` | C vs E (C=1,E=−1) | -0.502 | 0.0002 | 50 |
| `trench_zone_width_km` | C/N/E (C=1,N=0,E=−1) | -0.167 | 0.1129 | 91 |
| `trench_zone_width_km` | C vs E (C=1,E=−1) | -0.189 | 0.1035 | 75 |
| `dist_nearest_trench_edge_km` | C/N/E (C=1,N=0,E=−1) | -0.128 | 0.2273 | 91 |
| `dist_nearest_trench_edge_km` | C vs E (C=1,E=−1) | -0.132 | 0.2592 | 75 |

| UPS_3 | median age (Ma) | median zone width (km) | median dist-to-edge (km) |
|---|---|---|---|
| C | 38.289 Ma (n=18) | 4166.426 km (n=36) | 388.451 km (n=36) |
| N | 80.741 Ma (n=14) | 815.964 km (n=16) | 322.145 km (n=16) |
| E | 93.596 Ma (n=32) | 4365.731 km (n=39) | 555.311 km (n=39) |
| SS | 51.589 Ma (n=89) | 3860.945 km (n=120) | 515.762 km (n=120) |

Age, trench-zone width, and distance-to-edge are the independent H2/H3-style covariates
(slab age; slab/trench width; proximity to a slab edge). They are reported here against
SubMap UPS_3; they are not a test of SubMap M56 kinematics.

## Figures

- `/workspace/subduction-test/results/fig_submap_vcn_vs_UPS.png`
- `/workspace/subduction-test/results/fig_submap_vsn_vs_UPS.png`
- `/workspace/subduction-test/results/fig_submap_vdn_vs_UPS.png` (sanity; circular)
- `/workspace/subduction-test/results/fig_submap_Vupn_Vtn_vs_UPS.png` (H1/H2 absolute velocities from the full sheet)
- `/workspace/subduction-test/results/fig_age_width_vs_UPS.png`

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
