# Schellart width / edge test — Cenozoic, zone–time blocked

Built by `scripts/schellart_width_edge.py` from
`/workspace/subduction-test/results/trench_kinematics_0-60Ma_with_ups_proxy.csv`. Trench kinematics were **not** re-extracted.
Müller `v_T` vs `v_OP` is **not** used (collinear in this model).

- **n trench samples:** 76529
- **n times:** 61 (0–60 Ma, 1 Myr)
- **n zone–time blocks:** 2077
- **n bab_zone=1:** 187 (9.0%)
- **n ext_zone=1:** 442 (21.3%)
- **n blocks with a deforming-network sample:** 1084
- **n blocks with finite median age:** 1831

## Circularity (read this first)

Layer A (`bab_on` / `bab_zone`) is a back-arc MOR flag taken from the **same**
Müller et al. (2019) topologies that supply the trench catalog. Layer B is
dilatation inside those deforming networks. Both are reconstruction-internal.
The Cenozoic test asks whether *this model's* back-arc / extension geometry
follows Schellart et al. (2007) width and edge patterns, not whether nature
does. Present-day SubMap geodetic UPS remains the only UPS here that is
independent of the plate model; that 0 Ma table is a sanity check.

## Schellart prediction under test

Schellart, Freeman, Stegman, Moresi & May (2007, *Nature*): rollback and
back-arc opening concentrate on **narrow** slabs and **near a lateral slab
edge**; interiors of wide slabs (≳4000 km) are slow or advancing / more
compressive. They associated near-edge retreat with a length scale of order
**1200 km**.

Operationally, at zone–time:

- narrower `trench_zone_width_km` → higher P(bab_zone) and higher P(ext_zone)
  → Spearman ρ(**width**, extension flag) should be **negative**
- smaller `dist_nearest_trench_edge_km` → more extension
  → Spearman ρ(**edge distance**, extension flag) should be **negative**
- near-edge split: P(on | median edge < 1200 km) > P(on | ≥ 1200 km)

`trench_zone_width_km` is the along-trench length of a spatially connected
polyline sharing a subducting-plate ID (see `extract_trench_kinematics_0Ma.py`).
It is a trench-zone width, not a tomographic slab width. `zone_id` is unique
only within a time slice.

Zone-level median edge distance is mechanically related to width (a uniform
sample on a zone of length W has median distance-to-edge ≈ W/4).
Spearman ρ(width, edge) at zone–time = 0.729 (p=<1e-16, n=2077).
Treat the two predictors as overlapping geometric information.

## Inference unit and outcomes

**Unit:** `groupby(time_Ma, zone_id)` — **2077** blocks. Point-level n≈76k
is spatially autocorrelated along the trench and is not used for inference.

**Layer A — `bab_zone`.** 1 if any tessellated sample in the block has `bab_on=1`
(qualifying overriding-plate MOR; see `cenozoic_age_gate.md`). 0 otherwise.

**Layer B — `ext_zone`.** Among samples with `dilat_in_network==1`, take the
median `dilat_rate_1e-15_s`. `ext_zone = 1` if that median **> 0**
(threshold: strictly positive dilatation, units 10⁻¹⁵ s⁻¹). If the block has
no network samples (all rigid plates, dilatation recorded as zero /
`dilat_in_network==0`), `ext_zone = 0`. Negative median dilatation (network
shortening) is also `ext_zone = 0`. Hellenic / Aegean extension is a deforming
network in Müller2019, not a MOR: it should be A=0, B=1.

A network-only sensitivity drops `all_rigid` blocks so rigid zeros do not
dilute layer B.

Age is a covariate already known to be weak at this blocking
(`cenozoic_age_gate.md`); it is not the headline.

## 0 Ma region check (A vs B)

| region | n pts | n bab_on | bab_zone | ext_zone | med width (km) | med edge (km) | med dilat net |
|---|---:|---:|---:|---:|---:|---:|---:|
| Mariana | 40 | 39 | 1 | 0 | 5107 | 306 | 0.000 |
| Tonga | 32 | 28 | 1 | 0 | 3755 | 403 | -2.149 |
| Scotia | 17 | 17 | 1 | 0 | 929 | 222 | NA |
| Hellenic | 34 | 0 | 0 | 1 | 7890 | 556 | 0.123 |
| Andes | 66 | 0 | 0 | 1 | 6885 | 555 | 0.003 |
| Cascadia | 22 | 0 | 0 | 0 | 1290 | 311 | -0.143 |
| Andaman | 17 | 17 | 1 | 0 | 1460 | 494 | -0.010 |
| MAT_Mexico | 17 | 0 | 0 | 1 | 439 | 106 | 0.001 |

Expected: Mariana / Tonga / Scotia / Andaman **bab_zone=1**. Andes / Cascadia /
MAT_Mexico **bab_zone=0**. Hellenic **bab_zone=0** and **ext_zone=1**.
The `> 0` cut is literal: Andes and MAT_Mexico can land on `ext_zone=1` with
median network dilatation of order 10⁻³ (10⁻¹⁵ s⁻¹), which is not
Hellenic-style extension (Hellenic median is larger; see table). Continuous
median dilatation (below) does not use that binary cut.

0 Ma zone–time blocks: n=47, bab_zone=1 in 5, ext_zone=1 in 9.

## Layer A — MOR back-arc (`bab_zone`)

n_blocks = **2077**, n_on = **187**.

Median width bab_zone=1 vs 0: 3022 vs 1034 km
(U=253647.0, p=8.04e-23, n_on=187, n_off=1890).

Median edge distance bab_zone=1 vs 0: 445 vs 307 km
(U=226832.5, p=1.48e-10).

### Spearman / univariate logistic (zone–time)

| predictor | n | n_on | Spearman ρ | p | logistic OR (no time) | LRT p |
|---|---:|---:|---:|---:|---:|---:|
| trench-zone width (km) | 2077 | 187 | 0.216 | 2.57e-23 | 1.0003 | 3.91e-34 |
| distance to trench edge (km) | 2077 | 187 | 0.141 | 1.22e-10 | 1.0010 | 6.68e-06 |
| seafloor age (Ma) | 1831 | 177 | 0.042 | 0.0717 | 1.0014 | 0.4661 |

Width vs bab_zone: blocked numbers are not a Schellart confirmation (ρ=0.216, p=2.57e-23; sign is opposite Schellart (wider / farther-from-edge more extensional); ). No winner claimed.
Edge vs bab_zone: blocked numbers are not a Schellart confirmation (ρ=0.141, p=1.22e-10; sign is opposite Schellart (wider / farther-from-edge more extensional); ). No winner claimed.
Age covariate: ρ=0.042, p=0.0717 (n=1831).

### Logistic with time dummies

`y ~ predictor + time factor`. Coefficient is after absorbing time.

| predictor | n | n_on | coef | OR per unit | LRT vs time-only χ² | p |
|---|---:|---:|---:|---:|---:|---:|
| trench-zone width (km) | 2077 | 187 | 0.00027 | 1.00027 | 145.04 | 2.10e-33 |
| distance to trench edge (km) | 2077 | 187 | 0.00132 | 1.00132 | 26.64 | 2.45e-07 |
| seafloor age (Ma) | 1831 | 177 | 0.00165 | 1.00165 | 0.63 | 0.4257 |

Joint `bab_zone ~ width + edge + age + time dummies`: n=1831 (n_on=177); width coef=0.00026 (OR=1.00026), edge coef=0.00040 (OR=1.00040), age coef=-0.00156 (OR=0.99844); LRT vs time-only χ²=140.30 (df=3, p=3.26e-30).
Width and edge are collinear; joint coefficients are descriptive only.

### Quantile splits (Schellart geometry)

Width tertiles are pooled rank tertiles of zone–time median width (average-rank
ties). Narrow vs wide MH drops the middle tertile. Edge split is median
distance < 1200 km vs ≥ 1200 km.

| tertile | n blocks | n on | P(on) | width min–max (km) | median width (km) |
|---|---:|---:|---:|---|---:|
| narrow | 692 | 31 | 0.045 | 2–573 | 226 |
| mid | 692 | 36 | 0.052 | 574–1882 | 1108 |
| wide | 693 | 120 | 0.173 | 1882–21118 | 3272 |

Mantel–Haenszel OR (narrow T1 vs wide T3, y=1 more likely if exposed) = 0.204 (χ²=59.96, p=9.69e-15, 61 time strata); P(y=1|narrow T1)=0.045 (n=692, n_on=31) vs P(y=1|wide T3)=0.173 (n=693, n_on=120)

Edge <1200 km: P(bab_zone)=0.091 (n=2048, n_on=187);
≥1200 km: P(bab_zone)=0.000 (n=29, n_on=0).

Mantel–Haenszel near <1200 km vs far ≥1200 km: not defined (n_exp=2048, n_un=29)

| tertile | n blocks | n on | P(on) | edge min–max (km) | median edge (km) |
|---|---:|---:|---:|---|---:|
| near T1 | 678 | 39 | 0.058 | 0–222 | 101 |
| mid T2 | 684 | 46 | 0.067 | 222–445 | 317 |
| far T3 | 715 | 102 | 0.143 | 445–2222 | 608 |

An MH OR > 1 for *narrow* vs *wide* (or *near* vs *far*) would match Schellart
(exposed class more extensional).

## Layer B — dilatation (`ext_zone`)

Threshold: median network dilatation **> 0** × 10⁻¹⁵ s⁻¹.
n_blocks = **2077**, n_on = **442**, of which network blocks = 1084.

Median width ext_zone=1 vs 0: 1347 vs 1066 km
(U=392527.5, p=0.0053, n_on=442, n_off=1635).

Median edge distance ext_zone=1 vs 0: 336 vs 314 km
(U=371505.5, p=0.3632).

### Spearman / univariate logistic (all blocks; rigid → ext_zone=0)

| predictor | n | n_on | Spearman ρ | p | logistic OR (no time) | LRT p |
|---|---:|---:|---:|---:|---:|---:|
| trench-zone width (km) | 2077 | 442 | 0.061 | 0.0053 | 1.0000 | 0.1690 |
| distance to trench edge (km) | 2077 | 442 | 0.020 | 0.3633 | 1.0001 | 0.7206 |
| seafloor age (Ma) | 1831 | 373 | -0.059 | 0.0116 | 0.9988 | 0.3901 |

Width vs ext_zone: blocked numbers are not a Schellart confirmation (ρ=0.061, p=0.0053; sign is opposite Schellart (wider / farther-from-edge more extensional); ). No winner claimed.
Edge vs ext_zone: blocked numbers are not a Schellart confirmation (ρ=0.020, p=0.3633; not distinguishable from null at 0.05; ). No winner claimed.

### Logistic with time dummies

| predictor | n | n_on | coef | OR per unit | LRT vs time-only χ² | p |
|---|---:|---:|---:|---:|---:|---:|
| trench-zone width (km) | 2077 | 442 | 0.00003 | 1.00003 | 2.43 | 0.1189 |
| distance to trench edge (km) | 2077 | 442 | 0.00018 | 1.00018 | 0.83 | 0.3616 |
| seafloor age (Ma) | 1831 | 373 | -0.00215 | 0.99785 | 2.08 | 0.1493 |

Joint `ext_zone ~ width + edge + age + time dummies`: n=1831 (n_on=373); width coef=0.00004 (OR=1.00004), edge coef=-0.00005 (OR=0.99995), age coef=-0.00257 (OR=0.99743); LRT vs time-only χ²=5.87 (df=3, p=0.1182).

### Quantile splits

| tertile | n blocks | n on | P(on) | width min–max (km) | median width (km) |
|---|---:|---:|---:|---|---:|
| narrow | 692 | 152 | 0.220 | 2–573 | 226 |
| mid | 692 | 93 | 0.134 | 574–1882 | 1108 |
| wide | 693 | 197 | 0.284 | 1882–21118 | 3272 |

Mantel–Haenszel OR (narrow T1 vs wide T3, y=1 more likely if exposed) = 0.666 (χ²=10.44, p=0.0012, 61 time strata); P(y=1|narrow T1)=0.220 (n=692, n_on=152) vs P(y=1|wide T3)=0.284 (n=693, n_on=197)

Edge <1200 km: P(ext_zone)=0.212 (n=2048);
≥1200 km: P(ext_zone)=0.241 (n=29).

Mantel–Haenszel OR (near <1200 km vs far ≥1200 km, y=1 more likely if exposed) = 0.739 (χ²=0.43, p=0.5111, 61 time strata); P(y=1|near <1200 km)=0.212 (n=2048, n_on=435) vs P(y=1|far ≥1200 km)=0.241 (n=29, n_on=7)

### Continuous dilatation (network blocks; no binary cut)

Spearman of median network dilatation (signed, 10⁻¹⁵ s⁻¹) vs predictors.
Positive ρ means larger width / edge distance goes with more extension.

| predictor | n | Spearman ρ | p |
|---|---:|---:|---:|
| trench-zone width (km) | 1084 | 0.007 | 0.8057 |
| distance to trench edge (km) | 1084 | 0.037 | 0.2253 |
| seafloor age (Ma) | 927 | -0.012 | 0.7108 |

Schellart direction for this signed rate is still negative ρ (narrower /
nearer-edge → more positive dilatation).

### Sensitivity: network blocks only (drop all-rigid)

n=1084, n_on=442.

| predictor | n | n_on | Spearman ρ | p | logistic OR (no time) | LRT p |
|---|---:|---:|---:|---:|---:|---:|
| trench-zone width (km) | 1084 | 442 | -0.055 | 0.0687 | 0.9999 | 7.51e-05 |
| distance to trench edge (km) | 1084 | 442 | -0.031 | 0.3118 | 0.9997 | 0.0995 |
| seafloor age (Ma) | 927 | 373 | -0.041 | 0.2147 | 0.9992 | 0.5860 |

| tertile | n blocks | n on | P(on) | width min–max (km) | median width (km) |
|---|---:|---:|---:|---|---:|
| narrow | 329 | 152 | 0.462 | 2–571 | 220 |
| mid | 278 | 93 | 0.335 | 574–1882 | 1151 |
| wide | 477 | 197 | 0.413 | 1882–21118 | 3478 |

Mantel–Haenszel OR (narrow T1 vs wide T3, y=1 more likely if exposed) = 1.250 (χ²=2.22, p=0.1364, 61 time strata); P(y=1|narrow T1)=0.462 (n=329, n_on=152) vs P(y=1|wide T3)=0.413 (n=477, n_on=197)
Mantel–Haenszel OR (near <1200 km vs far ≥1200 km, y=1 more likely if exposed) = 2.364 (χ²=2.90, p=0.0885, 61 time strata); P(y=1|near <1200 km)=0.410 (n=1060, n_on=435) vs P(y=1|far ≥1200 km)=0.292 (n=24, n_on=7)

## 0 Ma SubMap sanity (width / edge vs UPS_3, independent of bab_on)

Join: `/workspace/subduction-test/results/joined_submap_muller2019_0Ma.csv` (nearest Müller trench sample, already computed).
Width and edge come from Müller2019; UPS_3 is SubMap geodetic strain.
This table does **not** use `bab_on`. SS held out of C vs E. Undefined dropped.
C = C+C-SS, E = E+E-SS. Score C=+1, E=−1.

- n defined UPS_3 after join: **211**
- C vs E n: **75** (C n=36, E n=39)

| predictor | Spearman ρ vs C vs E (C=+1, E=−1) | p | n |
|---|---:|---:|---:|
| trench-zone width (km) | -0.189 | 0.1035 | 75 |
| distance to trench edge (km) | -0.132 | 0.2592 | 75 |
| seafloor age (Ma) | -0.502 | 0.0002 | 50 |

| UPS_3 | n | median width (km) | median dist-to-edge (km) | median age (Ma) |
|---|---:|---|---|---|
| C | 36 | 4166.426 km (n=36) | 388.451 km (n=36) | 38.289 Ma (n=18) |
| N | 16 | 815.964 km (n=16) | 322.145 km (n=16) | 80.741 Ma (n=14) |
| E | 39 | 4365.731 km (n=39) | 555.311 km (n=39) | 93.596 Ma (n=32) |
| SS | 120 | 3860.945 km (n=120) | 515.762 km (n=120) | 51.589 Ma (n=89) |

Schellart at 0 Ma vs independent UPS: ρ(width, C=+1/E=−1) = -0.189 (p=0.1035, n=75). A Schellart width effect (wide → more C) would be **positive**. This 0 Ma width result is weak (not distinguishable from null at 0.05).
ρ(edge, C/E) = -0.132 (p=0.2592, n=75).
A Schellart edge effect (near-edge → more E) would be **positive** ρ vs C=+1
(farther from edge → more C). Do not treat this as a Cenozoic confirmation.

## Do width / edge survive blocking?

Blocked n = 2077 zone–time units, still reconstruction-internal.

**Layer A width does survive blocking as an association, but not as a Schellart test.** ρ=0.216, p=2.57e-23; time-dummy LRT p=2.10e-33. The sign is **opposite** Schellart: wide T3 P(bab_zone)=0.173 vs narrow T1 0.045; MH OR (narrow vs wide) = 0.204 (p=9.69e-15). Wider reconstructed zones are more likely to host a back-arc MOR in this model.

**Layer A edge association also survives blocking, also opposite Schellart** (ρ=0.141, time-dummy p=2.45e-07; far tertile P(bab_zone) higher than near). The Schellart 1200 km cut is **degenerate at zone-median** (n_far=29 of 2077); almost every zone has median edge distance < 1200 km, so that split is not an interior-vs-edge test.

**Layer B width/edge do not survive time dummies** (width LRT p=0.1189; edge LRT p=0.3616). Continuous median dilatation vs width among network blocks is null (ρ=0.007, p=0.8057).

0 Ma SubMap C vs E (independent of bab_on): width ρ=-0.189 (p=0.1035, n=75), matching the weak 0 Ma width result in `cenozoic_age_gate.md`.

No hypothesis winner. The Schellart *direction* does not survive this blocking.

## Figures

- `/workspace/subduction-test/results/fig_width_vs_bab.png`
- `/workspace/subduction-test/results/fig_edge_vs_bab.png`
- `/workspace/subduction-test/results/fig_dilatation_vs_width_edge.png`
- `/workspace/subduction-test/results/fig_width_edge_map_0Ma.png` (0 Ma map; teal rings = `bab_on`)

## How to rerun

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \
    /workspace/subduction-test/scripts/schellart_width_edge.py
```

## Caveats

- Reconstruction-internal. Opening histories of back-arc basins are encoded
  in the topologies that define `bab_on`.
- `bab_zone=1` if *any* point in the zone is `bab_on`, so a long connected
  trench with a local MOR (e.g. Mariana Trough on an Izu–Bonin–Mariana
  polyline) flags the whole zone.
- Median edge distance at zone level is not a within-zone interior-vs-edge
  contrast; it covaries with width.
- Hellenic / Okinawa-style continental extension is layer B, not A.
- Absolute frame / H1 vs H2 is not tested here. See
  `results/submap_h1h2_frames.md` for SubMap Vupn/Vtn across frames.
