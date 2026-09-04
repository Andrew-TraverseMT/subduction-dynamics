# SubMap H1/H2 across absolute frames

Sidecar to `results/submap_m56_h1h2.md`. Values are computed from
`/workspace/subduction-test/data/submap/submap_transects_full.csv`. Nothing is invented. UPS recode is the same function as
`scripts/submap_m56_h1h2.py`.

- **n rows in full sheet:** 260
- **Units:** SubMap stores mm/yr; medians below are **cm/yr** (÷ 10). Spearman
  is rank-based (unit-free). Sentinel −999 → NaN.
- **This sidecar tests absolute V_UP vs V_T only.** `vdn` vs UPS is circular
  and is not repeated. Müller `v_T` vs `v_OP` is not used.

## Frame names (do not invent)

Cerpa, Lallemand & Heuret (2025, §II.1) calculated absolute velocities in
**four** M56 reference frames, listed as:

1. MORVEL56-NNR (Argus et al., 2011)
2. spreading-alignment (Becker et al., 2015)
3. T25m (Wang et al., 2018)
4. GMHRF-1Ma (Doubrovine et al., 2012)

The SubMap dump labels columns `Vupn1`…`Vupn4` / `Vtn1`…`Vtn4` without a
name field. The existing analysis documents **frame 1 = MORVEL56-NNR**.
Cerpa does **not** map frames 2–4 onto column numbers, so those three stay
**Frame N (unknown)**. They are the remaining three names in the list above,
in unknown order.

Columns `Vupn5`…`Vupn9` / `Vtn5`…`Vtn9` sit next to the `N15_*` relative
velocities (NUVEL-1A-era kinematics in this dump). No N15 absolute-frame
names appear in `data/submap/README.md`, Cerpa fulltext, or the JSON sheets.
Frames 5–9 are **Frame N (unknown)**. HS3 / NNR / hotspot labels are **not**
assigned.

| n | family | columns | name | name known? | n finite both | n C vs E |
|---:|---|---|---|---|---:|---:|
| 1 | M56 | `Vupn1` / `Vtn1` | MORVEL56-NNR (Argus et al., 2011) | yes | 245 | 79 |
| 2 | M56 | `Vupn2` / `Vtn2` | Frame 2 (unknown) | no — unknown | 245 | 79 |
| 3 | M56 | `Vupn3` / `Vtn3` | Frame 3 (unknown) | no — unknown | 245 | 79 |
| 4 | M56 | `Vupn4` / `Vtn4` | Frame 4 (unknown) | no — unknown | 245 | 79 |
| 5 | N15 | `Vupn5` / `Vtn5` | Frame 5 (unknown) | no — unknown | 249 | 77 |
| 6 | N15 | `Vupn6` / `Vtn6` | Frame 6 (unknown) | no — unknown | 179 | 47 |
| 7 | N15 | `Vupn7` / `Vtn7` | Frame 7 (unknown) | no — unknown | 179 | 47 |
| 8 | N15 | `Vupn8` / `Vtn8` | Frame 8 (unknown) | no — unknown | 174 | 46 |
| 9 | N15 | `Vupn9` / `Vtn9` | Frame 9 (unknown) | no — unknown | 179 | 47 |

## Sample and recode

C = C + C-SS; E = E + E-SS; N = Neutral; SS = SS + SS-C + SS-E.
Undefined dropped. SS held out of C vs E. Score C=+1, E=−1 (shortening-positive,
matching SubMap `vdn`). A frame is used if C vs E has n≥10 and finite ρ for
both Vupn and Vtn.

Sign of absolute trench-normal components (from the frame-1 identity
`vdn ≈ Vupn1 − Vtn1` in the parent analysis): **positive Vupn / Vtn = oceanward**.

Frame 1 (MORVEL56-NNR) C vs E: ρ(Vupn)=0.498 (p=2.99e-06), ρ(Vtn)=-0.341 (p=0.0021), n=79, ρ(Vupn,Vtn)=0.389.

## Spearman C vs E, collinearity, ranking

| frame | family | n C/E | ρ(Vupn, C/E) | p | ρ(Vtn, C/E) | p | ρ(Vupn, Vtn) | p | |ρ_H1|−|ρ_H2| | stronger |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 MORVEL56-NNR (Argus et al., 2011) | M56 | 79 | 0.498 | 2.99e-06 | -0.341 | 0.0021 | 0.389 | 0.0004 | 0.157 | H1 Vupn |
| 2 Frame 2 (unknown) | M56 | 79 | 0.475 | 9.79e-06 | -0.501 | 2.51e-06 | 0.253 | 0.0247 | -0.027 | H2 Vtn |
| 3 Frame 3 (unknown) | M56 | 79 | 0.489 | 4.74e-06 | -0.488 | 5.04e-06 | 0.250 | 0.0263 | 0.001 | H1 Vupn |
| 4 Frame 4 (unknown) | M56 | 79 | 0.250 | 0.0263 | -0.621 | 1.04e-09 | 0.373 | 0.0007 | -0.371 | H2 Vtn |
| 5 Frame 5 (unknown) | N15 | 77 | 0.502 | 3.29e-06 | -0.032 | 0.7814 | 0.596 | 1.08e-08 | 0.470 | H1 Vupn |
| 6 Frame 6 (unknown) | N15 | 47 | 0.508 | 0.0003 | -0.140 | 0.3464 | 0.345 | 0.0175 | 0.368 | H1 Vupn |
| 7 Frame 7 (unknown) | N15 | 47 | 0.609 | 5.45e-06 | -0.108 | 0.4707 | 0.281 | 0.0553 | 0.502 | H1 Vupn |
| 8 Frame 8 (unknown) | N15 | 46 | 0.611 | 6.48e-06 | -0.100 | 0.5079 | 0.295 | 0.0466 | 0.511 | H1 Vupn |
| 9 Frame 9 (unknown) | N15 | 47 | 0.251 | 0.0882 | -0.256 | 0.0819 | 0.481 | 0.0006 | -0.005 | H2 Vtn |

ρ(Vupn, Vtn) on the C vs E rows is the collinearity check. Müller2019
`v_T ≈ −v_OP` is **not** this number; SubMap can separate the two series
when ρ is well below 1.

Sign of ρ(Vupn, C/E) stable across usable frames: yes. Sign of ρ(Vtn, C/E) stable: yes. Predicted signs (SubMap convention, + = oceanward): H1 ρ(Vupn)>0 (OP advance with C); H2 ρ(Vtn)<0 (trench retreat / rollback with E).

### Does H1 vs H2 ranking flip?

Ranking = which of |ρ(Vupn, C/E)| and |ρ(Vtn, C/E)| is larger.

H1 vs H2 **ranking flips** across frames: the larger-|ρ| predictor is frame 1 → H1 Vupn, frame 2 → H2 Vtn, frame 3 → H1 Vupn, frame 4 → H2 Vtn, frame 5 → H1 Vupn, frame 6 → H1 Vupn, frame 7 → H1 Vupn, frame 8 → H1 Vupn, frame 9 → H2 Vtn. Frames 3, 9 are near-ties (|Δρ|<0.02) and should not be read as a ranking. Unambiguous flips remain (e.g. frame 4 vs frame 1). That is the Schellart (2008) point: absolute-frame choice can reorder V_UP vs V_T. No hypothesis winner.

C/N/E ordinal (C=+1, N=0, E=−1) and collinearity on all defined rows:

| frame | ρ(Vupn, C/N/E) | p | n | ρ(Vtn, C/N/E) | p | n | ρ(Vupn, Vtn) all defined | p | n |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 M56 | 0.460 | 2.78e-06 | 95 | -0.299 | 0.0033 | 95 | 0.655 | 1.79e-29 | 229 |
| 2 M56 | 0.441 | 7.54e-06 | 95 | -0.456 | 3.49e-06 | 95 | 0.589 | 9.37e-23 | 229 |
| 3 M56 | 0.454 | 3.73e-06 | 95 | -0.447 | 5.66e-06 | 95 | 0.581 | 4.47e-22 | 229 |
| 4 M56 | 0.227 | 0.0266 | 95 | -0.565 | 2.44e-09 | 95 | 0.662 | 2.91e-30 | 229 |
| 5 N15 | 0.464 | 2.76e-06 | 93 | -0.016 | 0.8800 | 93 | 0.775 | 4.57e-45 | 219 |
| 6 N15 | 0.405 | 0.0012 | 61 | -0.143 | 0.2704 | 61 | 0.705 | 5.42e-26 | 164 |
| 7 N15 | 0.492 | 5.73e-05 | 61 | -0.111 | 0.3946 | 61 | 0.677 | 2.40e-23 | 164 |
| 8 N15 | 0.488 | 7.53e-05 | 60 | -0.100 | 0.4456 | 60 | 0.654 | 8.98e-21 | 159 |
| 9 N15 | 0.173 | 0.1835 | 61 | -0.268 | 0.0368 | 61 | 0.758 | 6.28e-32 | 164 |

## Median Vupn / Vtn by UPS_3 (cm/yr)

**Frame 1** — MORVEL56-NNR (Argus et al., 2011) (M56; n_defined=229)

| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |
|---|---:|---|---|
| C | 40 | 0.900 cm/yr (n=40) | -0.250 cm/yr (n=40) |
| N | 16 | -0.600 cm/yr (n=16) | -0.600 cm/yr (n=16) |
| E | 39 | -0.900 cm/yr (n=39) | 2.200 cm/yr (n=39) |
| SS | 134 | 0.650 cm/yr (n=134) | 0.800 cm/yr (n=134) |

**Frame 2** — Frame 2 (unknown) (M56; n_defined=229)

| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |
|---|---:|---|---|
| C | 40 | 1.200 cm/yr (n=40) | -0.150 cm/yr (n=40) |
| N | 16 | 0.100 cm/yr (n=16) | 0.100 cm/yr (n=16) |
| E | 39 | -0.400 cm/yr (n=39) | 2.700 cm/yr (n=39) |
| SS | 134 | 0.700 cm/yr (n=134) | 0.600 cm/yr (n=134) |

**Frame 3** — Frame 3 (unknown) (M56; n_defined=229)

| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |
|---|---:|---|---|
| C | 40 | 1.100 cm/yr (n=40) | -0.200 cm/yr (n=40) |
| N | 16 | 0.100 cm/yr (n=16) | 0.100 cm/yr (n=16) |
| E | 39 | -0.400 cm/yr (n=39) | 2.700 cm/yr (n=39) |
| SS | 134 | 0.700 cm/yr (n=134) | 0.600 cm/yr (n=134) |

**Frame 4** — Frame 4 (unknown) (M56; n_defined=229)

| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |
|---|---:|---|---|
| C | 40 | 0.600 cm/yr (n=40) | -0.450 cm/yr (n=40) |
| N | 16 | 0.400 cm/yr (n=16) | 0.400 cm/yr (n=16) |
| E | 39 | 0.000 cm/yr (n=39) | 2.700 cm/yr (n=39) |
| SS | 134 | 0.950 cm/yr (n=134) | 0.650 cm/yr (n=134) |

**Frame 5** — Frame 5 (unknown) (N15; n_defined=219)

| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |
|---|---:|---|---|
| C | 38 | 0.600 cm/yr (n=38) | -0.150 cm/yr (n=38) |
| N | 16 | -2.500 cm/yr (n=16) | -2.500 cm/yr (n=16) |
| E | 39 | -2.700 cm/yr (n=39) | 0.700 cm/yr (n=39) |
| SS | 126 | -0.400 cm/yr (n=126) | -0.100 cm/yr (n=126) |

**Frame 6** — Frame 6 (unknown) (N15; n_defined=164)

| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |
|---|---:|---|---|
| C | 17 | 1.100 cm/yr (n=17) | 0.800 cm/yr (n=17) |
| N | 14 | -0.200 cm/yr (n=14) | -0.200 cm/yr (n=14) |
| E | 30 | -0.500 cm/yr (n=30) | 1.750 cm/yr (n=30) |
| SS | 103 | 1.100 cm/yr (n=103) | 1.000 cm/yr (n=103) |

**Frame 7** — Frame 7 (unknown) (N15; n_defined=164)

| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |
|---|---:|---|---|
| C | 17 | 1.300 cm/yr (n=17) | 1.000 cm/yr (n=17) |
| N | 14 | -0.600 cm/yr (n=14) | -0.600 cm/yr (n=14) |
| E | 30 | -1.100 cm/yr (n=30) | 1.450 cm/yr (n=30) |
| SS | 103 | 1.200 cm/yr (n=103) | 1.000 cm/yr (n=103) |

**Frame 8** — Frame 8 (unknown) (N15; n_defined=159)

| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |
|---|---:|---|---|
| C | 17 | 1.200 cm/yr (n=17) | 1.000 cm/yr (n=17) |
| N | 14 | -0.950 cm/yr (n=14) | -0.950 cm/yr (n=14) |
| E | 29 | -1.500 cm/yr (n=29) | 1.100 cm/yr (n=29) |
| SS | 99 | 1.000 cm/yr (n=99) | 1.000 cm/yr (n=99) |

**Frame 9** — Frame 9 (unknown) (N15; n_defined=164)

| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |
|---|---:|---|---|
| C | 17 | 0.700 cm/yr (n=17) | 0.500 cm/yr (n=17) |
| N | 14 | 0.650 cm/yr (n=14) | 0.650 cm/yr (n=14) |
| E | 30 | 0.200 cm/yr (n=30) | 2.300 cm/yr (n=30) |
| SS | 103 | 0.900 cm/yr (n=103) | 0.800 cm/yr (n=103) |


## Figure

- `/workspace/subduction-test/results/fig_h1h2_by_frame.png` — forest / dot plot of ρ(Vupn) and ρ(Vtn) vs C vs E.
  Whiskers are approximate 95% Fisher-z intervals, not a bootstrap.

## How to rerun

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \
    /workspace/subduction-test/scripts/submap_h1h2_frames.py
```

## Caveats

- Geodetic UPS is classified from `Ovd` / `vd`. This sidecar does not treat
  `vdn` as a predictor.
- Absolute velocities change with the reference frame (Schellart 2008). A
  result that exists in only one named frame is fragile; unnamed frames
  cannot be interpreted as a specific hotspot or NNR model.
- N15 and M56 relative plate pairs differ; n finite is not the same in every
  column (some N15 frames have more −999).
- Pure C and pure E are small before recode; C-SS and E-SS dominate the
  C vs E n.
- Fisher-z intervals treat Spearman ρ like a Pearson correlation of ranks;
  they are a display aid.
