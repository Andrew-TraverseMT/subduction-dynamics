# Geological UPS test — Cenozoic, independent of Müller back-arc MORs

Built by `scripts/geological_ups_test.py`.

**Outcome:** published geological / marine-geophysical UPS intervals
(`data/geological_ups/cenozoic_ups_intervals.csv`).
**Predictors:** Müller et al. (2019) trench kinematics already in
`results/trench_kinematics_0-60Ma_Muller2019.csv` (seafloor age, trench-zone
width, distance to nearest trench edge). Trench kinematics were **not**
re-extracted.

This is **not** `bab_on` and **not** deforming-network dilatation. Those
proxies are reconstruction-internal (the topologies already know where they
opened back-arc basins). This test asks whether an *independent* geological
UPS series still shows a Sdrolias-style age gate or a Schellart-style
narrow / near-edge preference.

## Sparse by construction — read this first

- **n intervals:** 22 (E=16, C=2, N=4)
- **n named segments:** 14
- **n trench sample-times:** 76529 (61 times, 0–60 Ma, 1 Myr)
- **n classified sample-times:** 8094 (10.6%)
- **n Unknown:** 68435
- **n segment × time blocks (classified):** 295
- sample-level classified: E=4403, C=2270, N=1421
- block-level: E=177, C=56, N=62

Hand intervals cover well-dated back-arc basins plus Central Andean
shortening and Cascadia (no back-arc). Most of the global trench catalog
is `Unknown` and is **not** treated as C. Ages have published ranges;
midpoints are used. Where literature conflicts, the CSV `notes` field
records the conflict. **Do not read a p-value on 76k points as independent
information** — neighbouring tessellated samples along a trench are not
independent. The inference unit is **named segment × time**.

Cenozoic C is dominated by the Central Andes box (one wide, long-lived
shortening orogen). Cenozoic E is several short-lived, mostly western-Pacific
and Mediterranean basins. That imbalance is geological, not a sampling trick,
but it means E vs C is largely “Pacific BABs vs Andes”.

## How samples are classified

A trench sample is classified if its longitude/latitude at that reconstruction
time falls inside an active interval bbox (`t_end_Ma ≤ time ≤ t_start_Ma`).
Longitude wrap (`lon_min > lon_max`) is used for Tonga and Kermadec.
If two intervals hit the same sample, **E then C then N**. Otherwise Unknown.
Bboxes follow the trench in this catalog, not the basin floor; false-positive
risk is documented in `data/geological_ups/README.md`.

Aegean, Okinawa, Havre, and Tyrrhenian count as **E during rifting**, even
when there is no MOR. Cascadia is **N** (no back-arc; margin-parallel
shortening in Washington is not Andean C).

## Age — Sdrolias gate

Sdrolias & Müller (2006) associated back-arc opening with subducting
lithosphere older than ~55 Ma. Test among classified rows
with finite `seafloor_age_Ma`. Hellenic / Calabria ages in this catalog are
often Mesozoic (Tethyan) or NaN; they still enter E when classified.
Dropping those two segments (sensitivity below) makes E **younger**, not older.

The Cenozoic geological catalog does **not** reproduce the present-day
old-slab → back-arc association. West Philippine / early IBM extension
(55–30 Ma) sits on ~17–40 Ma Pacific lithosphere in these age grids,
while Central Andean shortening sits on older Nazca. Gate OR is < 1
(younger slabs more often E among classified blocks).

### Segment × time (preferred)

- Median age **E vs C**: 61.1 Ma (n=177) vs 80.4 Ma (n=56); Mann–Whitney p=0.0401
- Median age **E vs not-E** (C+N): 61.1 Ma (n=177) vs 80.0 Ma (n=118); p=0.3973
- Gate among classified blocks: P(E | age≥55) = 0.544 (98/180) vs P(E | age<55) = 0.687 (79/115); odds ratio = 0.5; Fisher p=0.0153
- Gate **E vs C only** (drop N): P(E | age≥55) = 0.658 (98/149) vs P(E | age<55) = 0.940 (79/84); OR = 0.1; p=3.74e-07
- Sensitivity **drop Hellenic+Calabria** (Tethyan ages): median age E vs C 51.9 vs 80.4 Ma (n=142/56); P(E|≥55)=0.434 (63/145) vs P(E|<55)=0.687 (79/115); OR=0.4; p=5.72e-05

### Sample-level (autocorrelated; supplementary)

- Median age **E vs C**: 65.3 Ma (n=3856) vs 80.8 Ma (n=2267); p=7.13e-24
- Gate: P(E | age≥55) = 0.442 (2168/4907) vs P(E | age<55) = 0.683 (1688/2470); OR = 0.4; p=5.85e-87
- Gate E vs C only: OR = 0.2; p=8.11e-168 (E|old 2168/4177; E|young 1688/1946)

Sample-level p-values will look decisive because n is large and points along
a trench are copies of each other. Believe the block-level n.

## Width and edge — Schellart direction

Schellart et al. (2007, *Nature*): rollback / back-arc opening concentrate
on **narrow** slabs and **near a lateral slab edge**; interiors of wide slabs
are more compressive. Operationally, at segment × time:

- narrower `trench_zone_width_km` → more E → Spearman ρ(width, is_E) **negative**
- smaller `dist_nearest_trench_edge_km` → more E → ρ(edge, is_E) **negative**
- near-edge split at 1200 km

`trench_zone_width_km` is the along-trench length of a connected polyline
sharing a subducting-plate ID in Müller2019, **not** a tomographic slab width.
Median edge distance is mechanically related to width (a uniform sample on a
zone of length W has median distance-to-edge ≈ W/4).

- Median width **E vs C**: 6981.2 km (n=177) vs 10270.4 km (n=56); p=1.52e-05
- Median edge **E vs C**: 500.4 km (n=177) vs 1501.1 km (n=56); p=4.55e-10
- Spearman ρ(width, is_E) = -0.039 (p=0.5015, n=295) — Schellart sign is negative
- Spearman ρ(edge, is_E) = -0.156 (p=0.0072, n=295) — Schellart sign is negative
- P(E | median edge < 1200 km) = 0.673 (140/208) vs P(E | ≥1200) = 0.425 (37/87)

Andes is a wide C orogen; Scotia / Manus / Andaman are short E segments.
A width contrast in the E-vs-C medians is **consistent with** Schellart
but is also the geometry of boxing those systems. Spearman ρ(width, is_E)
across *all* classified blocks (including N and IBM) is near zero because
Müller `trench_zone_width_km` often glues IBM–Japan–Kuril–Aleutian into a
10,000–20,000 km polyline — that is **not** a tomographic slab width and
is not the Schellart 2007 width of the Mariana or Tonga slab. Scotia
(~800 km) and Manus (~1,100 km) are the honest narrow E end-member.
Do not over-claim.

## 0 Ma sanity vs SubMap geodetic UPS

Geological UPS at 0 Ma assigned to SubMap transects using Müller coordinates in `joined_submap_muller2019_0Ma.csv`. Overlap classified: **86** of 241 joined transects; **83** with defined `UPS_3` (not Undefined).

Agreement:
- **strict** (E↔E, C↔C, N↔N): 37/83 = 0.45
- **relaxed** (E↔E, C↔C, N↔{N,SS}): 47/83 = 0.57

Relaxed N↔SS is the intended Cascadia rule (Wells et al. 1998;
SubMap Cascadia is SS, we labelled N). Andaman is an expected
**disagreement**: geological E (oblique Central Andaman spreading) vs
SubMap SS.

Crosstab geological UPS × SubMap `UPS_3` (defined only):

| geo\\SubMap | E | C | N | SS |
|---|---:|---:|---:|---:|
| E | 31 | 3 | 1 | 21 |
| C | 0 | 6 | 0 | 9 |
| N | 2 | 0 | 0 | 10 |

| segment | n | geo | SubMap UPS_3 | n relaxed-agree |
|---|---:|---|---|---:|
| Andaman | 5 | E | {'SS': 5} | 0 |
| Andes | 11 | C | {'SS': 6, 'C': 5} | 5 |
| Calabria | 2 | E | {'SS': 1, 'E': 1} | 1 |
| Cascadia | 5 | N | {'SS': 5} | 5 |
| Hellenic | 8 | E | {'E': 6, 'SS': 2} | 6 |
| Izu-Bonin | 7 | N | {'SS': 5, 'E': 2} | 5 |
| Japan | 4 | C | {'SS': 3, 'C': 1} | 1 |
| Kermadec | 7 | E | {'E': 6, 'SS': 1} | 6 |
| Manus | 4 | E | {'E': 2, 'SS': 1, 'C': 1} | 2 |
| Mariana | 7 | E | {'E': 4, 'SS': 3} | 4 |
| NewHebrides | 4 | E | {'SS': 2, 'C': 1, 'E': 1} | 1 |
| Ryukyu | 7 | E | {'E': 4, 'SS': 2, 'C': 1} | 4 |
| Scotia | 6 | E | {'SS': 3, 'E': 2, 'N': 1} | 2 |
| Tonga | 6 | E | {'E': 5, 'SS': 1} | 5 |

## Intervals used

Full citations: `data/geological_ups/README.md`.

| zone | segment | class | start Ma | end Ma | age basis (short) |
|---|---|---|---:|---:|---|
| Central_Andaman_Basin | Andaman | E | 4 | 0 | Kamesh Raju/Curray 4 Ma spreading; Jha et al. 2021 revise to 5.9 Ma; used 4.0 (range 4-6) |
| Central_Andes_shortening | Andes | C | 50 | 0 | shortening onset 50-45 Ma (range 60-40); continuous Cenozoic shortening to present; used 50-0 |
| Tyrrhenian_Calabria | Calabria | E | 10 | 0 | major rifting from ~10 Ma (Tortonian); Vavilov spreading ~4 Ma; Marsili ~2 Ma |
| Cascadia_no_backarc | Cascadia | N | 30 | 0 | Juan de Fuca-Cascadia configuration in this catalog from ~30 Ma; no Cenozoic back-arc spreading |
| Aegean_rifting | Hellenic | E | 23 | 0 | Ring et al. 2010 Aegean-wide extension from ~23 Ma; midpoint of 15 vs 35-45 Ma debate; used 23-0 |
| West_Philippine_Basin_IzuBonin | Izu-Bonin | E | 55 | 30 | same WPB spreading interval as Mariana row; northern proto-IBM may be sparsely sampled |
| Shikoku_Basin | Izu-Bonin | E | 29 | 15 | midpoint of 30 Ma rifting vs 27-26 Ma organized spreading; cessation 15 Ma |
| Izu_Bonin_post_Shikoku | Izu-Bonin | N | 15 | 0 | Shikoku cessation ~15 Ma; no subsequent back-arc spreading behind Izu-Bonin |
| Japan_Sea_opening | Japan | E | 28 | 14 | ODP 127/128 major opening 28-18 Ma; rifting from ~32 Ma; termination ~15-12 Ma; used 28-14 |
| Japan_post_opening_neutral | Japan | N | 14 | 4 | Sato 1994: after ~14 Ma termination a neutral stress regime until ~4 Ma |
| NE_Japan_inversion | Japan | C | 4 | 0 | Sato 1994: strong compression from ~4 Ma; inversion of Miocene normal faults |
| Havre_Trough | Kermadec | E | 5.5 | 0 | initial splitting ~5.5 Ma; southern Havre may be <2 Ma (conflict noted) |
| Manus_Basin | Manus | E | 3.5 | 0 | Manus Spreading Centre <3.5 Ma (Taylor 1979); East Manus rifting <0.78 Ma |
| West_Philippine_Basin_Mariana | Mariana | E | 55 | 30 | midpoint of published spreading range; rifting ~55 Ma, spreading ~54 to 33/30 Ma (used 55-30) |
| Parece_Vela_Basin | Mariana | E | 30 | 15 | initiation ~30 Ma; organized spreading by ~29-28 Ma (chron 10-9); cessation ~15 Ma |
| IBM_post_PV_gap | Mariana | N | 15 | 7 | bracketed by PV cessation (~15 Ma) and Mariana Trough rifting (~7 Ma); not a dated 'neutral' event |
| Mariana_Trough | Mariana | E | 7 | 0 | midpoint of rifting 8-6 Ma and spreading from ~6-5 Ma; used 7-0 |
| North_Fiji_Basin | NewHebrides | E | 12 | 0 | rifting/opening from ~12 Ma after Vitiaz polarity reversal; still spreading |
| Okinawa_Trough_phase1 | Ryukyu | E | 10 | 6 | Miki 1995 first phase 10-6 Ma; Sibuet et al. 1987 middle-late Miocene first rifting |
| Okinawa_Trough_phase2 | Ryukyu | E | 2 | 0 | second rifting phase from Plio-Pleistocene boundary ~2 Ma |
| East_Scotia_Ridge | Scotia | E | 16 | 0 | midpoint of 15-17 Ma spreading underway (Larter) vs possible C6 ~20 Ma start; used 16-0 |
| Lau_Basin | Tonga | E | 6 | 0 | opening ~6 Ma; spreading 6-5.5 Ma in the north, southward propagation |

## Segment × time coverage

| segment | n times | t range Ma | classes | nE | nC | nN | med age | med width km |
|---|---:|---|---|---:|---:|---:|---:|---:|
| Andaman | 5 | 4–0 | E | 5 | 0 | 0 | 79.7 | 2552.4 |
| Andes | 51 | 50–0 | C | 0 | 51 | 0 | 79.3 | 11304.6 |
| Calabria | 11 | 10–0 | E | 11 | 0 | 0 | 224.9 | 6848.8 |
| Cascadia | 31 | 30–0 | N | 0 | 0 | 31 | 8.2 | 1739.9 |
| Hellenic | 24 | 23–0 | E | 24 | 0 | 0 | 231.0 | 6281.4 |
| Izu-Bonin | 41 | 40–0 | E,N | 26 | 0 | 15 | 71.0 | 10308.5 |
| Japan | 21 | 20–0 | C,E,N | 7 | 5 | 9 | 94.6 | 10188.4 |
| Kermadec | 6 | 5–0 | E | 6 | 0 | 0 | 114.5 | 4114.1 |
| Manus | 4 | 3–0 | E | 4 | 0 | 0 | 27.9 | 1089.9 |
| Mariana | 56 | 55–0 | E,N | 49 | 0 | 7 | 52.7 | 10325.0 |
| NewHebrides | 13 | 12–0 | E | 13 | 0 | 0 | 38.8 | 2399.0 |
| Ryukyu | 8 | 10–0 | E | 8 | 0 | 0 | 61.9 | 2084.2 |
| Scotia | 17 | 16–0 | E | 17 | 0 | 0 | 46.4 | 803.1 |
| Tonga | 7 | 6–0 | E | 7 | 0 | 0 | 94.0 | 4427.6 |

## Reading the result (no winner claim unless n is honest)

This is a **sparse hand catalog** (22 intervals; 8094 classified sample-times; 295 segment×time blocks). Ages have published ranges. At segment×time the 55 Ma gate runs the **opposite** way of Sdrolias & Müller (2006): P(E | age≥55) < P(E | age<55) (OR<1). Median age of E is younger than C. A large part of that is West Philippine / early IBM geological E sitting on young (≲40 Ma) Pacific lithosphere in this age grid, plus Andes C on older Nazca. Do not declare a Sdrolias-gate winner; if anything the independent catalog does not reproduce the present-day old-slab → back-arc association through the Cenozoic. Width / edge do **not** give a strong, well-powered Schellart confirmation at this n. No hypothesis winner. 0 Ma relaxed agreement with SubMap geodetic UPS is 0.57 (47/83). Cascadia N↔SS is counted as agree; Andaman E vs SS is a known mismatch (oblique spreading). **No hypothesis is declared the winner.** The point of this file is an independent UPS series plus honest n, not a knockout of Heuret–Lallemand vs Schellart.

## Outputs

- `results/geological_ups_joined_summary.csv` — segment × time summary
- `results/fig_geo_ups_age.png`
- `results/fig_geo_ups_width_edge.png`
- `results/fig_geo_ups_timeline.png`
- `data/geological_ups/cenozoic_ups_intervals.csv`
- `data/geological_ups/README.md`

## Cite

- Geological ages: `data/geological_ups/README.md` (primary papers per basin).
- Müller, R. D., et al. (2019). *Tectonics*. https://doi.org/10.1029/2018TC005462
- Sdrolias, M., & Müller, R. D. (2006). *G3*. https://doi.org/10.1029/2005GC001090
- Schellart, W. P., et al. (2007). *Nature*. https://doi.org/10.1038/nature05737
- Schellart, W. P. (2024). *Earth-Sci. Rev.* https://doi.org/10.1016/j.earscirev.2024.104755
- SubMap: Cerpa, Lallemand & Heuret (2025). https://doi.org/10.31223/x56m8h

