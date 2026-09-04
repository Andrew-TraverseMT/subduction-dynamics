# Competing controls on upper-plate stress at subduction zones
### A reconstruction-based test of Heuret–Lallemand, Schellart, and Sdrolias–Müller

**Grok, Tectonics Researcher** · Advisor: Andrew Laskowski  
Repo: [Andrew-TraverseMT/subduction-dynamics](https://github.com/Andrew-TraverseMT/subduction-dynamics)

---

<!-- Slide: Title context -->
## One-sentence takeaway

Present-day correlations between subduction drivers and upper-plate stress look tidy. Once we add **time**, **spatial blocking**, **absolute-frame choice**, and an **independent geological UPS**, those tidy stories mostly fall apart — including the famous ~55 Ma age gate, which runs the **other way** in this Cenozoic catalog.

---

## Outline

1. Hypothesis and why it matters  
2. Background: three competing stories  
3. Methods: tools, data, tests  
4. Results: 0 Ma, Cenozoic proxy, frames, geological UPS  
5. Interpretation and what is still open  

---

## The problem

Since Uyeda and Kanamori (1979), the first-order puzzle is why some overriding plates are **Chilean** (compressive) and others **Mariana** (extensional).

Drivers people keep nominating:

- oceanic crust age / slab pull  
- convergence rate and obliquity  
- upper-plate absolute motion  
- trench rollback / slab width / distance to a slab edge  

Those predictors co-vary geographically today (old, steep, retreating western Pacific vs young, flat, compressive eastern Pacific). Pairwise 0 Ma scatterplots cannot isolate a driver.

---

## Hypotheses we tested

| ID | Story | Prediction |
|---|---|---|
| **H1** | Heuret and Lallemand (2005) | Upper-plate absolute motion leads: OP retreat → extension; advance → compression |
| **H2** | Schellart et al. (2007, 2008, 2024) | Trench rollback / slab width leads: narrow slabs and near-edge segments extend; wide interiors compress |
| **H3** | Molnar and Atwater (1978); Sdrolias and Müller (2006) | Old subducting lithosphere (≥ ~55 Ma) gates back-arc basins |
| **H4** | Convergence / obliquity | Faster trench-normal coupling favors compression (often fails in the literature) |
| **H5** | Interactions | No univariate winner; age × kinematics or width × rollback |

---

## What is new here

Not a new cartoon of subduction.

A **state-of-the-art empirical horse-race**:

1. Dense trench kinematics from **pyGPlates / GPlately** (Müller et al. 2019)  
2. Independent present-day UPS from **SubMap v7** geodetic classes (Cerpa, Lallemand and Heuret 2025)  
3. Time: **0–60 Ma** at 1 Myr steps  
4. Multiple absolute frames for H1 vs H2  
5. Honest sample units: **zone × time** blocks, not 76k autocorrelated trench points  

---

## Methods — kinematic engine

- **Model:** Müller et al. (2019) deforming plate model  
- **Tools:** pyGPlates 1.0, GPlately 2.0  
- **Sampling:** ~0.5° (~56 km) along every trench, every 1 Myr from 0–60 Ma → **76,529** samples  
- **Predictors kept:** seafloor age, trench-zone width, distance to trench edge, convergence  
- **Not usable as independent H1/H2 in this model:** \(v_{T\perp} \approx -v_{\mathrm{OP}\perp}\) (trench glued to overriding plate; ρ ≈ −0.999)

---

## Methods — present-day UPS (independent)

**SubMap v7** geodetic UPS (260 transects):

- Classes: C, C-SS, SS-C, SS, SS-E, E-SS, E, Neutral, Undefined  
- Recode for tests: **C** = C+C-SS, **E** = E+E-SS, SS held separate  
- Absolute velocities in up to 9 frames (`Vupn1`…`Vtn9`); frame 1 = MORVEL56-NNR  

Joined 241 / 260 to Müller trenches (median distance **25 km**).

---

## Methods — Cenozoic UPS outcomes (two generations)

**Generation 1 (reconstruction-internal, circular by design):**

- Layer A: landward back-arc MOR on the overriding plate (`bab_on`)  
- Layer B: deforming-network dilatation  

Useful as a Sdrolias-style kinematic association. **Not** independent geology.

**Generation 2 (independent geological intervals):**

- **22** published intervals on **14** named segments (16 E, 2 C, 4 N)  
- Opening/cessation ages from the literature (WPB, Parece Vela, Shikoku, Mariana Trough, Lau/Havre, Scotia, Andaman, Japan Sea, Okinawa, Aegean, Tyrrhenian, Manus, North Fiji; Central Andes shortening; Cascadia N)  
- Joined by geographic bbox to trench samples through time  
- Only **10.6%** of sample-times classified; rest Unknown (sparse by design)

---

## Methods — statistics (honesty rules)

- Inference unit for Cenozoic = **zone × time** (n ≈ 2077), not points  
- Univariate Spearman + logistic with **time dummies**  
- Mantel–Haenszel / tertile splits for gates  
- No winner claimed unless ranking is strong **and** stable across frames / blocking  
- Circularity stated whenever the outcome shares DNA with the plate model  

---

## Results — 0 Ma SubMap H1 vs H2 (frame 1)

C vs E (n = 79), MORVEL56-NNR:

| Predictor | Spearman ρ | p |
|---|---:|---:|
| \(V_{\mathrm{UPn}}\) (H1) | **+0.50** | 3×10⁻⁶ |
| \(V_{\mathrm{Tn}}\) (H2) | **−0.34** | 0.002 |
| Seafloor age | **−0.50** | 0.0002 |

Compression goes with oceanward upper plate. Extension goes with oceanward trench. Overlap is large. Age also separates C (median ~38 Ma) from E (~94 Ma) at 0 Ma.

---

## Results — absolute frame choice rearranges the horse-race

Signs are stable (H1 ρ > 0, H2 ρ < 0).

**Ranking is not:**

- Frame 1: H1 stronger (|ρ| 0.50 vs 0.34)  
- Frame 4: H2 stronger (0.62 vs 0.25)  

That is Schellart’s (2008) point, live: absolute-frame choice can manufacture a winner.

---

## Results — Cenozoic age gate (reconstruction-internal)

At **0 Ma**, old slabs look special: median age 97 Ma (`bab_on`) vs 44 Ma (`bab_off`); P(on\|≥55) = 0.25 vs 0.02.

After **zone × time** blocking:

- Mantel–Haenszel OR for the 55 Ma gate ≈ **1.21**, p ≈ 0.24  
- Age logistic after time dummies: **null** (p ≈ 0.43)  
- Gate **flips** in several intervals  

**The 0 Ma age headline does not survive Cenozoic zone-blocking.**

---

## Results — Schellart width / edge (reconstruction-internal)

Schellart prediction: **narrow** and **near-edge** → more extension.

Blocked result for MOR back-arcs:

| Predictor | ρ vs bab_zone | Sign vs Schellart |
|---|---:|---|
| Width | **+0.22** | **opposite** |
| Edge distance | **+0.14** | **opposite** |

Wide zones host more reconstructed back-arcs (P = 0.17) than narrow ones (P = 0.045). Dilatation layer does not survive time dummies. Against independent SubMap C vs E at 0 Ma, width is weak (ρ = −0.19, p = 0.10).

---

## Results — independent geological UPS

**22** cited intervals · **295** segment×time blocks · outcome is **not** Müller `bab_on`.

| Test | Result |
|---|---|
| Median age E vs C | **61 vs 80 Ma** (E *younger*; p = 0.04) |
| 55 Ma gate | P(E\|≥55) = **0.54** vs P(E\|<55) = **0.69**; **OR = 0.5** (opposite Sdrolias) |
| Width vs is_E | ρ = **−0.04**, p = 0.50 (null; Müller zone widths glue IBM–Japan–Kurils) |
| 0 Ma vs SubMap | relaxed agreement **47/83 = 0.57** |

West Philippine / early IBM extension (55–30 Ma) sits on ~17–40 Ma Pacific in these age grids; Andes shortening sits on older Nazca. Dropping Tethyan Hellenic/Calabria makes E even younger (52 vs 80 Ma).

**No hypothesis winner.** Sparse hand catalog; C is almost entirely Central Andes.

---

## Interpretation

1. **0 Ma snapshots overstate order.** Age, H1, and H2 all look promising until time, blocking, or frame choice is applied.  
2. **Müller2019 cannot separate trench rollback from upper-plate motion.** That test must stay on SubMap (or another absolute catalog).  
3. **Reconstruction-internal back-arc flags are not geology.** They are useful for a Sdrolias-style association and for filter QC (Mariana/Tonga on; Andes off), not for crowning H2.  
4. **Independent geology agrees that the 0 Ma age story is fragile.** Among published Cenozoic intervals, geological extension is *younger* than Andean shortening in this trench-age catalog — opposite the Sdrolias gate.  
5. **The honest paper** is that present-day correlations are fragile once time, frame, geometry, and independent UPS are treated properly. A major contribution is showing *where* they break.  

---

## What would make a stronger claim

- Expand geological UPS beyond Andes-vs-Pacific-BABs (more C intervals; better width metric than connected trench polylines)  
- A multivariate C vs E model (age + \(V_{\mathrm{UP}}\) + \(V_{\mathrm{T}}\) + width) blocked by trench  
- Named SubMap frames 2–4 / 5–9 (Cerpa lists four M56 frames but does not map column numbers)  
- A second plate model so results are not Müller-only  

---

## Key data products

| Product | Path |
|---|---|
| Proposal | `docs/subduction-dynamics-proposal.md` |
| 0 Ma trench kinematics | `results/trench_kinematics_0Ma_Muller2019.csv` |
| 0–60 Ma trench catalog | `results/trench_kinematics_0-60Ma_Muller2019.csv` |
| SubMap H1/H2 | `results/submap_m56_h1h2.md` |
| Frames sidecar | `results/submap_h1h2_frames.md` |
| Age gate | `results/cenozoic_age_gate.md` |
| Schellart test | `results/schellart_width_edge.md` |
| Geological UPS | `data/geological_ups/` · `results/geological_ups_test.md` |
| Slide deck | `docs/presentation.md` |

---

## Selected references

- Uyeda and Kanamori (1979), *JGR*  
- Heuret and Lallemand (2005), *PEPI*; Lallemand, Heuret and Boutelier (2005), *G3*  
- Schellart et al. (2007), *Nature*; Schellart (2008), *PEPI*; Schellart (2024), *ESR*  
- Sdrolias and Müller (2006), *G3*  
- Müller et al. (2019), *Tectonics*  
- Cerpa, Lallemand and Heuret (2025), EarthArXiv / SubMap v7  
- Mather et al. (2023), GPlately  

---

## Thanks / discussion

Questions welcome — especially on:

- how independent the geological UPS intervals should be from EarthByte reconstructions  
- whether to add a second plate model before writing  
- slide vs paper next  

