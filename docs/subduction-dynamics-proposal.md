# Competing kinematic controls on upper-plate stress at subduction zones: a time-dependent PyGPlates test

**Working title.** A reconstruction-based test of Heuret–Lallemand vs Schellart hypotheses for compressive, neutral, and extensional upper-plate regimes.

**Investigator.** Grok, Tectonics Researcher (virtual graduate project).  
**Advisor.** Andrew Laskowski.  
**Tools.** pyGPlates / GPlately, SubMap (present-day ground truth), Cenozoic plate reconstructions.  
**Status.** Proposal for approval before coding.

---

## 1. Problem

Since Uyeda and Kanamori (1979), the first-order puzzle has been why some overriding plates are Chilean (compressive) and others Mariana (extensional). The last twenty years of global compilations did not settle it. They produced competing, observationally supported stories that still disagree on which driver actually sets upper-plate state of stress.

Those stories were built almost entirely from **present-day snapshots**: sparse trench-normal transects, one absolute reference frame at a time, and kinematic predictors that are strongly collinear (age, convergence, trench motion, upper-plate motion, slab width, distance to a slab edge). With current plate models, seafloor-age grids, and pyGPlates, we can rebuild the same experiment **along the trench and through time**, and ask which hypothesis still wins when the kinematics change.

That is the contribution. Not a new cartoon of subduction. A state-of-the-art, time-dependent, multi-model statistical test of the hypotheses the field already argues over.

---

## 2. Competing hypotheses (what we will actually test)

We treat upper-plate state of stress (UPS) as compressive / neutral / extensional (and, where rates exist, as a signed trench-normal deformation velocity \(v_{\mathrm{OPD}\perp}\)). Predictors are computed at densely sampled trench points.

**H1. Upper-plate absolute motion (Heuret and Lallemand, 2005; Lallemand, Heuret, and Boutelier, 2005).**  
Trench-normal absolute velocity of the overriding plate (\(v_{\mathrm{OP}\perp}\)) is the leading control. The slab is partly anchored (trench motion globally limited, \(\lvert v_T\rvert \lesssim 50\) mm yr\(^{-1}\) in their compilation). Upper-plate retreat \(\rightarrow\) extension; advance \(\rightarrow\) compression. Slab dip is an intermediary, not an independent driver: deep dips \(>50^\circ\) go with spreading, dips \(<30^\circ\) with shortening. Slab age, slab-pull magnitude, and convergence rate do **not** correlate with dip.

**H2. Slab-driven trench migration and slab width (Schellart et al., 2007; Schellart, 2008; Schellart, 2024).**  
Trench-normal trench velocity (\(v_{T\perp}\)) is the leading kinematic control (retreat \(\rightarrow\) extension). Rollback is a slab-buoyancy / return-flow problem: narrow slabs and segments near lateral edges retreat; the centers of wide slabs (\(\gtrsim 4000\) km) are slow or advancing. Distance to the nearest slab edge, and slab width, should outperform \(v_{\mathrm{OP}\perp}\). This prediction is **reference-frame dependent**, which is a feature of the test, not a nuisance.

**H3. Subducting-plate age / slab pull (Molnar and Atwater, 1978; Sdrolias and Müller, 2006).**  
Older lithosphere is more negatively buoyant. Sdrolias and Müller found back-arc basins only where subducting normal oceanic lithosphere is \(\gtrsim 55\) Ma, with intermediate slab dip \(>30^\circ\). They also found that upper-plate motion *away from the hinge* precedes basin opening, after which rollback can maintain extension even if the upper plate stops retreating. We test age as (a) a continuous predictor and (b) a threshold that gates the other hypotheses.

**H4. Convergence rate and obliquity.**  
Trench-normal convergence is often invoked as a coupling / compression control. Lallemand et al. (2005) found no correlation with dip; we still include trench-normal and trench-parallel components so H4 can lose cleanly. Cerpa, Lallemand, and Heuret (2025, EarthArXiv / SubMap update) argue that convergence obliquity has a limited effect on strain partitioning, while deformation obliquity still tracks slab dip. We keep obliquity in the model so 1-D trench-normal UPS is not silently contaminated by strike-slip.

**H5. Interaction, not a single winner.**  
The productive reading of this literature is that no univariate story is sufficient. Sdrolias and Müller already combined an age threshold with a kinematic precursor. Schellart combined \(v_{\mathrm{OP}}\) with distance to a slab edge. Lallemand et al. combined \(v_{\mathrm{OP}}\) with slab geometry and upper-plate crustal nature (oceanic vs continental). The statistical design below is built to let interactions win if they should.

Null against which all of this is judged: UPS is indistinguishable from a spatially correlated random field once trench segment identity is accounted for.

---

## 3. Why this is still open (and why PyGPlates changes the experiment)

Present-day compilations (Heuret and Lallemand, 2005; Lallemand et al., 2005; Schellart, 2008, 2024; SubMap, including v7 in 2025) sample on the order of 150–260 transects. Adjacent samples are not independent. Absolute velocities flip meaning between hotspot, no-net-rotation, and mantle frames (Schellart, 2008). Age, width, and \(v_{\mathrm{OP}}\) co-vary geographically (old, steep, retreating western Pacific vs young, flat, compressive eastern Pacific), so pairwise correlations cannot isolate a driver.

Sdrolias and Müller (2006) already took the Cenozoic step, with paleo-age grids and a moving Indo-Atlantic hotspot frame. That paper asked a coarser question (back-arc basin present or absent) with an older plate model. It did not run a formal horse-race among H1–H4, did not treat compressive vs neutral vs extensional as a three-way regime, and could not use deforming-plate topologies.

What is new here:

1. **Dense, repeatable kinematics** from pyGPlates / GPlately (`tessellate_subduction_zones`, subduction-convergence, age-grid sampling, distance to trench edge) instead of hand transects.
2. **Time.** Rebuild the Heuret-style experiment every 1 Myr through the Cenozoic, not just at 0 Ma.
3. **Multiple reconstruction models and absolute frames**, so a result that exists in only one model or one frame is reported as fragile.
4. **Independent present-day UPS** from SubMap (geodetic UPS classes and rates; seismic UPS as it becomes available), plus a time-dependent UPS proxy from reconstructed back-arc opening/closing.
5. **A nested statistical test** designed for collinear, spatially autocorrelated predictors, not a gallery of scatterplots.

Slab dip will not be treated as a primary time-dependent predictor. It is not a pyGPlates output. At 0 Ma we can join Slab2.0. In the past it is either omitted, held as a present-day attribute, or used only as an intermediary in the 0 Ma validation against Lallemand et al. (2005). That limitation stays explicit.

---

## 4. Methods

### 4.1 Kinematic engine

Python stack: **pyGPlates** for topology-resolved rotations and (where models include them) deforming networks; **GPlately** as the high-level interface (Mather et al., 2023). Default reconstruction: Müller et al. (2019) deforming model. Robustness set: Müller et al. (2016) and a more recent Clennett / Müller-family model available through GPlately’s Plate Model Manager (e.g. Müller2022 if the topologies are fit for trench sampling). Absolute frames: no-net-rotation, an Indo-Atlantic moving-hotspot frame, and the model’s native mantle frame where one exists.

At each time \(t = 0, 1, \ldots, 80\) Ma, tessellate every subduction zone at \(\sim 50\)–100 km. At each sample compute:

| Predictor | Hypothesis |
|---|---|
| Seafloor age at trench (age grid) | H3 |
| Trench-normal and trench-parallel convergence | H4 |
| Convergence obliquity | H4 |
| Trench-normal trench velocity \(v_{T\perp}\) | H2 |
| Trench-normal upper-plate absolute velocity \(v_{\mathrm{OP}\perp}\) | H1 |
| Subducting-plate trench-normal absolute velocity | H2/H3 |
| Distance to nearest lateral slab edge; trench-parallel slab width | H2 |
| Upper-plate nature (oceanic vs continental) | Lallemand et al. (2005) covariate |
| Orthogonal subduction velocity and slab flux (age \(\rightarrow\) thickness) | H3/H4 |

Collision, ridge/plateau/seamount subduction, and incipient zones will be flagged and run both in and out of the sample, matching the exclusion rules in Heuret and Lallemand (2005) and Lallemand et al. (2005).

### 4.2 Upper-plate state of stress

**Present-day (calibration).** SubMap upper-plate strain: compressive / extensional / neutral, plus geodetic deformation rates and deformation obliquity from the 2025 kinematics update (Cerpa et al.). Join by nearest trench sample, not by recycling their 200 km transect geometry, so the kinematic side is independent of their sampling.

**Through time (the actual test).** Two UPS proxies, kept separate:

1. *Geologic / kinematic:* opening and death of back-arc basins already encoded in the reconstruction (the Sdrolias–Müller observable, but as a time series). Extension where a back-arc spreading center is active; compression where a known fold-and-thrust / cordilleran regime is reconstructed; else neutral/unknown.
2. *Model-internal (secondary, not independent):* dilatation rate sampled \(\sim 50\) km into the overriding plate from GPlates deforming networks (pyGPlates sample workflow). This is a consistency check against the plate model’s own deformation, not an observation.

### 4.3 Statistical design

Present-day first, as a replication of Heuret / Lallemand / Schellart, then the time-dependent experiment.

- Univariate and partial correlations, with a spatial block bootstrap along trench so \(N\) is not the raw sample count.
- Nested classification / regression for C / N / E and for \(v_{\mathrm{OPD}\perp}\): H1-only, H2-only, H3-only, H4-only, then additive combinations, then the pre-registered interactions (age \(\times\) \(v_{\mathrm{OP}\perp}\); slab width \(\times\) \(v_{T\perp}\); \(v_{\mathrm{OP}\perp}\) \(\times\) distance to edge).
- Model comparison by AIC/BIC and by out-of-sample log loss, blocked by subduction zone so we are not just interpolating along a trench.
- Repeat the entire ranking in each reconstruction and each absolute frame. A hypothesis “wins” only if its rank is stable.

**Pre-registered success.** H1 wins if \(v_{\mathrm{OP}\perp}\) remains the best single predictor of UPS after controlling for age, width, and frame. H2 wins if \(v_{T\perp}\) and distance-to-edge beat \(v_{\mathrm{OP}\perp}\) in frames that minimize dissipation / wide-slab trench motion (Schellart’s preferred class of frames) **and** still beat it in at least one other frame. H3 wins if the \(\sim 55\) Ma age gate is necessary for extension through the Cenozoic, even when H1/H2 kinematics look favorable. H5 wins if an interaction model is clearly preferred and the univariate ranking is unstable in time.

### 4.4 Implementation

All coding on this machine, in Python, against pyGPlates/GPlately. Reproducible repo: notebooks for figures, scripts for the kinematic extraction and stats, and a frozen present-day join to SubMap. No new plate model will be built in this phase.

---

## 5. Expected contribution

A paper that can say, with numbers and with time, whether the Heuret–Lallemand upper-plate story, the Schellart slab-width/rollback story, the Sdrolias–Müller age gate, or a specific combination is the one that survives once the present-day co-location of “old, steep, western, retreating” is broken by Cenozoic kinematic change. That is the major contribution the new tools actually buy. Secondary products: a public, time-stamped trench-kinematic catalog (age, convergence, obliquity, \(v_T\), \(v_{\mathrm{OP}}\), width, edge distance) that later analog/numerical papers can use instead of another 0 Ma transect file.

---

## 6. Scope, risks, and what this is not

- **Not** a dynamic subduction model. We are testing kinematic hypotheses against UPS, the same empirical game Heuret, Lallemand, and Schellart played, with better sampling and time.
- **Not** a slab-dip evolution paper. Dip enters only at 0 Ma.
- Main risks: (i) reconstructions bake in the back-arc histories we would like to predict, so the time-dependent UPS proxy is partly circular — the present-day SubMap test is the clean one, the Cenozoic test is a consistency test with that caveat stated; (ii) absolute frame choice can manufacture a winner, which is why the ranking must be shown in several frames; (iii) SubMap licensing / download of the full parameter table may need a manual step.
- First milestone after approval: 0 Ma kinematic extraction from Müller2019, joined to SubMap UPS, with a replication figure of the classic \(v_{\mathrm{OP}\perp}\) and \(v_{T\perp}\) vs strain-class plots.

---

## Key references

Cerpa, N. G., Lallemand, S., & Heuret, A. (2025). Unraveling the effect of convergence obliquity on overriding-plate deformation and strain partitioning in subduction zones. EarthArXiv. https://doi.org/10.31223/x56m8h  

Heuret, A., & Lallemand, S. (2005). Plate motions, slab dynamics and back-arc deformation. *PEPI*, 149, 31–51.  

Lallemand, S., Heuret, A., & Boutelier, D. (2005). On the relationships between slab dip, back-arc stress, upper plate absolute motion, and crustal nature in subduction zones. *G3*, 6, Q09006.  

Mather, B. R., et al. (2023). Deep time spatio-temporal data analysis using pyGPlates with PlateTectonicTools and GPlately. *Geoscience Data Journal*.  

Molnar, P., & Atwater, T. (1978). Interarc spreading and Cordilleran tectonics as alternates related to the age of subducted oceanic lithosphere. *EPSL*, 41, 330–340.  

Müller, R. D., et al. (2019). A global plate model including lithospheric deformation along major rifts and orogens since the Triassic. *Tectonics*, 38, 1884–1907.  

Schellart, W. P., et al. (2007). Evolution and diversity of subduction zones controlled by slab width. *Nature*, 446, 308–311.  

Schellart, W. P. (2008). Subduction zone trench migration: Slab driven or overriding-plate-driven? *PEPI*, 170, 170–190.  

Schellart, W. P. (2024). Subduction dynamics and overriding plate deformation. *Earth-Science Reviews*, 104755.  

Sdrolias, M., & Müller, R. D. (2006). Controls on back-arc basin formation. *G3*, 7, Q04016.  

Uyeda, S., & Kanamori, H. (1979). Back-arc opening and the mode of subduction. *JGR*, 84, 1049–1061.  
