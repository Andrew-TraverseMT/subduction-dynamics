# Geological upper-plate strain (UPS) intervals — Cenozoic (0–60 Ma)

Hand-compiled from **published geological and marine-geophysical opening/cessation ages** of back-arc basins and from standard reviews of Central Andean shortening. Ages are **not** taken from Müller et al. (2019) topologies.

File: `cenozoic_ups_intervals.csv` (22 intervals, 14 named trench segments).

This table is sparse and uncertain by construction. Use midpoints stored in `t_start_Ma` / `t_end_Ma`; ranges and conflicts live in `age_basis` and `notes`. An interval is active at reconstruction time \(t\) when `t_end_Ma ≤ t ≤ t_start_Ma`. If two intervals hit the same trench sample, the join script prefers **E then C then N**.

## What this is / is not

- **Is:** an independent geological UPS time series for a handful of well-dated upper plates, joined to Müller2019 *trench kinematics* (age, width, edge) as predictors.
- **Is not:** a global UPS catalog. Most trench sample-times stay `Unknown`.
- **Is not:** Müller back-arc MORs (`bab_on`) or deforming-network dilatation. Those are reconstruction-internal and are not used as the outcome.

`ups_class`: **E** = documented back-arc spreading or rifting; **C** = documented cordilleran / arc shortening; **N** = documented absence of a back-arc, or a published neutral-stress interval. Cascadia is **N**, not Andean **C** (see below).

Bboxes are in degrees, east-longitude (−180 to 180). If `lon_min > lon_max` the box wraps across 180° (Tonga, Kermadec). Bboxes are meant to follow the *trench* in Müller2019 `results/trench_kinematics_0-60Ma_Muller2019.csv`, not the basin floor. False-positive risk is noted per interval.

## Named segments

Andes, Andaman, Calabria, Cascadia, Hellenic, Izu-Bonin, Japan, Kermadec, Manus, Mariana, NewHebrides, Ryukyu, Scotia, Tonga.

These are the independent units in `scripts/geological_ups_test.py`. They are **not** Müller `zone_id`.

---

## Extension (E)

### West Philippine Basin — segments Mariana + Izu-Bonin — 55–30 Ma

Spreading behind proto-IBM (Palau–Kyushu Ridge).

| Choice | 55–30 Ma |
|---|---|
| Primary | Deschamps & Lallemand (2002): rifting ~55 Ma, spreading 54 to 33/30 Ma |
| Compilation | Sdrolias & Müller (2006): 55 to 30–33 Ma |
| Conflict | Hilde & Lee (1984): 58–33 Ma (two phases). Sasaki et al. (2014): cessation ~36 Ma |

- Deschamps, A., & Lallemand, S. (2002). *J. Geophys. Res.*, 107(B12). https://doi.org/10.1029/2001JB001706
- Hilde, T. W. C., & Lee, C.-S. (1984). *Tectonophysics*, 102, 85–104. https://doi.org/10.1016/0040-1951(84)90009-X
- Sdrolias, M., & Müller, R. D. (2006). *Geochem. Geophys. Geosyst.*, 7, Q04016. https://doi.org/10.1029/2005GC001090

Bboxes are west-shifted relative to the modern Mariana Trench because proto-IBM in this catalog sits near 135–143°E at 50–55 Ma. Not matched to Ryukyu.

### Parece Vela Basin — Mariana — 30–15 Ma

- Sdrolias, M., Roest, W. R., & Müller, R. D. (2004). *Tectonophysics*, 394, 69–86. https://doi.org/10.1016/j.tecto.2004.07.003 — rifting ~30 Ma, spreading by chron 9o (28 Ma) / possibly 10o (29 Ma), cessation ~15 Ma.
- Okino, K., Kasuga, S., & Ohara, Y. (1998). *Mar. Geophys. Res.*, 20, 21–40. https://doi.org/10.1023/A:1004377422118
- Mrozowski, C. L., & Hayes, D. E. (1979). Magnetic chrons 10–5D (~30–17 Ma).

Conflict: Maffione et al. (2022, IODP U1438) delay *true* spreading to ~24–22 Ma after 30–24 Ma rifting. Rifting + spreading are both E here.

### Shikoku Basin — Izu-Bonin — 29–15 Ma

Twin of Parece Vela. Chamot-Rooke et al. (1987) / Okino et al. (1994) often quote 27–15 or 26–15 Ma spreading; Sdrolias et al. (2004) start rifting at ~30 Ma. Used 29–15. Nankai is excluded (`lon_min = 138`).

- Chamot-Rooke, N., Renard, V., & Huchon, P. (1987). *Earth Planet. Sci. Lett.*, 83, 214–228. https://doi.org/10.1016/0012-821X(87)90071-9

### IBM gap after Parece Vela — Mariana — N, 15–7 Ma

Not a dated “neutral orogeny”. Bracketed by PV cessation (~15 Ma) and Mariana Trough rifting (~7 Ma). Labelled N (no active back-arc), not C.

### Mariana Trough — Mariana — 7–0 Ma

| Choice | 7–0 Ma (rifting + spreading) |
|---|---|
| DSDP 60 | Hussong & Uyeda (1981): rifting latest Miocene ~6 Ma |
| Marine geophysics | Martinez et al. (1995): ~7 Ma; Yamazaki et al. (2003): spreading at 18°N from ~6 Ma |
| Conflict | Some reviews put spreading at 3–4 Ma; Stern et al. (2003) rifting after ~10 Ma, spreading ~5 Ma |

- Martinez, F., Fryer, P., Baker, N. A., & Yamazaki, T. (1995). *J. Geophys. Res.*, 100, 3807–3827. https://doi.org/10.1029/94JB02466
- Yamazaki, T., Seama, N., Okino, K., et al. (2003). *Geochem. Geophys. Geosyst.*, 4, 1073. https://doi.org/10.1029/2002GC000492

Modern box `142–152°E, 10–24°N` excludes Yap (~138°E) and Izu-Bonin (lat ≥ 24°).

### Izu-Bonin after Shikoku — Izu-Bonin — N, 15–0 Ma

No organized back-arc after Shikoku cessation. Sumisu/Torishima rifts are nascent intra-arc rifts, not treated as a BAB.

### Lau Basin — Tonga — 6–0 Ma

- Parson, L. M., & Hawkins, J. W. (1994). In: *Proc. ODP, Sci. Results*, 135.
- Ruellan, E., Delteil, J., Wright, I., & Matsumoto, T. (2003). *Geochem. Geophys. Geosyst.*, 4, 9105. https://doi.org/10.1029/2001GC000261 — spreading 6–5.5 Ma in the north, southward propagation.
- Taylor, B., Zellmer, K., Martinez, F., & Goodliffe, A. (1996). Three-plate kinematics (Brunhes).
- Hawkins, J. W. (1995). ODP Leg 135.

Wrap box across 180°. Split from Kermadec at 26°S.

### Havre Trough — Kermadec — 5.5–0 Ma (rifting)

- Caratori Tontini, F., Bassett, D., de Ronde, C. E. J., Timm, C., & Wysoczanski, R. (2019). *Nat. Geosci.* https://doi.org/10.1038/s41561-019-0473-9 — split ~5.5 Ma; early spreading burst then disorganised rifting.
- Wysoczanski, R. J., et al. (2019). *N. Z. J. Geol. Geophys.* https://doi.org/10.1080/00288306.2019.1602059 — southern Havre Ar–Ar <2 Ma, coeval with the Taupo Volcanic Zone.

Conflict: the south may be much younger than 5.5 Ma. Used 5.5–0 for the Kermadec upper plate as a whole. This is rifting, not a mature MOR.

### East Scotia Ridge — Scotia — 16–0 Ma

- Larter, R. D., Vanneste, L. E., Morris, P., & Smythe, D. K. (2003). *Geochem. Geophys. Geosyst.*, 4, 1070. https://doi.org/10.1029/2001GC000238 — spreading underway by 15–17 Ma; possibly C6 (~20 Ma).
- Eagles, G., Livermore, R. A., Fairhead, J. D., & Morris, P. (2005). *J. Geophys. Res.*, 110, B02401. https://doi.org/10.1029/2004JB003154 — West Scotia Ridge ~26.5–6 Ma (not assigned here).
- Eagles, G., & Jokat, W. (2014). East Scotia Ridge onset ~17 Ma.

Used 16–0 (midpoint of 15–17 vs ~20). West Scotia opening is *not* mapped onto the Sandwich trench box. In Müller2019 the Sandwich trench itself only appears after ~15 Ma.

### Central Andaman Basin — Andaman — 4–0 Ma

| Choice | 4–0 Ma (range 4–6) |
|---|---|
| Spreading | Kamesh Raju et al. (2004); Curray (2005): ~4 Ma |
| Revision | Jha et al. (2021): chron C3An.1ny = 5.894 Ma |
| Conflict | Morley & Alvey (2015): sedimentation incompatible with continuous 4 Ma spreading |

- Kamesh Raju, K. A., et al. (2004). *Earth Planet. Sci. Lett.*, 221, 145–162. https://doi.org/10.1016/j.epsl.2004.02.022
- Curray, J. R. (2005). *J. Asian Earth Sci.*, 25, 187–232. https://doi.org/10.1016/j.jseaes.2004.09.001
- Jha, P., et al. (2021). *Tectonophysics*. https://doi.org/10.1016/j.tecto.2021.229085

Spreading is a leaky, highly oblique transform. **SubMap geodetic UPS is SS, not E** — an expected 0 Ma mismatch, not a join bug.

### Okinawa Trough — Ryukyu — E 10–6 Ma and E 2–0 Ma

Continental rifting, not MOR (Arai et al., 2017: still rifting).

- Sibuet, J.-C., et al. (1987). *J. Geophys. Res.*, 92, 14041–14063. https://doi.org/10.1029/JB092iB13p14041 — first phase middle–late Miocene, ceased in Pliocene; second phase from ~2 Ma.
- Letouzey, J., & Kimura, M. (1986). *Mar. Geol.*, 71, 57–91. https://doi.org/10.1016/0025-3227(86)90080-5
- Miki, M. (1995). Two-phase model: 10–6 Ma then ~1 Ma.

Gap 6–2 Ma is left **unclassified** (not forced to N). Timing is debated (Kimura 1996: 6–4 Ma).

### Japan Sea — Japan — E 28–14 Ma

- Tamaki, K., Suyehiro, K., Allan, J., Ingle, J. C., Jr., & Pisciotto, K. A. (1992). *Proc. ODP, Sci. Results*, 127/128 — extending 32–10 Ma; major opening 28–18 Ma; Ar–Ar basement 24–17 Ma.
- Jolivet, L., Tamaki, K., & Fournier, M. (1994). *J. Geophys. Res.*, 99, 22237–22259. https://doi.org/10.1029/93JB03463 — stopped ~15–12 Ma.
- Sato, H. (1994). *J. Geophys. Res.*, 99, 22261–22274. https://doi.org/10.1029/94JB00854 — extension from 32 Ma; termination ~14 Ma.

Used 28–14. **Sparsity:** Müller2019 Japan Trench samples are almost absent before ~15 Ma, so this E interval is poorly sampled in the kinematics catalog.

### Aegean rifting — Hellenic — 23–0 Ma (not a MOR)

Onset is genuinely contested:

| School | Onset |
|---|---|
| Le Pichon & Angelier (1979) | ~13–15 Ma (basins) |
| Ring et al. (2010) | ~23–19 Ma Aegean-wide lithospheric extension |
| Jolivet & Brun (2010), Brun & Sokoutis (2010) | ~45–35 Ma rollback |
| Some 2022 papers | regional extension only after ~15 Ma |

Used **23–0** (Ring et al. 2010) as a cited midpoint of the 15 vs 45 Ma debate. Most of the *magnitude* of extension is post-15 Ma.

- Ring, U., Glodny, J., Will, T., & Thomson, S. (2010). *Annu. Rev. Earth Planet. Sci.*, 38, 45–76. https://doi.org/10.1146/annurev.earth.012809.104804
- Jolivet, L., & Brun, J.-P. (2010). *Int. J. Earth Sci.*, 99, 109–124. https://doi.org/10.1007/s00531-008-0356-6
- Le Pichon, X., & Angelier, J. (1979). *Tectonophysics*, 60, 1–42.

Calabria is a separate segment (`lon < 18.5`).

### Tyrrhenian — Calabria — 10–0 Ma

- Rosenbaum, G., & Lister, G. S. (2004). *Tectonics*, 23, TC1018. https://doi.org/10.1029/2003TC001518 — major rifting from ~10–9 Ma (Tortonian).
- Kastens, K., et al. (1988). *GSA Bull.*, 100, 1140–1156. https://doi.org/10.1130/0016-7606(1988)100<1140:OLITTS>2.3.CO;2 — ODP Leg 107; Vavilov before Marsili.
- Malinverno, A., & Ryan, W. B. F. (1986). Rollback model.

Vavilov spreading ~4.3–2.6 Ma, Marsili ~2–1 Ma (still rifting+spreading = E from 10 Ma).

### Manus Basin — Manus (New Britain Trench) — 3.5–0 Ma

- Taylor, B. (1979). *Geology* — MSC/ETZ spreading <3.5 Ma.
- Martinez, F., & Taylor, B. (1996). *Mar. Geophys. Res.*, 18, 203–224. https://doi.org/10.1007/BF00286078 — Brunhes reorganisation; East Manus rifting <0.78 Ma.

### North Fiji Basin — NewHebrides — 12–0 Ma

- Auzende, J.-M., Lafoy, Y., & Marsset, B. (1988); Auzende et al. (1995) in *The North Fiji Basin*. Opening from ~12 Ma in three stages after Vitiaz lockup / polarity reversal.
- Martin, A. K. (2013). *Tectonophysics*. https://doi.org/10.1016/j.tecto.2013.06.019

The d’Entrecasteaux collision locally compresses central Vanuatu; that along-strike exception is not split out.

---

## Compression (C)

### Central Andes shortening — Andes — 50–0 Ma

Altiplano–Puna latitudes only (`−32 to −12°`, west coast of South America).

- Oncken, O., et al. (2006). In: *The Andes — Active Subduction Orogeny*. Springer. https://doi.org/10.1007/978-3-540-48684-8_1 — shortening from 45–50 Ma, rate increase ~32–27 Ma, still active.
- Oncken, O., Boutelier, D., Dresen, G., & Schemmann, K. (2012). *Geochem. Geophys. Geosyst.*, 13, Q12007. https://doi.org/10.1029/2012GC004280
- Horton, B. K. (2018). *Annu. Rev. Earth Planet. Sci.*, 46, 213–244. https://doi.org/10.1146/annurev-earth-063016-020612 — Paleogene onset, eastward migration.
- Garzione, C. N., et al. (2008). *Science*, 320, 1304–1307. https://doi.org/10.1126/science.1148615 — rapid paleoelevation rise is *not* the UPS classifier; shortening is.
- Schellart, W. P. (2024). *Earth-Sci. Rev.*, 253, 104755. https://doi.org/10.1016/j.earscirev.2024.104755 — Andes as the wide-slab shortening end-member vs Scotia extension.

Incaic (late Eocene) and Quechua (Miocene) pulses are phases of one Cenozoic shortening, not a switch to E. Conflict: some authors put most *uplift* after ~10–6 Ma (Garzione) or most *growth* after ~25/~15 Ma; the onset of shortening is still Paleogene (Oncken, Horton, Barnes & Ehlers 2009). Used 50–0 (onset range 60–40 Ma).

Northern Andes (Colombia/Ecuador) and Patagonia/south-central Chile are outside the box.

### NE Japan inversion — Japan — C 4–0 Ma

- Sato, H. (1994). https://doi.org/10.1029/94JB00854 — strong E–W compression from ~4 Ma; inversion of Miocene Japan Sea normal faults.

This is arc/back-arc inversion, not Andean plateau building. Some papers place localized inversion as young as ~1.8 Ma.

Between Japan Sea opening and inversion, Sato (1994) describes a **neutral** stress regime: `Japan_post_opening_neutral` N, 14–4 Ma.

---

## Neutral (N), labelled carefully

### Cascadia — N 30–0 Ma — no back-arc

**Not C.** Cascadia is not an Andean retroarc.

- Wells, R. E., Weaver, C. S., & Blakely, R. J. (1998). *Geology*, 26, 759–762. https://doi.org/10.1130/0091-7613(1998)026<0759:FAMICA>2.3.CO;2 — Oregon rotates clockwise; western Washington shortens N–S against the Canadian Coast Mountains (margin-parallel, not trench-normal).
- Johnson, S. Y., et al. (2004). *Tectonics*, 23, TC1011. https://doi.org/10.1029/2003TC001507

Southern Cascadia’s volcanic arc sits on the trailing edge of the rotating Oregon block (Basin and Range extension). That is **not** a Cascadia back-arc basin. SubMap `Geodetic_UPS` for Cascadia is **SS**, which we treat as agreeing with N, not with C.

Müller2019 only has a Cascadia trench from ~30 Ma in this catalog.

---

## How to match trench samples

`scripts/geological_ups_test.py` assigns a sample if its (`lon`, `lat`, `time_Ma`) falls in an active interval bbox. Otherwise `Unknown`. Named `segment` comes from the matched interval.

False-positive watch-list:

| Risk | Mitigation |
|---|---|
| Yap labelled as Mariana Trough | Mariana 7–0 Ma box starts at 142°E |
| Nankai labelled as Izu-Bonin / Shikoku | Izu-Bonin `lon_min` 138–139 |
| Ryukyu labelled as WPB | WPB is Pacific-facing proto-IBM; Ryukyu is a separate segment |
| Northern Sumatra as Andaman | `lat_min = 5` |
| Colombia or Patagonia as Central Andes | Andes box `−32 to −12°` |
| West Scotia as Sandwich back-arc | only East Scotia is assigned |
| Gulf of California as MAT back-arc | MAT is not in this table |

## Cite (minimum)

Geological ages: papers listed above (prefer the primary marine-geophysical / structural paper; Sdrolias & Müller 2006 is a secondary compilation).

Trench kinematics (predictors only): Müller, R. D., et al. (2019). *Tectonics*, 38, 1884–1907. https://doi.org/10.1029/2018TC005462

Present-day geodetic UPS (0 Ma sanity only): Cerpa, N., Lallemand, S., & Heuret, A. (2025). EarthArXiv. https://doi.org/10.31223/x56m8h — SubMap v7.
