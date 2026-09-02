# Subduction-zone upper-plate strain (UPS) data bundle

Downloaded 2026-09-02 (MDT). Nothing here is fabricated: values come from the sources named below. Missing numeric fields use SubMap’s sentinel `-999`.

## Which source succeeded (present-day global UPS table)

**SubMap Sub-DATA v7 REST API succeeded** (source 2). It is the present-day Heuret & Lallemand / SubMap transect database, including **geodetic UPS**, trench-point coordinates, `vd`, and deformation obliquity `Ovd`. Kinematics in this dump are the **M56+** (MORVEL56 + arc-block) set described by Cerpa, Lallemand & Heuret (EarthArXiv 2025 / Tektonika revision).

The SubMap **website UI is a clicky Vue GIS**. There is **no single static global CSV/SHP file** advertised for anonymous download. Map/section PNG/EPS/PDF and rupture SHP exports are generated after interactive map jobs. Sub-DATA’s “Download .CSV FILE” button is **per transect** (client-side JSON→CSV of one sheet).

A browser is **not** required for the global table: the same JSON the UI loads is public, no login:

- `GET https://submap.gm.umontpellier.fr/api/transects/` — 260 short names + transect endpoints (Lon1/Lat1–Lon2/Lat2)
- `GET https://submap.gm.umontpellier.fr/api/transects/{Short_Name}/` — full sheet (~221 fields), including `Geodetic_UPS`, `M56_vd`, `M56_Ovd`

Those 260 sheets were pulled and compiled into the CSVs below.

## Primary files (use these)

| Path | What | n |
|---|---|---|
| `submap_geodetic_UPS.csv` | Compact UPS / kinematics table | **260 rows** + header |
| `submap_transects_full.csv` | All 221 Sub-DATA fields | **260 rows** |
| `submap_transects_full.json` | Same records as JSON array | 260 |
| `submap_transect_index.json` | List endpoint (name + endpoints only) | 260 |
| `submap_api_raw/{Short_Name}.json` | One raw API sheet per transect | 260 files |

### `submap_geodetic_UPS.csv` columns (34)

`Short_Name, Zone, id, Trench_Name, Lon, Lat, Lon1, Lat1, Lon2, Lat2, Az, UPN, Geodetic_UPS, Seismic_UPS, M56_SP_Name, M56_UP_Name, M56_ARC_Name, M56_ARC_ref, M56_vc, M56_Azvc, M56_Ovc, M56_vcn, M56_vs, M56_Azvs, M56_Ovs, M56_vsn, M56_vd, M56_Azvd, M56_Ovd, M56_vdn, alphaS, alphaD, A, Phi`

Meanings (from SubMap / Cerpa et al. 2025):

- `Lon, Lat` — trench point at the transect (degrees; lon 0–360 in some zones)
- `Lon1,Lat1,Lon2,Lat2,Az` — transect endpoints and azimuth
- `Geodetic_UPS` — present-day upper-plate strain class from `Ovd` / `vd`
- `M56_vd, M56_Azvd, M56_Ovd, M56_vdn` — deformation-velocity magnitude, azimuth, obliquity, trench-normal component (mm/yr, degrees)
- `M56_vc*` / `M56_vs*` — convergence / subduction velocities (same convention)
- `Seismic_UPS` — **all `-999`** in this dump (seismic UPS “available soon”, Bonnamy et al. in review)

### Geodetic_UPS counts (n=260)

| Class | n | notes |
|---|---:|---|
| C | 7 | pure compression (`Ovd` small) |
| C-SS | 33 | dominant compression + strike-slip |
| SS-C | 59 | dominant strike-slip + compression |
| SS | 37 | pure strike-slip |
| SS-E | 38 | dominant strike-slip + extension |
| E-SS | 32 | dominant extension + strike-slip |
| E | 7 | pure extension |
| Neutral | 16 | all have `M56_vd = 0` |
| Undefined | 31 | all have `M56_vd = -999` |

This matches Cerpa et al. 2025: **213** transects with constrained non-zero `vd` (260 − 16 − 31), ~90% of active zones, 9-class scheme (C / C-SS / SS-C / SS / SS-E / E-SS / E / Neutral / Undefined).

Zones in the dump: SEA 90, SAM 72, NPA 47, SWP 33, MED 18.

Conventions (SubMap database page): shortening ⇒ `0 ≤ Ovd < 90°`, `vdn > 0`; extension ⇒ `90° < Ovd ≤ 180°`, `vdn < 0`.

## Source-by-source log

### 1. Cerpa, Lallemand & Heuret 2025 EarthArXiv — PDF yes, no supplementary CSV/XLSX

- Preprint: https://eartharxiv.org/repository/view/9423/ DOI [10.31223/x56m8h](https://doi.org/10.31223/x56m8h)
- PDF: `cerpa_lallemand_heuret_2025_eartharxiv.pdf` (55 pp, 5.8 MB) from `https://eartharxiv.org/repository/object/9423/download/17567/`
- EarthArXiv hosts a **single PDF**; no linked Zenodo/Figshare/CSV. Data-availability statement in the PDF: *“Our database will be made available through our webtool submap (www.submap.fr).”*
- Tables in the PDF are **figures** (not selectable text):
  - `cerpa2025/table1_arc_blocks.png` — Table 1, 60 arc blocks with constrained `vd` (11 M56 + 35 rigid + 14 GPS-deformed areas). **No transect lon/lat or UPS class.**
  - `cerpa2025/table2_unconstrained_vd.png` — Table 2, 13 areas with unconstrained `vd`.
  - `cerpa2025/table1_ocr.txt` and `table2_ocr.txt` are noisy Tesseract dumps of those images; **do not treat OCR as authoritative**.
- Full extracted text: `cerpa2025/fulltext.txt`.

### 2. SubMap website — global table obtained via public API (see primary files)

- https://submap.gm.umontpellier.fr/ (SPA; JS bundle `/assets/index-*.js`, API base `https://submap.gm.umontpellier.fr/api`)
- Sub-DATA UI: `/data-index` and `/data-index/:short_name`
- How to cite (CMS, retrieved 2026-09-02): mention https://submap.fr plus the parameter-specific references. For M56+ kinematics / geodetic UPS cite Cerpa et al. 2025 (below). General SubMap parameters: Lallemand & Heuret (2017), doi:10.1016/B978-0-12-409548-9.09495-1. GMT maps: Wessel et al. (2019).
- Site is freely accessible; publication use requires those citations (no separate license text beyond “How to cite”).
- **SHP of UPS transects:** not found as a static public file. SHP download in the UI is for Subquake rupture shapes after a map generate, not the global UPS table.

### 3. Schellart 2024 Earth-Science Reviews — supplementary 258-segment table **not** publicly downloadable

- doi:10.1016/j.earscirev.2024.104755 (open-access article describing 28 zones / **258** ~200 km segments, including `v_OPD⊥`, age, slab width, etc.)
- Elsevier `mmc1`/`mmc2` URLs for pdf/xlsx/csv/docx all **404**. ScienceDirect lists supplementary **videos**, not a data spreadsheet.
- No VU/author CSV found. **Not downloaded.** Do not confuse this 258-segment compilation with SubMap’s 260 transects (different segmentation and UPS definition).

### 4. Heuret & Lallemand 2005 PEPI — PDF posted; Lallemand et al. 2005 G3 — no standalone supplement retrieved

**Heuret & Lallemand (2005) PEPI** doi:10.1016/j.pepi.2004.08.022  
Posted course PDF (Clint Conrad, UH): `other_sources/heuret_lallemand_2005_PEPI_clintconrad.pdf` (21 pp). Elsevier paywall still applies for the publisher version; this is a copy that is actually posted.

Table 1 (absolute motions, `Vdn`, **Df class** E3–C3, age) is in the PDF (journal pp. 36–41 ≈ PDF pp. 7–10):

- `other_sources/heuret_lallemand_2005_PEPI_Table1.txt` (pdftotext -layout)
- `other_sources/heuret_2005_table-0{7,8,9,10}.png` and `heuret_2005_table_p*.txt`

Two-column layout prevented a clean lossless CSV without row-merging errors, so **no machine CSV was kept**. Df class is the **2005 focal-mechanism UPS**, superseded by SubMap geodetic UPS.

**Lallemand, Heuret & Boutelier (2005) G3** doi:10.1029/2005GC000917  
Table 1 (159 transects, columns include lat/lon/UPS) is **inside the article**, not a separate supplementary CSV. HAL (`hal-01261567`) returned a bot-check HTML page; AGU `pdfdirect` was not retrieved here. **No G3 file saved.**

## Recommended citations

Cerpa, N. G., Lallemand, S., & Heuret, A. (2025). Unraveling the effect of convergence obliquity on overriding-plate deformation and strain partitioning in subduction zones. EarthArXiv. https://doi.org/10.31223/X56M8H (submitted to *Tektonika*; SubMap v7 M56+ kinematics / geodetic UPS).

SubMap web-tool: https://submap.gm.umontpellier.fr/ (https://submap.fr). Lallemand, S., & Heuret, A. (2017). Subduction Zones Parameters. In *Reference Module in Earth Systems and Environmental Sciences*. Elsevier. https://doi.org/10.1016/B978-0-12-409548-9.09495-1

Heuret, A., & Lallemand, S. (2005). Plate motions, slab dynamics and back-arc deformation. *PEPI*, 149, 31–51. https://doi.org/10.1016/j.pepi.2004.08.022 (historical Df-class table only).
