# Competing kinematic controls on upper-plate stress

PyGPlates / GPlately test of Heuret–Lallemand vs Schellart hypotheses for compressive, neutral, and extensional upper-plate regimes at subduction zones.

Investigator: Grok, Tectonics Researcher. Advisor: Andrew Laskowski.

## Status

Present-day (0 Ma) milestone is in. Cenozoic trench catalog (0–60 Ma, 1 Myr) is extracted. Müller et al. (2019) **cannot** separate trench rollback from overriding-plate motion (`v_T ≈ −v_OP`); use SubMap M56 `Vupn1` / `Vtn1` for that 0 Ma H1 vs H2 test.

## Layout

- `docs/subduction-dynamics-proposal.md` — approved proposal
- `scripts/` — extraction, SubMap join, H1/H2 stats
- `results/` — CSVs, figures, notes
- `data/submap/` — SubMap v7 geodetic UPS table (cite Cerpa, Lallemand & Heuret 2025)
- `data/Muller2019/` — rotations, topologies, polygons (age grids are gitignored)

## Setup

Python 3.11 env with pygplates 1.0 and gplately 2.0 (conda-forge). See `ENV.md`.

```bash
# Müller2019 age grids (0–60 Ma) if you need the Cenozoic extract
# PlateModelManager downloads them into data/Muller2019/Rasters/AgeGrids/
```

## Key 0 Ma results (short)

- 241 / 260 SubMap transects join to Müller2019 trenches (median 25 km).
- SubMap M56, C vs E (n=79): H1 `Vupn1` ρ = 0.50; H2 `Vtn1` ρ = −0.34; slab age median 38 Ma (C) vs 94 Ma (E). Not a knockout winner.
- Width and distance-to-edge do not separate C from E.

## Cite

- Müller et al. (2019), *Tectonics* https://doi.org/10.1029/2018TC005462
- SubMap: https://submap.fr — kinematics / geodetic UPS: Cerpa, Lallemand & Heuret (2025) https://doi.org/10.31223/x56m8h
- Heuret & Lallemand (2005), Lallemand et al. (2005), Schellart (2008, 2024)
