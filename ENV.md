# pyGPlates / GPlately environment (subduction-test)

Installed 2026-09-02 ~12:53 MDT on this Linux box. System Python is 3.13.5; the stack lives in a dedicated micromamba env with Python 3.11.

## Activate

Micromamba binary: `/workspace/subduction-test/.micromamba/bin/micromamba` (version 2.9.0)  
Root prefix: `/workspace/subduction-test/.micromamba`  
Env name: `pygplates`  
Env prefix: `/workspace/subduction-test/.micromamba/envs/pygplates`  
Python: `/workspace/subduction-test/.micromamba/envs/pygplates/bin/python`

Interactive:

```bash
export MAMBA_ROOT_PREFIX=/workspace/subduction-test/.micromamba
eval "$(/workspace/subduction-test/.micromamba/bin/micromamba shell hook -s bash)"
micromamba activate pygplates
```

One-shot (no shell init):

```bash
export MAMBA_ROOT_PREFIX=/workspace/subduction-test/.micromamba
/workspace/subduction-test/.micromamba/bin/micromamba run -n pygplates --root-prefix "$MAMBA_ROOT_PREFIX" python
```

Direct interpreter:

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python
```

Import check:

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python -c \
  "import pygplates, gplately; print(pygplates.__version__ if hasattr(pygplates,'__version__') else 'ok', gplately.__version__)"
```

Output: `1.0.0 2.0.0`

## Versions

| package | version |
|---|---|
| python | 3.11.16 (conda-forge, GCC 14.4.0) |
| pygplates | 1.0.0 (`py311h795f763_7`, conda-forge) |
| gplately | 2.0.0 |
| plate-model-manager | 1.3.2 |
| cartopy | 0.25.0 |
| geopandas | 1.1.4 |
| pandas | 3.0.5 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| matplotlib | 3.11.1 |
| shapely | 2.1.2 |
| rasterio | 1.4.4 |
| pyproj | 3.7.2 |
| xarray | 2026.7.0 |
| netCDF4 | 1.7.4 |

Channel: conda-forge. Env size ~2.9G.

## Data: Müller et al. (2019)

GPlately name: `Muller2019` (global deforming plate model, 0–250 Ma, mantle reference frame).  
Downloaded with `gplately.PlateModelManager().get_model("Muller2019", data_dir="/workspace/subduction-test/data")`.  
DOI / version in model metadata: https://doi.org/10.5281/zenodo.10525286 , `10.5281/zenodo.11601026`.

Canonical directory (requested path):

`/workspace/subduction-test/data/Muller2019`

PlateModelManager looks for lowercase `muller2019`; a symlink is in place:

`/workspace/subduction-test/data/muller2019` → `Muller2019`

Top-level listing:

```
ContinentalPolygons/
Rasters/
Rotations/
StaticPolygons/
Topologies/
readme.txt
.metadata.json
```

Key files:

| role | path |
|---|---|
| rotations | `Rotations/Muller_etal_2019_CombinedRotations.rot` |
| topologies (incl. deforming networks) | `Topologies/Muller_etal_2019_PlateBoundaries_DeformingNetworks.gpmlz` |
| static polygons | `StaticPolygons/Muller_etal_2019_Global_StaticPlatePolygons.gpmlz` |
| continental polygons | `ContinentalPolygons/Global_PresentDay_ContPolygons_2019_v1.shp` (+ `.dbf/.shx/.prj/.xml`) |
| present-day seafloor age grid (0 Ma) | `Rasters/AgeGrids/Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc` (~8 MB, raster shape 1801×3601) |

Only the 0 Ma age grid was downloaded (not the full 0–250 Ma stack). Other PMM layers (coastlines, COBs, hotspots, LIPs, sediment thickness) were not fetched.

## Smoke test (0 Ma tessellated subduction zones)

Command:

```bash
export MAMBA_ROOT_PREFIX=/workspace/subduction-test/.micromamba
/workspace/subduction-test/.micromamba/bin/micromamba run -n pygplates --root-prefix "$MAMBA_ROOT_PREFIX" python << 'PY'
import numpy as np
import gplately
from gplately import PlateModelManager

m = PlateModelManager().get_model("Muller2019", data_dir="/workspace/subduction-test/data")
recon = gplately.PlateReconstruction(
    m.get_rotation_model(), m.get_topologies(), m.get_static_polygons()
)
arr = np.asarray(
    recon.tessellate_subduction_zones(
        0, tessellation_threshold_radians=np.radians(2.0), ignore_warnings=True
    )
)
print("shape:", arr.shape)
print(arr[:5])
print("age:", m.get_age_grid(0))
PY
```

Output (2026-09-02):

```
shape: (335, 10)
[[   144.3749     41.1318      9.6717     31.1002      4.8435   -107.6549      1.        276.6068    901.     601118.    ]
 [   144.0718     39.1453      9.4947     18.98        3.6431   -124.6006      2.        278.8283    901.     601116.    ]
 [   143.4293     37.2294      9.3577      5.3753      3.7504   -133.4061      2.        294.5663    901.     601114.    ]
 [   142.2342     35.4825      8.738     -14.3462      3.3397   -144.6334      1.3156    314.2349    901.     601112.    ]
 [   140.0469     34.7519      6.0964    -57.9722      1.9373    127.981       2.6844      3.9214    659.     601127.    ]]
age: /workspace/subduction-test/data/muller2019/Rasters/AgeGrids/Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc
```

(`muller2019` in the printed age path is the symlink to `Muller2019`.)

Column meanings for `tessellate_subduction_zones` (10-col float64 array):

| col | meaning |
|---|---|
| 0 | longitude of sampled trench point (deg) |
| 1 | latitude of sampled trench point (deg) |
| 2 | subducting convergence (relative to trench) velocity magnitude (cm/yr) |
| 3 | subducting convergence velocity obliquity angle (deg; trench-normal vs convergence) |
| 4 | trench absolute (relative to anchor plate) velocity magnitude (cm/yr) |
| 5 | trench absolute velocity obliquity angle (deg) |
| 6 | length of arc segment that current point is on (deg) |
| 7 | trench normal azimuth (clockwise from North, 0–360 deg) |
| 8 | subducting plate ID |
| 9 | trench plate ID |

Smoke-test tessellation used a coarse 2° threshold so only a few hundred points are returned. Default `tessellation_threshold_radians=0.001` is much finer.

## Blockers

None. pygplates imported and reconstructed at 0 Ma successfully.
