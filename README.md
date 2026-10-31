# fastSnappyHexMesh

A faster `snappyHexMesh` for OpenFOAM-13 (openfoam.org).

It reads the same `system/snappyHexMeshDict`, takes the same options and, in
serial, writes **the same mesh, byte for byte**, only about **3x faster** and
with **about 20 % less memory**. Castellation, snapping and layer addition
are all faster.

## Speed

Stock `snappyHexMesh` vs `fastSnappyHexMesh`, serial, same machine
(Intel i7-10850H laptop, 6 cores), identical meshes:

| case | snappyHexMesh | fastSnappyHexMesh | speed-up | peak memory |
|---|---:|---:|---:|---:|
| motorBike tutorial, full run | 707 s | 196 s | **3.6x** | 830 → 673 MB |
| Ahmed body, full run (2.8M cells) | 676 s | 244 s | **2.8x** | 4.5 → 3.7 GB |
| Ahmed body, 14 layers for y+ = 1 (7.1M cells) | 518 s | 135 s | **3.8x** | 7.5 → 5.1 GB |
| Ahmed body, castellation only (3.5M cells) | 412 s | 130 s | **3.1x** | 4.8 → 4.0 GB |

In serial, fastSnappyHexMesh uses all physical cores through threads.
In parallel (MPI) runs the gain is smaller (about 1.4-2x on 4 ranks).

## What is different

- Faster algorithms for cell search, region splitting, mesh topology changes
  and mesh quality checks
- Threads for ray casting, mesh geometry, displacement smoothing, wall
  distance, nearest-surface searches and the mesh quality checks
- Less copying when the mesh changes, old data freed earlier
- The tcmalloc memory allocator is linked into the executable
  (no extra installation needed to run it)

Two changes are intended to alter **parallel** results slightly:

- a fix for a refinement inconsistency at processor boundaries (refinement
  could creep on for many extra iterations and add cells), and
- no full rebalance for a refinement step that adds less than 5 % cells.

Set `FS_CREEPFIX=0 FS_BALRULE=0` to switch both off; the parallel mesh is then
the same as the one of stock snappyHexMesh as well.

## Requirements

- OpenFOAM-13 from openfoam.org (Linux, or WSL on Windows)
- Optional, recommended: `libgoogle-perftools-dev` (only needed to build;
  about 15 % faster)

```sh
sudo apt install libgoogle-perftools-dev
```

## Build

```sh
. /opt/openfoam13/etc/bashrc
git clone https://github.com/geth03/fastSnappyHexMesh.git
cd fastSnappyHexMesh
./Allwmake
```

This installs `fastSnappyHexMesh` and `libfastSnappyHexMesh.so` into your
OpenFOAM user directories. The stock snappyHexMesh is not touched.

## Use

Exactly like snappyHexMesh, in any existing case:

```sh
blockMesh
surfaceFeatures                       # if the case uses feature edges
fastSnappyHexMesh
```

Parallel:

```sh
decomposePar
mpirun -np 4 fastSnappyHexMesh -parallel
```

`system/decomposeParDict` needs both entries:

```
decomposer      scotch;     // initial split, used by decomposePar
distributor     ptscotch;   // load balancing during refinement
```

`fastSnappyHexMesh -help` lists everything, including recommended settings
for thin boundary layers.

### Example

```sh
cd tutorials/ahmedBody
./Allrun          # mesh with 3 boundary layers, open ahmedBody.foam in ParaView
./Allclean
```

### Options (environment variables)

| variable | effect |
|---|---|
| `FS_THREADS=n` | threads per process (default: physical cores / MPI ranks) |
| `FS_LOWMEM=1` | lower peak memory (about -28 %), slower (about +26 %) |
| `FS_PROFILE=0` | no timing profile (default: written to `fsProfile/` in the case) |
| `FS_CREEPFIX=0`, `FS_BALRULE=0` | parallel: same mesh as stock snappyHexMesh (see above) |

In parallel, pass them with `mpirun -x`, e.g.
`mpirun -np 2 -x FS_THREADS=3 fastSnappyHexMesh -parallel`.

## Thin boundary layers (y+ = 1)

Some of the default mesh quality checks judge the cell *shape*, not its
validity, and reject very thin layer cells (`minDeterminant 0.001` rejects
cells flatter than about 1:165). Keep the usual values for snapping and let
layer addition use validity checks only:

```
addLayersControls
{
    ...
    nRelaxedIter    0;      // layers use 'relaxed' from the start
}

meshQualityControls
{
    ...                     // snapping: usual values
    relaxed                 // layers: reject invalid cells only
    {
        maxNonOrtho         85;
        maxBoundarySkewness 40;
        maxInternalSkewness 20;
        minVol              1e-30;
        minTwist            0.01;
        minDeterminant      1e-30;
        minFaceWeight       0.01;
        minVolRatio         0.001;
    }
}
```

Ahmed body with 14 layers from 0.0095 mm: layers on 14.7 % of the faces with
the default settings, 99.9 % with these.

## License

GNU General Public License v3 or later, see [LICENSE](LICENSE).

fastSnappyHexMesh is a modified version of snappyHexMesh from OpenFOAM-13,
Copyright (C) 2011-2025 OpenFOAM Foundation. Changes and new code
Copyright (C) 2026 Seyit Temiz.

This offering is not approved or endorsed by OpenCFD Limited, producer and
distributor of the OpenFOAM software via www.openfoam.com, and owner of the
OPENFOAM® and OpenCFD® trade marks.
