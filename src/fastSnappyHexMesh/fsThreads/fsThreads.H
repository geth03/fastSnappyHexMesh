/*---------------------------------------------------------------------------*\
  fastSnappyHexMesh
  Copyright (C) 2026 Seyit Temiz
-------------------------------------------------------------------------------
License
    This file is part of fastSnappyHexMesh.

    fastSnappyHexMesh is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    fastSnappyHexMesh is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fastSnappyHexMesh.  If not, see
    <http://www.gnu.org/licenses/>.

Namespace
    Foam::fsThreads

Description
    Shared-memory threading for independent per-item loops (rays).

    nThreads(): env FS_THREADS if set, else the physical cores left per MPI
    rank (hardware threads / 2 / nProcs, at least 1). With one rank per core
    this is 1, i.e. no threads compete with the MPI ranks.

    findLineAny(): threaded version of searchableSurface::findLineAny for
    plain triangulated surfaces (searchableSurfaces::triSurface). Every ray
    is an independent const octree query; the octree is built and the
    search tolerance set in the calling thread, exactly as
    triSurfaceSearch::findLineAny does. Results are identical to the
    serial call. Other surface types (analytic, distributed) fall back to
    their own findLineAny.

SourceFiles
    fsThreads.C

\*---------------------------------------------------------------------------*/

#ifndef fsThreads_H
#define fsThreads_H

#include "searchableSurface.H"
#include "pointIndexHit.H"

#include <thread>
#include <vector>

namespace Foam
{
namespace fsThreads
{

//- Number of threads to use per rank
label nThreads();

//- Run f(i0, i1) over [0, n) in contiguous chunks, at most nThreads()
//  threads with at least minChunk items each; inline when that is one
template<class F>
void forRange(const label n, const F& f, const label minChunk = 10000)
{
    const label nt = min(nThreads(), max(label(1), n/minChunk));

    if (nt <= 1)
    {
        f(label(0), n);
        return;
    }

    const label chunk = (n + nt - 1)/nt;
    std::vector<std::thread> workers;
    for (label t = 0; t < nt; t++)
    {
        const label i0 = t*chunk;
        const label i1 = min(n, i0 + chunk);
        if (i0 >= i1)
        {
            break;
        }
        workers.emplace_back([&f, i0, i1]() { f(i0, i1); });
    }
    for (std::thread& w : workers)
    {
        w.join();
    }
}

//- Number of chunks for n items with at least minChunk items each
inline label nChunks(const label n, const label minChunk)
{
    return min(nThreads(), max(label(1), n/minChunk));
}

//- Run f(t, i0, i1) on the nt contiguous chunks of [0, n), chunk t on its
//  own thread (inline when nt is one). For results combined per chunk.
template<class F>
void forChunks(const label n, const label nt, const F& f)
{
    if (nt <= 1)
    {
        f(label(0), label(0), n);
        return;
    }

    const label chunk = (n + nt - 1)/nt;
    std::vector<std::thread> workers;
    for (label t = 0; t < nt; t++)
    {
        const label i0 = min(n, t*chunk);
        const label i1 = min(n, i0 + chunk);
        workers.emplace_back([&f, t, i0, i1]() { f(t, i0, i1); });
    }
    for (std::thread& w : workers)
    {
        w.join();
    }
}

//- See description
void findLineAny
(
    const searchableSurface& geom,
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
);

//- Threaded searchableSurface::findLine, as findLineAny
void findLine
(
    const searchableSurface& geom,
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
);

//- Threaded searchableSurface::findNearest, as findLineAny
void findNearest
(
    const searchableSurface& geom,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
);

//- searchableSurfacesQueries::findNearest on the threaded findNearest
void findNearest
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearestSurfaces,
    List<pointIndexHit>& nearestInfo
);

} // End namespace fsThreads
} // End namespace Foam

#endif

// ************************************************************************* //
