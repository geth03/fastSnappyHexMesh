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

  fastSnappyHexMesh - see fsThreads.H
\*---------------------------------------------------------------------------*/

#include "fsThreads.H"
#include "triSurface_searchableSurface.H"
#include "triSurfaceSearch.H"
#include "indexedOctree.H"
#include "treeDataTriSurface.H"
#include "Pstream.H"

#include <thread>
#include <vector>
#include <cstdlib>
#include <typeinfo>

// * * * * * * * * * * * * * * * Functions  * * * * * * * * * * * * * * * * //

Foam::label Foam::fsThreads::nThreads()
{
    static label n = -1;

    if (n < 0)
    {
        if (const char* e = std::getenv("FS_THREADS"))
        {
            n = max(1, std::atoi(e));
        }
        else
        {
            const label hw = max(1, label(std::thread::hardware_concurrency()));
            const label physical = max(1, hw/2);
            n = max(1, physical/max(1, Pstream::nProcs()));
        }
    }

    return n;
}


void Foam::fsThreads::findLineAny
(
    const searchableSurface& geom,
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
)
{
    const label nRays = start.size();
    const label nt = min(nThreads(), max(1, nRays/2000));

    if (nt <= 1 || typeid(geom) != typeid(searchableSurfaces::triSurface))
    {
        geom.findLineAny(start, end, info);
        return;
    }

    const triSurfaceSearch& search = dynamic_cast<const triSurfaceSearch&>
    (
        dynamic_cast<const searchableSurfaces::triSurface&>(geom)
    );

    // As triSurfaceSearch::findLineAny: build tree, set tolerance (static)
    // in this thread; the workers only read them.
    const indexedOctree<treeDataTriSurface>& octree = search.tree();

    info.setSize(nRays);

    const scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = search.tolerance();

    const label chunk = (nRays + nt - 1)/nt;

    std::vector<std::thread> workers;
    workers.reserve(nt);
    for (label t = 0; t < nt; t++)
    {
        const label i0 = t*chunk;
        const label i1 = min(nRays, i0 + chunk);
        if (i0 >= i1)
        {
            break;
        }
        workers.emplace_back
        (
            [&octree, &start, &end, &info, i0, i1]()
            {
                for (label i = i0; i < i1; i++)
                {
                    info[i] = octree.findLineAny(start[i], end[i]);
                }
            }
        );
    }
    for (std::thread& w : workers)
    {
        w.join();
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::fsThreads::findLine
(
    const searchableSurface& geom,
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
)
{
    const label nRays = start.size();

    if
    (
        min(nThreads(), max(1, nRays/2000)) <= 1
     || typeid(geom) != typeid(searchableSurfaces::triSurface)
    )
    {
        geom.findLine(start, end, info);
        return;
    }

    // As triSurfaceSearch::findLine
    const triSurfaceSearch& search = dynamic_cast<const triSurfaceSearch&>
    (
        dynamic_cast<const searchableSurfaces::triSurface&>(geom)
    );

    const indexedOctree<treeDataTriSurface>& octree = search.tree();

    info.setSize(nRays);

    const scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = search.tolerance();

    forRange
    (
        nRays,
        [&](const label i0, const label i1)
        {
            for (label i = i0; i < i1; i++)
            {
                info[i] = octree.findLine(start[i], end[i]);
            }
        },
        2000
    );

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::fsThreads::findNearest
(
    const searchableSurface& geom,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
)
{
    const label nSamples = samples.size();

    if
    (
        min(nThreads(), max(1, nSamples/2000)) <= 1
     || typeid(geom) != typeid(searchableSurfaces::triSurface)
    )
    {
        geom.findNearest(samples, nearestDistSqr, info);
        return;
    }

    // As triSurfaceSearch::findNearest
    const triSurfaceSearch& search = dynamic_cast<const triSurfaceSearch&>
    (
        dynamic_cast<const searchableSurfaces::triSurface&>(geom)
    );

    const scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = search.tolerance();

    const indexedOctree<treeDataTriSurface>& octree = search.tree();

    info.setSize(nSamples);

    forRange
    (
        nSamples,
        [&](const label i0, const label i1)
        {
            for (label i = i0; i < i1; i++)
            {
                info[i] = octree.findNearest
                (
                    samples[i],
                    nearestDistSqr[i],
                    treeDataTriSurface::findNearestOp(octree)
                );
            }
        },
        2000
    );

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::fsThreads::findNearest
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearestSurfaces,
    List<pointIndexHit>& nearestInfo
)
{
    // Find nearest. Return -1 or nearest point

    // Initialise
    nearestSurfaces.setSize(samples.size());
    nearestSurfaces = -1;
    nearestInfo.setSize(samples.size());

    // Work arrays
    scalarField minDistSqr(nearestDistSqr);
    List<pointIndexHit> hitInfo(samples.size());

    forAll(surfacesToTest, testI)
    {
        findNearest
        (
            allSurfaces[surfacesToTest[testI]],
            samples,
            minDistSqr,
            hitInfo
        );

        // Update minDistSqr and arguments
        forAll(hitInfo, pointi)
        {
            if (hitInfo[pointi].hit())
            {
                minDistSqr[pointi] = magSqr
                (
                    hitInfo[pointi].hitPoint()
                  - samples[pointi]
                );
                nearestInfo[pointi] = hitInfo[pointi];
                nearestSurfaces[pointi] = testI;
            }
        }
    }
}


// ************************************************************************* //
