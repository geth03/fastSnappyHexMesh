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

  fastSnappyHexMesh - see fsPatchDist.H and fsPointEdgeWave.H
\*---------------------------------------------------------------------------*/

#include "fsPatchDist.H"
#include "fsPointEdgeWave.H"
#include "fsProfiler.H"
#include "pointEdgePoint.H"

#include <cstdlib>
#include <string>

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::fsWave::usable(const polyMesh& mesh)
{
    {
        const char* e = std::getenv("FS_FASTWAVE");
        if (e && std::string(e) == "0")
        {
            return false;
        }
    }

    if (Pstream::parRun() || fsThreads::nThreads() <= 1)
    {
        return false;
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchi)
    {
        if (patches[patchi].coupled())
        {
            return false;
        }
    }

    return true;
}


bool Foam::fsPatchDist::edgePatchDist
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    scalarField& edgeDist
)
{
    if (!fsWave::usable(mesh))
    {
        return false;
    }

    FS_PROF("fsPatchDist");

    // As snappySnapDriver::edgePatchDist
    List<pointEdgePoint> wallInfo(pp.nPoints());

    forAll(pp.localPoints(), ppi)
    {
        wallInfo[ppi] = pointEdgePoint(pp.localPoints()[ppi], 0.0);
    }

    List<pointEdgePoint> allPointInfo(mesh.nPoints());
    List<pointEdgePoint> allEdgeInfo(mesh.nEdges());

    const label maxIter = mesh.globalData().nTotalPoints();
    int dummyTrackData = 0;

    const label iter = fsWave::pointEdgeWave
    (
        mesh,
        pp.meshPoints(),
        wallInfo,
        allPointInfo,
        allEdgeInfo,
        maxIter,
        dummyTrackData
    );

    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter." << endl
            << "    maxIter:" << maxIter << endl
            << exit(FatalError);
    }

    edgeDist.setSize(mesh.nEdges());
    fsThreads::forRange
    (
        mesh.nEdges(),
        [&](const label e0, const label e1)
        {
            for (label edgei = e0; edgei < e1; edgei++)
            {
                edgeDist[edgei] = Foam::sqrt(allEdgeInfo[edgei].distSqr());
            }
        }
    );

    return true;
}


// ************************************************************************* //
