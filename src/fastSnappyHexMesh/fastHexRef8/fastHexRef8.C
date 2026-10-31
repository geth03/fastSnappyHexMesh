/*---------------------------------------------------------------------------*\
  fastSnappyHexMesh
  Copyright (C) 2026 Seyit Temiz
  Based on OpenFOAM: Copyright (C) 2011-2025 OpenFOAM Foundation
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

  fastSnappyHexMesh - see fastHexRef8.H
  Based on hexRef8.C of OpenFOAM-13 (GPL v3).
\*---------------------------------------------------------------------------*/

#include "fastHexRef8.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "FaceCellWave.H"
#include "refinementDistanceData.H"

#include <cstdlib>
#include <cstring>

// * * * * * * * * * * * * * * * Functions  * * * * * * * * * * * * * * * * //

bool Foam::fastHexRef8::creepFixEnabled()
{
    const char* e = std::getenv("FS_CREEPFIX");
    return !(e && std::strcmp(e, "0") == 0);
}


Foam::label Foam::fastHexRef8::faceConsistentRefinement
(
    const fsHexRef8& cutter,
    const bool maxSet,
    PackedBoolList& refineCell
)
{
    const polyMesh& mesh = cutter.mesh();
    const labelList& cellLevel = cutter.cellLevel();

    label nChanged = 0;

    // Internal faces.
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        label own = mesh.faceOwner()[facei];
        label ownLevel = cellLevel[own] + refineCell.get(own);

        label nei = mesh.faceNeighbour()[facei];
        label neiLevel = cellLevel[nei] + refineCell.get(nei);

        if (ownLevel > (neiLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(nei);
            }
            else
            {
                refineCell.unset(own);
            }
            nChanged++;
        }
        else if (neiLevel > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
            }
            else
            {
                refineCell.unset(nei);
            }
            nChanged++;
        }
    }


    // Coupled faces. Swap owner level to get neighbouring cell level.
    labelList neiLevel(mesh.nFaces()-mesh.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh.faceOwner()[i+mesh.nInternalFaces()];

        neiLevel[i] = cellLevel[own] + refineCell.get(own);
    }

    syncTools::swapBoundaryFaceList(mesh, neiLevel);

    forAll(neiLevel, i)
    {
        label own = mesh.faceOwner()[i+mesh.nInternalFaces()];
        label ownLevel = cellLevel[own] + refineCell.get(own);

        if (ownLevel > (neiLevel[i]+1))
        {
            if (!maxSet)
            {
                refineCell.unset(own);
                nChanged++;
            }
        }
        else if (neiLevel[i] > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
                nChanged++;
            }
        }
    }

    return nChanged;
}


Foam::labelList Foam::fastHexRef8::consistentSlowRefinement2
(
    const fsHexRef8& cutter,
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck
)
{
    const polyMesh& mesh = cutter.mesh();
    const labelList& cellLevel = cutter.cellLevel();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const label nInt = mesh.nInternalFaces();
    const label nBnd = mesh.nFaces() - nInt;

    if (maxFaceDiff <= 0)
    {
        FatalErrorInFunction
            << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }

    const scalar level0Size = 2*maxFaceDiff*cutter.level0EdgeLength();

    List<refinementDistanceData> allCellInfo(mesh.nCells());
    List<refinementDistanceData> allFaceInfo(mesh.nFaces());
    int dummyTrackData = 0;

    const pointField& cc = mesh.cellCentres();

    // Mark cells with wanted refinement level
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli] = refinementDistanceData
        (
            level0Size,
            cc[celli],
            cellLevel[celli]+1             // wanted refinement
        );
    }
    // Mark all others with existing refinement level
    forAll(allCellInfo, celli)
    {
        if (!allCellInfo[celli].valid(dummyTrackData))
        {
            allCellInfo[celli] = refinementDistanceData
            (
                level0Size,
                cc[celli],
                cellLevel[celli]           // wanted refinement
            );
        }
    }

    // fastSnappy: level and cell centre of the cell across each coupled face
    boolList isCoupledBFace(nBnd, false);
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];
            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    isCoupledBFace[pp.start() + i - nInt] = true;
                }
            }
        }
    }

    labelList cellOriginLevel(mesh.nCells());
    forAll(allCellInfo, celli)
    {
        cellOriginLevel[celli] = allCellInfo[celli].originLevel();
    }
    labelList bNeiLevel;
    syncTools::swapBoundaryCellList(mesh, cellOriginLevel, bNeiLevel);
    pointField bNeiCc;
    syncTools::swapBoundaryCellPositions(mesh, cc, bNeiCc);


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementDistanceData> seedFacesInfo(mesh.nFaces()/100);

    const scalar tol =
        FaceCellWave<refinementDistanceData, int>::propagationTol();

    forAll(facesToCheck, i)
    {
        label facei = facesToCheck[i];

        if (allFaceInfo[facei].valid(dummyTrackData))
        {
            FatalErrorInFunction
                << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << facei << " occurs at positions "
                << findIndices(facesToCheck, facei)
                << abort(FatalError);
        }

        label own = faceOwner[facei];

        label ownLevel =
        (
            allCellInfo[own].valid(dummyTrackData)
          ? allCellInfo[own].originLevel()
          : cellLevel[own]
        );

        const bool coupled =
            !mesh.isInternalFace(facei) && isCoupledBFace[facei - nInt];

        if (!mesh.isInternalFace(facei) && !coupled)
        {
            // Real boundary: do as if it had a neighbour with one higher
            // refinement level.
            const point& fc = mesh.faceCentres()[facei];

            refinementDistanceData neiData
            (
                level0Size,
                2*fc - cc[own],    // est'd cell centre
                ownLevel+1
            );

            allFaceInfo[facei].updateFace
            (
                mesh, facei, own, neiData, tol, dummyTrackData
            );
        }
        else
        {
            // Internal face, or coupled face treated like one
            label neiLevel;
            point neiCc;
            if (coupled)
            {
                neiLevel = bNeiLevel[facei - nInt];
                neiCc = bNeiCc[facei - nInt];
            }
            else
            {
                label nei = faceNeighbour[facei];
                neiLevel =
                (
                    allCellInfo[nei].valid(dummyTrackData)
                  ? allCellInfo[nei].originLevel()
                  : cellLevel[nei]
                );
                neiCc = cc[nei];
            }

            if (ownLevel == neiLevel)
            {
                // Fake as if nei>own or own>nei (whichever one 'wins')
                allFaceInfo[facei].updateFace
                (
                    mesh, facei, own,
                    refinementDistanceData(level0Size, neiCc, neiLevel+1),
                    tol, dummyTrackData
                );
                allFaceInfo[facei].updateFace
                (
                    mesh, facei, own,
                    refinementDistanceData(level0Size, cc[own], ownLevel+1),
                    tol, dummyTrackData
                );
            }
            else
            {
                // Difference in level anyway.
                allFaceInfo[facei].updateFace
                (
                    mesh, facei, own,
                    refinementDistanceData(level0Size, neiCc, neiLevel),
                    tol, dummyTrackData
                );
                allFaceInfo[facei].updateFace
                (
                    mesh, facei, own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    tol, dummyTrackData
                );
            }
        }
        seedFaces.append(facei);
        seedFacesInfo.append(allFaceInfo[facei]);
    }


    // Seed all faces in between different refinement levels: internal faces
    // (as the original) and coupled faces (fastSnappy).
    for (label facei = 0; facei < mesh.nFaces(); facei++)
    {
        const bool internal = facei < nInt;
        if (!internal && !isCoupledBFace[facei - nInt])
        {
            continue;
        }

        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];

            label ownLevel =
            (
                allCellInfo[own].valid(dummyTrackData)
              ? allCellInfo[own].originLevel()
              : cellLevel[own]
            );

            label neiLevel;
            point neiCc;
            if (internal)
            {
                label nei = faceNeighbour[facei];
                neiLevel =
                (
                    allCellInfo[nei].valid(dummyTrackData)
                  ? allCellInfo[nei].originLevel()
                  : cellLevel[nei]
                );
                neiCc = cc[nei];
            }
            else
            {
                neiLevel = bNeiLevel[facei - nInt];
                neiCc = bNeiCc[facei - nInt];
            }

            if (ownLevel > neiLevel)
            {
                seedFaces.append(facei);
                allFaceInfo[facei].updateFace
                (
                    mesh, facei, own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    tol, dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[facei]);
            }
            else if (neiLevel > ownLevel)
            {
                seedFaces.append(facei);
                allFaceInfo[facei].updateFace
                (
                    mesh, facei, own,
                    refinementDistanceData(level0Size, neiCc, neiLevel),
                    tol, dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[facei]);
            }
        }
    }

    seedFaces.shrink();
    seedFacesInfo.shrink();

    // face-cell-face transport engine
    FaceCellWave<refinementDistanceData, int> levelCalc
    (
        mesh,
        seedFaces,
        seedFacesInfo,
        allFaceInfo,
        allCellInfo,
        mesh.globalData().nTotalCells()+1,
        dummyTrackData
    );


    // Convert back to labelList of cells to refine.

    // 1. Force original refinement cells to be picked up
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli].originLevel() = sizeof(label)*8-2;
        allCellInfo[celli].origin() = cc[celli];
    }

    // 2. Extend to 2:1
    PackedBoolList refineCell(mesh.nCells());
    forAll(allCellInfo, celli)
    {
        label wanted = allCellInfo[celli].wantedLevel(cc[celli]);

        if (wanted > cellLevel[celli]+1)
        {
            refineCell.set(celli);
        }
    }
    faceConsistentRefinement(cutter, true, refineCell);

    while (true)
    {
        label nChanged = faceConsistentRefinement(cutter, true, refineCell);

        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }
    }

    // 3. Convert back to labelList.
    label nRefined = 0;
    forAll(refineCell, celli)
    {
        if (refineCell[celli])
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;
    forAll(refineCell, celli)
    {
        if (refineCell[celli])
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    return newCellsToRefine;
}


// ************************************************************************* //
