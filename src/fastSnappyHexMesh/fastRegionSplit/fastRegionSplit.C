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

  fastSnappyHexMesh - see fastRegionSplit.H
\*---------------------------------------------------------------------------*/

#include "fastRegionSplit.H"
#include "syncTools.H"
#include "Map.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fastRegionSplit::propagateMin
(
    const boolList& blockedFace,
    labelList& regionLabel
) const
{
    const labelList& cellRegion = *this;
    const labelList& own = mesh_.faceOwner();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label nInt = mesh_.nInternalFaces();

    // Unblocked coupled faces
    DynamicList<label> cFaces;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (pp.coupled())
        {
            forAll(pp, i)
            {
                const label facei = pp.start() + i;
                if (blockedFace.empty() || !blockedFace[facei])
                {
                    cFaces.append(facei);
                }
            }
        }
    }

    if (!returnReduce(cFaces.size() > 0, orOp<bool>()))
    {
        return;
    }

    labelList bv(mesh_.nFaces() - nInt);

    while (true)
    {
        bv = labelMax;
        forAll(cFaces, i)
        {
            const label facei = cFaces[i];
            bv[facei - nInt] = regionLabel[cellRegion[own[facei]]];
        }

        syncTools::swapBoundaryFaceList(mesh_, bv);

        bool changed = false;
        forAll(cFaces, i)
        {
            const label facei = cFaces[i];
            const label nbr = bv[facei - nInt];
            label& mine = regionLabel[cellRegion[own[facei]]];
            if (nbr < mine)
            {
                mine = nbr;
                changed = true;
            }
        }

        if (!returnReduce(changed, orOp<bool>()))
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fastRegionSplit::fastRegionSplit
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections
)
:
    labelList(mesh.nCells(), -1),
    mesh_(mesh)
{
    labelList& cellRegion = *this;

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const cellList& cells = mesh.cells();
    const label nInt = mesh.nInternalFaces();

    // Explicit face-face connections (baffle pairs). As in FaceCellWave the
    // connection only carries information if both faces are unblocked.
    Map<label> partner(2*explicitConnections.size());
    forAll(explicitConnections, i)
    {
        const label a = explicitConnections[i].first();
        const label b = explicitConnections[i].second();
        const bool aOpen = blockedFace.empty() || !blockedFace[a];
        const bool bOpen = blockedFace.empty() || !blockedFace[b];
        if (aOpen && bOpen)
        {
            partner.insert(a, b);
            partner.insert(b, a);
        }
    }


    // 1. Local flood fill
    // ~~~~~~~~~~~~~~~~~~~

    label nLocal = 0;
    DynamicList<label> stack(1024);

    forAll(cellRegion, seedi)
    {
        if (cellRegion[seedi] != -1)
        {
            continue;
        }

        cellRegion[seedi] = nLocal;
        stack.append(seedi);

        while (stack.size())
        {
            const label celli = stack.remove();
            const cell& cFaces = cells[celli];

            forAll(cFaces, i)
            {
                const label facei = cFaces[i];

                if (blockedFace.size() && blockedFace[facei])
                {
                    continue;
                }

                label other = -1;
                if (facei < nInt)
                {
                    other = (own[facei] == celli ? nei[facei] : own[facei]);
                }
                else if (partner.size())
                {
                    Map<label>::const_iterator fnd = partner.find(facei);
                    if (fnd != partner.end())
                    {
                        const label pf = fnd();
                        other =
                        (
                            pf < nInt
                          ? (own[pf] == celli ? nei[pf] : own[pf])
                          : own[pf]
                        );
                    }
                }

                if (other != -1 && cellRegion[other] == -1)
                {
                    cellRegion[other] = nLocal;
                    stack.append(other);
                }
            }
        }

        nLocal++;
    }


    // 2. Merge regions across coupled faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const globalIndex globalRegions(nLocal);

    labelList regionLabel(nLocal);
    forAll(regionLabel, regioni)
    {
        regionLabel[regioni] = globalRegions.toGlobal(regioni);
    }

    propagateMin(blockedFace, regionLabel);


    // 3. Compact: the processor holding a region's minimum label numbers it
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList localToCompact(nLocal, -1);
    label nCompact = 0;
    forAll(regionLabel, regioni)
    {
        const label l = regionLabel[regioni];
        if (globalRegions.isLocal(l))
        {
            const label li = globalRegions.toLocal(l);
            if (localToCompact[li] == -1)
            {
                localToCompact[li] = -2;    // used, numbered below
            }
        }
    }
    forAll(localToCompact, li)
    {
        if (localToCompact[li] == -2)
        {
            localToCompact[li] = nCompact++;
        }
    }

    globalNumberingPtr_.reset(new globalIndex(nCompact));
    const globalIndex& compact = globalNumberingPtr_();

    labelList compactLabel(nLocal, labelMax);
    forAll(regionLabel, regioni)
    {
        const label l = regionLabel[regioni];
        if (globalRegions.isLocal(l))
        {
            compactLabel[regioni] =
                compact.toGlobal(localToCompact[globalRegions.toLocal(l)]);
        }
    }

    propagateMin(blockedFace, compactLabel);

    forAll(cellRegion, celli)
    {
        cellRegion[celli] = compactLabel[cellRegion[celli]];
    }
}


// ************************************************************************* //
