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

  fastSnappyHexMesh - see fsMeshCheck.H

  The check bodies are copies of OpenFOAM-13 meshCheck (polyMeshCheckQuality.C,
  checkMesh.C) without the report output, which the fallback handles.
\*---------------------------------------------------------------------------*/

#include "fsMeshCheck.H"
#include "fsThreads.H"
#include "fsProfiler.H"
#include "meshCheck.H"
#include "syncTools.H"
#include "triPointRef.H"
#include "tetPointRef.H"
#include "pyramidPointFaceRef.H"
#include "polyMeshTetDecomposition.H"
#include "unitConversion.H"
#include "IOmanip.H"
#include "fsGeometry.H"

#include <thread>
#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

// Faces each thread would insert, in its own loop order
typedef List<DynamicList<label>> insertLists;

// Run f(i0, i1, inserts) on chunks of [0, n), one insert list per chunk
template<class F>
void runChunks(const label n, insertLists& parts, const F& f)
{
    const label nt = max(1, min(fsThreads::nThreads(), n/5000));
    parts.setSize(nt);

    if (nt == 1)
    {
        f(0, n, parts[0]);
        return;
    }

    const label chunk = (n + nt - 1)/nt;
    std::vector<std::thread> workers;
    for (label t = 0; t < nt; t++)
    {
        const label i0 = min(n, t*chunk);
        const label i1 = min(n, i0 + chunk);
        DynamicList<label>& part = parts[t];
        workers.emplace_back([&f, i0, i1, &part]() { f(i0, i1, part); });
    }
    for (std::thread& w : workers)
    {
        w.join();
    }
}


// Insert the collected faces in chunk order: the insertion sequence of the
// serial loop
void mergeInto(const insertLists& parts, labelHashSet& set)
{
    forAll(parts, t)
    {
        const DynamicList<label>& part = parts[t];
        forAll(part, i)
        {
            set.insert(part[i]);
        }
    }
}


void appendOnce(DynamicList<label>& part, const label facei)
{
    if (part.empty() || part.last() != facei)
    {
        part.append(facei);
    }
}


void checkFaceOrthogonality
(
    const scalar orthWarn,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& set
)
{
    FS_PROF("checkFaceOrthogonality");
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const scalar severeNonorthogonalityThreshold = ::cos(orthWarn);

    pointField neiCc(mesh.nFaces() - mesh.nInternalFaces());
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    auto nonOrth = [&](const vector& s, const vector& d)
    {
        const scalar dDotS = (d & s)/(mag(d)*mag(s) + vSmall);
        return dDotS < severeNonorthogonalityThreshold;
    };

    insertLists parts;
    runChunks
    (
        checkFaces.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const label facei = checkFaces[i];
                const point& ownCc = cellCentres[own[facei]];

                if (mesh.isInternalFace(facei))
                {
                    if
                    (
                        nonOrth
                        (
                            faceAreas[facei],
                            cellCentres[nei[facei]] - ownCc
                        )
                    )
                    {
                        part.append(facei);
                    }
                }
                else if (patches[patches.whichPatch(facei)].coupled())
                {
                    if
                    (
                        nonOrth
                        (
                            faceAreas[facei],
                            neiCc[facei-mesh.nInternalFaces()] - ownCc
                        )
                    )
                    {
                        part.append(facei);
                    }
                }
            }
        }
    );
    mergeInto(parts, set);

    forAll(baffles, i)
    {
        const label face0 = baffles[i].first();
        const label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];

        if (nonOrth(faceAreas[face0], cellCentres[own[face1]] - ownCc))
        {
            set.insert(face0);
        }
    }
}


void checkFaceAngles
(
    const scalar maxConcave,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet& set
)
{
    FS_PROF("checkFaceAngles");
    if (maxConcave < -small || maxConcave > degToRad(180)+small)
    {
        FatalErrorInFunction
            << "maxConcave should be [0..180] degrees but is "
            << radToDeg(maxConcave) << abort(FatalError);
    }

    const scalar maxSin = Foam::sin(maxConcave);

    const faceList& fcs = mesh.faces();

    insertLists parts;
    runChunks
    (
        checkFaces.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const label facei = checkFaces[i];

                const face& f = fcs[facei];

                vector faceNormal = faceAreas[facei];
                faceNormal /= mag(faceNormal) + vSmall;

                // Get edge from f[0] to f[size-1];
                vector ePrev(p[f.first()] - p[f.last()]);
                scalar magEPrev = mag(ePrev);
                ePrev /= magEPrev + vSmall;

                forAll(f, fp0)
                {
                    // Get vertex after fp
                    const label fp1 = f.fcIndex(fp0);

                    // Normalised vector between two consecutive points
                    vector e10(p[f[fp1]] - p[f[fp0]]);
                    const scalar magE10 = mag(e10);
                    e10 /= magE10 + vSmall;

                    if (magEPrev > small && magE10 > small)
                    {
                        vector edgeNormal = ePrev ^ e10;
                        const scalar magEdgeNormal = mag(edgeNormal);

                        if (magEdgeNormal >= maxSin)
                        {
                            edgeNormal /= magEdgeNormal;

                            if ((edgeNormal & faceNormal) < small)
                            {
                                appendOnce(part, facei);
                            }
                        }
                    }

                    ePrev = e10;
                    magEPrev = magE10;
                }
            }
        }
    );
    mergeInto(parts, set);
}


void checkFaceSkewness
(
    const scalar internalSkew,
    const scalar boundarySkew,
    const polyMesh& mesh,
    const pointField& points,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& set
)
{
    FS_PROF("checkFaceSkewness");
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    pointField neiCc;
    syncTools::swapBoundaryCellPositions(mesh, cellCentres, neiCc);

    insertLists parts;
    runChunks
    (
        checkFaces.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const label facei = checkFaces[i];

                if (mesh.isInternalFace(facei))
                {
                    const scalar skewness = meshCheck::faceSkewness
                    (
                        mesh,
                        points,
                        faceCentres,
                        faceAreas,

                        facei,
                        cellCentres[own[facei]],
                        cellCentres[nei[facei]]
                    );

                    if (skewness > internalSkew)
                    {
                        part.append(facei);
                    }
                }
                else if (patches[patches.whichPatch(facei)].coupled())
                {
                    const scalar skewness = meshCheck::faceSkewness
                    (
                        mesh,
                        points,
                        faceCentres,
                        faceAreas,

                        facei,
                        cellCentres[own[facei]],
                        neiCc[facei-mesh.nInternalFaces()]
                    );

                    if (skewness > internalSkew)
                    {
                        part.append(facei);
                    }
                }
                else
                {
                    const scalar skewness = meshCheck::boundaryFaceSkewness
                    (
                        mesh,
                        points,
                        faceCentres,
                        faceAreas,

                        facei,
                        cellCentres[own[facei]]
                    );

                    if (skewness > boundarySkew)
                    {
                        part.append(facei);
                    }
                }
            }
        }
    );
    mergeInto(parts, set);

    forAll(baffles, i)
    {
        const label face0 = baffles[i].first();
        const label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];
        const point& neiCc = cellCentres[own[face1]];

        const scalar skewness = meshCheck::faceSkewness
        (
            mesh,
            points,
            faceCentres,
            faceAreas,

            face0,
            ownCc,
            neiCc
        );

        if (skewness > internalSkew)
        {
            set.insert(face0);
        }
    }
}


void checkFaceWeights
(
    const scalar warnWeight,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& set
)
{
    FS_PROF("checkFaceWeights");
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    pointField neiCc(mesh.nFaces()-mesh.nInternalFaces());
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    insertLists parts;
    runChunks
    (
        checkFaces.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const label facei = checkFaces[i];

                const point& fc = faceCentres[facei];
                const vector& fa = faceAreas[facei];

                const scalar dOwn = mag(fa & (fc-cellCentres[own[facei]]));

                if (mesh.isInternalFace(facei))
                {
                    const scalar dNei = mag(fa & (cellCentres[nei[facei]]-fc));
                    const scalar weight =
                        min(dNei, dOwn)/(dNei + dOwn + vSmall);

                    if (weight < warnWeight)
                    {
                        part.append(facei);
                    }
                }
                else if (patches[patches.whichPatch(facei)].coupled())
                {
                    const scalar dNei =
                        mag(fa & (neiCc[facei-mesh.nInternalFaces()]-fc));
                    const scalar weight =
                        min(dNei, dOwn)/(dNei + dOwn + vSmall);

                    if (weight < warnWeight)
                    {
                        part.append(facei);
                    }
                }
            }
        }
    );
    mergeInto(parts, set);

    forAll(baffles, i)
    {
        const label face0 = baffles[i].first();
        const label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];
        const point& fc = faceCentres[face0];
        const vector& fa = faceAreas[face0];

        const scalar dOwn = mag(fa & (fc-ownCc));
        const scalar dNei = mag(fa & (cellCentres[own[face1]]-fc));
        const scalar weight = min(dNei, dOwn)/(dNei + dOwn + vSmall);

        if (weight < warnWeight)
        {
            set.insert(face0);
        }
    }
}


void checkVolRatio
(
    const scalar warnRatio,
    const polyMesh& mesh,
    const scalarField& cellVolumes,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& set
)
{
    FS_PROF("checkVolRatio");
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    scalarField neiVols(mesh.nFaces()-mesh.nInternalFaces());
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiVols[facei-mesh.nInternalFaces()] = cellVolumes[own[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiVols);

    insertLists parts;
    runChunks
    (
        checkFaces.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const label facei = checkFaces[i];

                const scalar ownVol = mag(cellVolumes[own[facei]]);

                scalar neiVol = -great;

                if (mesh.isInternalFace(facei))
                {
                    neiVol = mag(cellVolumes[nei[facei]]);
                }
                else if (patches[patches.whichPatch(facei)].coupled())
                {
                    neiVol = mag(neiVols[facei-mesh.nInternalFaces()]);
                }

                if (neiVol >= 0)
                {
                    const scalar ratio =
                        min(ownVol, neiVol)/(max(ownVol, neiVol) + vSmall);

                    if (ratio < warnRatio)
                    {
                        part.append(facei);
                    }
                }
            }
        }
    );
    mergeInto(parts, set);

    forAll(baffles, i)
    {
        const label face0 = baffles[i].first();
        const label face1 = baffles[i].second();

        const scalar ownVol = mag(cellVolumes[own[face0]]);
        const scalar neiVol = mag(cellVolumes[own[face1]]);

        if (neiVol >= 0)
        {
            const scalar ratio =
                min(ownVol, neiVol)/(max(ownVol, neiVol) + vSmall);

            if (ratio < warnRatio)
            {
                set.insert(face0);
            }
        }
    }
}


void checkFaceTwist
(
    const scalar minTwist,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet& set
)
{
    FS_PROF("checkFaceTwist");
    if (minTwist < -1-small || minTwist > 1+small)
    {
        FatalErrorInFunction
            << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    pointField neiCc(mesh.nFaces()-mesh.nInternalFaces());
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    insertLists parts;
    runChunks
    (
        checkFaces.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const label facei = checkFaces[i];

                const face& f = fcs[facei];

                if (f.size() > 3)
                {
                    vector nf(Zero);

                    if (mesh.isInternalFace(facei))
                    {
                        nf = cellCentres[nei[facei]] - cellCentres[own[facei]];
                        nf /= mag(nf) + vSmall;
                    }
                    else if (patches[patches.whichPatch(facei)].coupled())
                    {
                        nf =
                            neiCc[facei-mesh.nInternalFaces()]
                          - cellCentres[own[facei]];
                        nf /= mag(nf) + vSmall;
                    }
                    else
                    {
                        nf = faceCentres[facei] - cellCentres[own[facei]];
                        nf /= mag(nf) + vSmall;
                    }

                    if (nf != vector::zero)
                    {
                        const point& fc = faceCentres[facei];

                        forAll(f, fpI)
                        {
                            vector triArea
                            (
                                triPointRef
                                (
                                    p[f[fpI]],
                                    p[f.nextLabel(fpI)],
                                    fc
                                ).area()
                            );

                            const scalar magTri = mag(triArea);

                            if
                            (
                                magTri > vSmall
                             && ((nf & triArea/magTri) < minTwist)
                            )
                            {
                                part.append(facei);
                                break;
                            }
                        }
                    }
                }
            }
        }
    );
    mergeInto(parts, set);
}


// Determinant check of one cell, as meshCheck::checkCellDeterminant
inline bool cellDetBad
(
    const cell& cFaces,
    const vectorField& faceAreas,
    const scalar warnDet
)
{
    tensor areaSum(Zero);
    scalar magAreaSum = 0;

    forAll(cFaces, cFacei)
    {
        const label facei = cFaces[cFacei];

        const scalar magArea = mag(faceAreas[facei]);

        magAreaSum += magArea;
        areaSum += faceAreas[facei]*(faceAreas[facei]/(magArea + vSmall));
    }

    // Scalad by 1/(cube root(1/3)) = 27
    const scalar scaledDet = 27*det(areaSum/(magAreaSum + vSmall));

    return scaledDet < warnDet;
}


void checkCellDeterminant
(
    const scalar warnDet,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    labelHashSet& set
)
{
    FS_PROF("checkCellDeterminant");
    const cellList& cells = mesh.cells();

    // As meshCheck::getAffectedCells: the hash table order of the cells is
    // the order in which their faces are inserted.
    // A label hashes to itself and the table iterates its buckets in
    // ascending order, so when every cell label is below the final table
    // size the order is ascending: then a mark list gives the same list
    // without the hash table.
    labelList affectedCells;
    {
        FS_PROF("affectedCells");

        const labelList& own = mesh.faceOwner();
        const labelList& nei = mesh.faceNeighbour();

        List<char> isAffected(mesh.nCells(), char(0));
        label nAffected = 0;
        label maxCelli = -1;
        auto mark = [&](const label celli)
        {
            if (!isAffected[celli])
            {
                isAffected[celli] = 1;
                nAffected++;
                maxCelli = max(maxCelli, celli);
            }
        };
        forAll(checkFaces, i)
        {
            const label facei = checkFaces[i];

            mark(own[facei]);

            if (mesh.isInternalFace(facei))
            {
                mark(nei[facei]);
            }
        }

        // Final size of labelHashSet(2*checkFaces.size()) after nAffected
        // inserts: doubled while the load is above 0.8
        label tableSize = HashTableCore::canonicalSize(2*checkFaces.size());
        if (tableSize == 0 && nAffected > 0)
        {
            tableSize = 2;
        }
        while (tableSize > 0 && double(nAffected)/tableSize > 0.8)
        {
            tableSize *= 2;
        }

        if (maxCelli < tableSize)
        {
            affectedCells.setSize(nAffected);
            label n = 0;
            for (label celli = 0; celli <= maxCelli; celli++)
            {
                if (isAffected[celli])
                {
                    affectedCells[n++] = celli;
                }
            }
        }
        else
        {
            labelHashSet cellSet(2*checkFaces.size());

            forAll(checkFaces, i)
            {
                const label facei = checkFaces[i];

                cellSet.insert(own[facei]);

                if (mesh.isInternalFace(facei))
                {
                    cellSet.insert(nei[facei]);
                }
            }

            affectedCells = cellSet.toc();
        }
    }

    insertLists parts;
    runChunks
    (
        affectedCells.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const cell& cFaces = cells[affectedCells[i]];

                if (cellDetBad(cFaces, faceAreas, warnDet))
                {
                    forAll(cFaces, cFacei)
                    {
                        part.append(cFaces[cFacei]);
                    }
                }
            }
        }
    );
    mergeInto(parts, set);
}


void checkFacePyramids
(
    const scalar minPyrVol,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& set
)
{
    FS_PROF("checkFacePyramids");

    // check whether face area vector points to the cell with higher label
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    const faceList& f = mesh.faces();

    insertLists parts;
    runChunks
    (
        checkFaces.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const label facei = checkFaces[i];

                // pyramidPointFaceRef(face, apex).mag(p) is
                // (1/3)*(face.area(p) & (apex - face.centre(p))): the face
                // area and centre are the same for both pyramids
                const vector fArea = f[facei].area(p);
                const point fCentre = f[facei].centre(p);

                // Create the owner pyramid - it will have negative volume
                const scalar pyrVol =
                    (1.0/3.0)*(fArea & (cellCentres[own[facei]] - fCentre));

                if (pyrVol > -minPyrVol)
                {
                    appendOnce(part, facei);
                }

                if (mesh.isInternalFace(facei))
                {
                    // Create the neighbour pyramid - it will have positive
                    // volume
                    const scalar pyrVol =
                        (1.0/3.0)
                       *(fArea & (cellCentres[nei[facei]] - fCentre));

                    if (pyrVol < minPyrVol)
                    {
                        appendOnce(part, facei);
                    }
                }
            }
        }
    );
    mergeInto(parts, set);

    forAll(baffles, i)
    {
        const label face0 = baffles[i].first();
        const label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];

        // Create the owner pyramid - it will have negative volume
        const scalar pyrVolOwn = pyramidPointFaceRef
        (
            f[face0],
            ownCc
        ).mag(p);

        if (pyrVolOwn > -minPyrVol)
        {
            set.insert(face0);
        }

        // Create the neighbour pyramid - it will have positive volume
        const scalar pyrVolNbr =
            pyramidPointFaceRef(f[face0], cellCentres[own[face1]]).mag(p);

        if (pyrVolNbr < minPyrVol)
        {
            set.insert(face0);
        }
    }
}


// meshCheck::checkFaceTet (namespace scope in OpenFOAM) without report
bool checkFaceTet
(
    const primitiveMesh& mesh,
    const scalar minTetQuality,
    const pointField& p,
    const label facei,
    const point& fc,    // face centre
    const point& cc     // cell centre
)
{
    const face& f = mesh.faces()[facei];

    forAll(f, fp)
    {
        const scalar tetQual = tetPointRef
        (
            p[f[fp]],
            p[f.nextLabel(fp)],
            fc,
            cc
        ).quality();

        if (tetQual < minTetQuality)
        {
            return true;
        }
    }

    return false;
}


void checkFaceTets
(
    const scalar minTetQuality,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& set
)
{
    FS_PROF("checkFaceTets");

    // check whether decomposing each cell into tets results in
    // positive volume, non-flat tets
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces() - mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei - mesh.nInternalFaces()] = cellCentres[own[facei]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    insertLists parts;
    runChunks
    (
        checkFaces.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const label facei = checkFaces[i];

                // Create the owner pyramid - note: exchange cell and face
                // centre to get positive volume.
                if
                (
                    checkFaceTet
                    (
                        mesh,
                        minTetQuality,
                        p,
                        facei,
                        cellCentres[own[facei]],    // face centre
                        faceCentres[facei]          // cell centre
                    )
                )
                {
                    appendOnce(part, facei);
                }

                if (mesh.isInternalFace(facei))
                {
                    // Create the neighbour tets - they will have positive
                    // volume
                    if
                    (
                        checkFaceTet
                        (
                            mesh,
                            minTetQuality,
                            p,
                            facei,
                            faceCentres[facei],         // face centre
                            cellCentres[nei[facei]]     // cell centre
                        )
                    )
                    {
                        appendOnce(part, facei);
                    }

                    if
                    (
                        polyMeshTetDecomposition::findSharedBasePoint
                        (
                            mesh,
                            facei,
                            minTetQuality,
                            false
                        ) == -1
                    )
                    {
                        appendOnce(part, facei);
                    }
                }
                else
                {
                    const label patchi = patches.whichPatch(facei);

                    if (patches[patchi].coupled())
                    {
                        if
                        (
                            polyMeshTetDecomposition::findSharedBasePoint
                            (
                                mesh,
                                facei,
                                neiCc[facei - mesh.nInternalFaces()],
                                minTetQuality,
                                false
                            ) == -1
                        )
                        {
                            appendOnce(part, facei);
                        }
                    }
                    else
                    {
                        if
                        (
                            polyMeshTetDecomposition::findBasePoint
                            (
                                mesh,
                                facei,
                                minTetQuality,
                                false
                            ) == -1
                        )
                        {
                            appendOnce(part, facei);
                        }
                    }
                }
            }
        }
    );
    mergeInto(parts, set);

    forAll(baffles, i)
    {
        const label face0 = baffles[i].first();
        const label face1 = baffles[i].second();

        if
        (
            checkFaceTet
            (
                mesh,
                minTetQuality,
                p,
                face0,
                cellCentres[own[face0]],    // face centre
                faceCentres[face0]          // cell centre
            )
        )
        {
            set.insert(face0);
        }

        // Create the neighbour tets - they will have positive volume
        if
        (
            checkFaceTet
            (
                mesh,
                minTetQuality,
                p,
                face0,
                faceCentres[face0],         // face centre
                cellCentres[own[face1]]     // cell centre
            )
        )
        {
            set.insert(face0);
        }

        if
        (
            polyMeshTetDecomposition::findSharedBasePoint
            (
                mesh,
                face0,
                cellCentres[own[face1]],
                minTetQuality,
                false
            ) == -1
        )
        {
            set.insert(face0);
        }
    }
}


void checkFaceFlatness
(
    const scalar minFlatness,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet& set
)
{
    FS_PROF("checkFaceFlatness");

    if (minFlatness < -small || minFlatness > 1+small)
    {
        FatalErrorInFunction
            << "minFlatness should be [0..1] but is now " << minFlatness
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    insertLists parts;
    runChunks
    (
        checkFaces.size(),
        parts,
        [&](const label i0, const label i1, DynamicList<label>& part)
        {
            for (label i = i0; i < i1; i++)
            {
                const label facei = checkFaces[i];

                const face& f = fcs[facei];

                if (f.size() > 3)
                {
                    const point& fc = faceCentres[facei];

                    // Sum triangle areas
                    scalar sumArea = 0.0;

                    forAll(f, fp)
                    {
                        sumArea += triPointRef
                        (
                            p[f[fp]],
                            p[f.nextLabel(fp)],
                            fc
                        ).mag();
                    }

                    if (sumArea/mag(faceAreas[facei]) < minFlatness)
                    {
                        part.append(facei);
                    }
                }
            }
        }
    );
    mergeInto(parts, set);
}


// The two versions call different collectives, so all ranks must agree
bool useOriginal(const bool report)
{
    const char* e = std::getenv("FS_FASTCHECK");
    bool orig = report || (e && std::string(e) == "0");
    reduce(orig, orOp<bool>());
    return orig;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::fsMeshCheck::checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& wrongFaces
)
{
    if (useOriginal(report))
    {
        return meshCheck::checkMesh
        (
            report,
            mesh,
            dict,
            checkFaces,
            baffles,
            wrongFaces
        );
    }

    const scalar maxNonOrtho
    (
        dict.lookup<scalar>("maxNonOrtho", unitDegrees, true)
    );
    const scalar minVol
    (
        dict.lookup<scalar>("minVol", true)
    );
    const scalar minTetQuality
    (
        dict.lookup<scalar>("minTetQuality", true)
    );
    const scalar maxConcave
    (
        dict.lookup<scalar>("maxConcave", unitDegrees, true)
    );
    const scalar maxIntSkew
    (
        dict.lookup<scalar>("maxInternalSkewness", true)
    );
    const scalar maxBounSkew
    (
        dict.lookup<scalar>("maxBoundarySkewness", true)
    );
    const scalar minWeight
    (
        dict.lookup<scalar>("minFaceWeight", true)
    );
    const scalar minVolRatio
    (
        dict.lookup<scalar>("minVolRatio", true)
    );
    const scalar minTwist
    (
        dict.lookup<scalar>("minTwist", true)
    );
    const scalar minFaceFlatness
    (
        dict.found("minFaceFlatness")
      ? dict.lookup<scalar>("minFaceFlatness", true)
      : -1.0
    );
    const scalar minDet
    (
        dict.lookup<scalar>("minDeterminant", true)
    );

    // Lazily built mesh data the threads read: build it here
    mesh.faces();
    mesh.faceOwner();
    mesh.faceNeighbour();
    mesh.cells();

    // Counter for bad faces
    label nWrongFaces = 0;

    Info<< "Checking faces in error :" << endl;

    if (maxNonOrtho < degToRad(180.0)-small)
    {
        checkFaceOrthogonality
        (
            maxNonOrtho,
            mesh,
            mesh.cellCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    non-orthogonality > "
            << setw(3) << radToDeg(maxNonOrtho)
            << " degrees                        : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVol > -great)
    {
        const scalar refVol = pow3(mesh.bounds().minDim());

        checkFacePyramids
        (
            minVol*refVol,
            mesh,
            mesh.cellCentres(),
            mesh.points(),
            checkFaces,
            baffles,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face pyramid volume < "
            << setw(5) << minVol << "                 : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTetQuality > -great)
    {
        checkFaceTets
        (
            minTetQuality,
            mesh,
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            baffles,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face-decomposition tet quality < "
            << setw(5) << minTetQuality << "      : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxConcave < degToRad(180.0)-small)
    {
        checkFaceAngles
        (
            maxConcave,
            mesh,
            mesh.faceAreas(),
            mesh.points(),
            checkFaces,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with concavity > "
            << setw(3) << radToDeg(maxConcave)
            << " degrees                     : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxIntSkew > 0 || maxBounSkew > 0)
    {
        checkFaceSkewness
        (
            maxIntSkew,
            maxBounSkew,
            mesh,
            mesh.points(),
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with skewness > "
            << setw(3) << maxIntSkew
            << " (internal) or " << setw(3) << maxBounSkew
            << " (boundary) : " << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minWeight >= 0 && minWeight < 1)
    {
        checkFaceWeights
        (
            minWeight,
            mesh,
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with interpolation weights (0..1)  < "
            << setw(5) << minWeight
            << "       : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVolRatio >= 0)
    {
        checkVolRatio
        (
            minVolRatio,
            mesh,
            mesh.cellVolumes(),
            checkFaces,
            baffles,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with volume ratio of neighbour cells < "
            << setw(5) << minVolRatio
            << "     : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTwist > -1)
    {
        checkFaceTwist
        (
            minTwist,
            mesh,
            mesh.cellCentres(),
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face twist < "
            << setw(5) << minTwist
            << "                          : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minFaceFlatness > -small)
    {
        checkFaceFlatness
        (
            minFaceFlatness,
            mesh,
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with flatness < "
            << setw(5) << minFaceFlatness
            << "                      : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minDet > -1)
    {
        checkCellDeterminant
        (
            minDet,
            mesh,
            mesh.faceAreas(),
            checkFaces,
            wrongFaces
        );

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces on cells with determinant < "
            << setw(5) << minDet << "                : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    return nWrongFaces > 0;
}


bool Foam::fsMeshCheck::checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    labelHashSet& wrongFaces
)
{
    const List<labelPair> emptyBaffles;

    return checkMesh
    (
        report,
        mesh,
        dict,
        checkFaces,
        emptyBaffles,
        wrongFaces
    );
}


bool Foam::fsMeshCheck::checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    labelHashSet& wrongFaces
)
{
    return checkMesh
    (
        report,
        mesh,
        dict,
        identityMap(mesh.nFaces()),
        wrongFaces
    );
}



// * * * * * * * * * * * * * * incrementalCheck  * * * * * * * * * * * * * * //

namespace
{

// "0": off, "verify": compare with the full calculation
std::string incrMode()
{
    const char* e = std::getenv("FS_INCRCHECK");
    return e ? std::string(e) : std::string();
}


bool incrUsable()
{
    return !Pstream::parRun() && incrMode() != "0";
}


// Settings of the motion check in checkMesh's order
scalarList checkSettings(const dictionary& dict)
{
    scalarList s(11);
    s[0] = dict.lookup<scalar>("maxNonOrtho", unitDegrees, true);
    s[1] = dict.lookup<scalar>("minVol", true);
    s[2] = dict.lookup<scalar>("minTetQuality", true);
    s[3] = dict.lookup<scalar>("maxConcave", unitDegrees, true);
    s[4] = dict.lookup<scalar>("maxInternalSkewness", true);
    s[5] = dict.lookup<scalar>("maxBoundarySkewness", true);
    s[6] = dict.lookup<scalar>("minFaceWeight", true);
    s[7] = dict.lookup<scalar>("minVolRatio", true);
    s[8] = dict.lookup<scalar>("minTwist", true);
    s[9] =
    (
        dict.found("minFaceFlatness")
      ? dict.lookup<scalar>("minFaceFlatness", true)
      : -1.0
    );
    s[10] = dict.lookup<scalar>("minDeterminant", true);
    return s;
}


inline bool samePoint(const point& a, const point& b)
{
    return std::memcmp(&a, &b, sizeof(point)) == 0;
}


bool sameValues(const scalarList& a, const scalarList& b)
{
    return
        a.size() == b.size()
     && (a.empty() || std::memcmp(a.cdata(), b.cdata(), a.byteSize()) == 0);
}


template<class FieldType>
label nBitDiff(const FieldType& a, const FieldType& b)
{
    label n = 0;
    forAll(a, i)
    {
        if (std::memcmp(&a[i], &b[i], sizeof(a[i])) != 0)
        {
            n++;
        }
    }
    return n;
}

} // End anonymous namespace


Foam::fsMeshCheck::incrementalCheck::incrementalCheck()
:
    nPoints_(-1),
    nFaces_(-1),
    nCells_(-1),
    facesData_(nullptr),
    givenFCtrs_(nullptr),
    haveGeometry_(false),
    all_(true),
    refVol_(-1),
    haveResults_(false)
{}


void Foam::fsMeshCheck::incrementalCheck::clearGeometry()
{
    fCtrs_.clear();
    fAreas_.clear();
    magfAreas_.clear();
    cCtrs_.clear();
    cVols_.clear();
    haveGeometry_ = false;
}


void Foam::fsMeshCheck::incrementalCheck::clear()
{
    clearGeometry();
    givenFCtrs_ = nullptr;
    nPoints_ = -1;
    nFaces_ = -1;
    nCells_ = -1;
    facesData_ = nullptr;
    points_.clear();
    offsets_.clear();
    cellFaces_.clear();
    recheckFaces_.clear();
    recheckCells_.clear();
    all_ = true;
    settings_.clear();
    refVol_ = -1;
    baffles_.clear();
    faceBad_.clear();
    cellBad_.clear();
    haveResults_ = false;
}


void Foam::fsMeshCheck::incrementalCheck::beforeMove(const polyMesh& mesh)
{
    if (!incrUsable())
    {
        clear();
        return;
    }

    if
    (
        mesh.nPoints() != nPoints_
     || mesh.nFaces() != nFaces_
     || mesh.nCells() != nCells_
     || mesh.faces().cdata() != facesData_
    )
    {
        clear();
        nPoints_ = mesh.nPoints();
        nFaces_ = mesh.nFaces();
        nCells_ = mesh.nCells();
        facesData_ = mesh.faces().cdata();
    }

    // Keep the geometry if it is the one handed to the mesh at the last
    // check (it belongs to points_)
    clearGeometry();
    if (givenFCtrs_)
    {
        if (fsGeometry::take(mesh, fCtrs_, fAreas_, magfAreas_, cCtrs_, cVols_))
        {
            if (fCtrs_.operator->() == givenFCtrs_)
            {
                haveGeometry_ = true;
            }
            else
            {
                fsGeometry::give
                (
                    mesh,
                    fCtrs_,
                    fAreas_,
                    magfAreas_,
                    cCtrs_,
                    cVols_
                );
            }
        }
        givenFCtrs_ = nullptr;
    }
}


void Foam::fsMeshCheck::incrementalCheck::afterMove(const polyMesh& mesh)
{
    if (!incrUsable())
    {
        fsGeometry::ensure(mesh);
        return;
    }

    const pointField& p = mesh.points();
    const faceList& fs = mesh.faces();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    // Faces with a moved point, cells with such a face, faces of those
    // cells
    labelList dirtyFaces;
    all_ = true;

    if (points_.size() == p.size() && nFaces_ == mesh.nFaces())
    {
        List<char> moved(p.size());
        fsThreads::forRange
        (
            p.size(),
            [&](const label i0, const label i1)
            {
                for (label pointi = i0; pointi < i1; pointi++)
                {
                    moved[pointi] = !samePoint(p[pointi], points_[pointi]);
                }
            }
        );

        List<char> faceDirty(fs.size());
        fsThreads::forRange
        (
            fs.size(),
            [&](const label i0, const label i1)
            {
                for (label facei = i0; facei < i1; facei++)
                {
                    const face& f = fs[facei];
                    char d = 0;
                    forAll(f, fp)
                    {
                        if (moved[f[fp]])
                        {
                            d = 1;
                            break;
                        }
                    }
                    faceDirty[facei] = d;
                }
            }
        );

        List<char> cellDirty(mesh.nCells(), char(0));
        forAll(faceDirty, facei)
        {
            if (faceDirty[facei])
            {
                cellDirty[own[facei]] = 1;
                if (facei < mesh.nInternalFaces())
                {
                    cellDirty[nei[facei]] = 1;
                }
            }
        }

        DynamicList<label> dFaces;
        DynamicList<label> rFaces;
        forAll(faceDirty, facei)
        {
            if (faceDirty[facei])
            {
                dFaces.append(facei);
            }
            if
            (
                faceDirty[facei]
             || cellDirty[own[facei]]
             || (facei < mesh.nInternalFaces() && cellDirty[nei[facei]])
            )
            {
                rFaces.append(facei);
            }
        }
        DynamicList<label> rCells;
        forAll(cellDirty, celli)
        {
            if (cellDirty[celli])
            {
                rCells.append(celli);
            }
        }

        dirtyFaces.transfer(dFaces);
        recheckFaces_.transfer(rFaces);
        recheckCells_.transfer(rCells);
        all_ = false;
    }

    // Geometry: update the faces and cells that changed, or full
    if (!all_ && haveGeometry_)
    {
        fsGeometry::updateFaces
        (
            mesh,
            dirtyFaces,
            fCtrs_(),
            fAreas_(),
            magfAreas_()
        );
        if (offsets_.empty())
        {
            fsGeometry::cellFaceOrder(mesh, offsets_, cellFaces_);
        }
        fsGeometry::updateCells
        (
            mesh,
            offsets_,
            cellFaces_,
            recheckCells_,
            fCtrs_(),
            fAreas_(),
            cCtrs_(),
            cVols_()
        );
        givenFCtrs_ = fCtrs_.operator->();
        fsGeometry::give(mesh, fCtrs_, fAreas_, magfAreas_, cCtrs_, cVols_);

        if (incrMode() == "verify")
        {
            // Compare with a full calculation
            autoPtr<vectorField> a, b, d;
            autoPtr<scalarField> c, e;
            fsGeometry::take(mesh, a, b, c, d, e);
            fsGeometry::ensure(mesh);
            const label nDiff =
                nBitDiff(a(), mesh.faceCentres())
              + nBitDiff(b(), mesh.faceAreas())
              + nBitDiff(c(), mesh.magFaceAreas())
              + nBitDiff(d(), mesh.cellCentres())
              + nBitDiff(e(), mesh.cellVolumes());
            Info<< "fsIncrementalCheck geometry: " << dirtyFaces.size()
                << " faces, " << recheckCells_.size() << " cells updated, "
                << nDiff << " values differ from the full calculation"
                << endl;
        }
    }
    else
    {
        clearGeometry();
        fsGeometry::ensure(mesh);

        // Remember the arrays the mesh now holds
        if (fsGeometry::take(mesh, fCtrs_, fAreas_, magfAreas_, cCtrs_, cVols_))
        {
            givenFCtrs_ = fCtrs_.operator->();
            fsGeometry::give
            (
                mesh,
                fCtrs_,
                fAreas_,
                magfAreas_,
                cCtrs_,
                cVols_
            );
        }
    }
    haveGeometry_ = false;

    points_ = p;
}


bool Foam::fsMeshCheck::incrementalCheck::check
(
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& wrongFaces
)
{
    // Only for checkFaces = all faces, serial, unique baffle faces
    bool usable = incrUsable() && checkFaces.size() == mesh.nFaces();
    if (usable)
    {
        forAll(checkFaces, i)
        {
            if (checkFaces[i] != i)
            {
                usable = false;
                break;
            }
        }
    }
    if (usable && baffles.size())
    {
        labelHashSet face0s(2*baffles.size());
        forAll(baffles, i)
        {
            if (!face0s.insert(baffles[i].first()))
            {
                usable = false;
                break;
            }
        }
    }

    if (!usable)
    {
        haveResults_ = false;
        return fsMeshCheck::checkMesh
        (
            false,
            mesh,
            dict,
            checkFaces,
            baffles,
            wrongFaces
        );
    }

    FS_PROF("incrementalCheck");

    const label capacity0 = wrongFaces.capacity();
    const scalarList settings(checkSettings(dict));
    const scalar refVol = pow3(mesh.bounds().minDim());

    const bool incremental =
        haveResults_
     && !all_
     && sameValues(settings, settings_)
     && baffles == baffles_;

    const label nFaces = mesh.nFaces();

    // Lazily built mesh data the threads read: build it here
    const faceList& fs = mesh.faces();
    mesh.faceOwner();
    mesh.faceNeighbour();
    const cellList& cells = mesh.cells();
    const pointField& p = mesh.points();
    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();
    const vectorField& faceAreas = mesh.faceAreas();
    const scalarField& cellVolumes = mesh.cellVolumes();
    (void)fs;

    const labelList allFaces(identityMap(nFaces));

    if (!incremental)
    {
        faceBad_.setSize(9);
        forAll(faceBad_, k)
        {
            faceBad_[k].setSize(nFaces);
        }
        cellBad_.setSize(mesh.nCells());
    }

    // Evaluate check k on the faces to evaluate (all, or the faces changed
    // since the last check), into faceBad_[k]
    auto evaluate = [&]
    (
        const label k,
        const bool full,
        const auto& checkFn
    )
    {
        const labelList& faces = (full ? allFaces : recheckFaces_);
        PackedBoolList& bad = faceBad_[k];

        labelHashSet set;
        checkFn(faces, List<labelPair>(), set);

        if (full)
        {
            bad.reset();
            bad.setSize(nFaces);
        }
        else
        {
            forAll(faces, i)
            {
                bad.unset(faces[i]);
            }
        }
        forAllConstIter(labelHashSet, set, iter)
        {
            bad.set(iter.key());
        }
    };

    // Insert check k's wrong faces (faces ascending, then the baffles) as
    // the full check does
    auto insert = [&](const label k, const auto& checkFn)
    {
        const PackedBoolList& bad = faceBad_[k];
        for (label facei = 0; facei < nFaces; facei++)
        {
            if (bad.get(facei))
            {
                wrongFaces.insert(facei);
            }
        }

        if (baffles.size())
        {
            labelHashSet set;
            checkFn(labelList(), baffles, set);
            forAll(baffles, i)
            {
                if (set.found(baffles[i].first()))
                {
                    wrongFaces.insert(baffles[i].first());
                }
            }
        }
    };

    const scalar maxNonOrtho = settings[0];
    const scalar minVol = settings[1];
    const scalar minTetQuality = settings[2];
    const scalar maxConcave = settings[3];
    const scalar maxIntSkew = settings[4];
    const scalar maxBounSkew = settings[5];
    const scalar minWeight = settings[6];
    const scalar minVolRatio = settings[7];
    const scalar minTwist = settings[8];
    const scalar minFaceFlatness = settings[9];
    const scalar minDet = settings[10];

    // Counter for bad faces
    label nWrongFaces = 0;

    Info<< "Checking faces in error :" << endl;

    if (maxNonOrtho < degToRad(180.0)-small)
    {
        auto fn = [&]
        (
            const labelList& faces,
            const List<labelPair>& bfs,
            labelHashSet& set
        )
        {
            checkFaceOrthogonality
            (
                maxNonOrtho,
                mesh,
                cellCentres,
                faceAreas,
                faces,
                bfs,
                set
            );
        };
        evaluate(0, !incremental, fn);
        insert(0, fn);

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    non-orthogonality > "
            << setw(3) << radToDeg(maxNonOrtho)
            << " degrees                        : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVol > -great)
    {
        auto fn = [&]
        (
            const labelList& faces,
            const List<labelPair>& bfs,
            labelHashSet& set
        )
        {
            checkFacePyramids(minVol*refVol, mesh, cellCentres, p, faces, bfs, set);
        };
        // The threshold scales with the mesh bounds
        evaluate(1, !incremental || refVol != refVol_, fn);
        insert(1, fn);

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face pyramid volume < "
            << setw(5) << minVol << "                 : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTetQuality > -great)
    {
        auto fn = [&]
        (
            const labelList& faces,
            const List<labelPair>& bfs,
            labelHashSet& set
        )
        {
            checkFaceTets
            (
                minTetQuality,
                mesh,
                cellCentres,
                faceCentres,
                p,
                faces,
                bfs,
                set
            );
        };
        evaluate(2, !incremental, fn);
        insert(2, fn);

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face-decomposition tet quality < "
            << setw(5) << minTetQuality << "      : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxConcave < degToRad(180.0)-small)
    {
        auto fn = [&]
        (
            const labelList& faces,
            const List<labelPair>&,
            labelHashSet& set
        )
        {
            checkFaceAngles(maxConcave, mesh, faceAreas, p, faces, set);
        };
        evaluate(3, !incremental, fn);
        insert(3, fn);

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with concavity > "
            << setw(3) << radToDeg(maxConcave)
            << " degrees                     : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxIntSkew > 0 || maxBounSkew > 0)
    {
        auto fn = [&]
        (
            const labelList& faces,
            const List<labelPair>& bfs,
            labelHashSet& set
        )
        {
            checkFaceSkewness
            (
                maxIntSkew,
                maxBounSkew,
                mesh,
                p,
                cellCentres,
                faceCentres,
                faceAreas,
                faces,
                bfs,
                set
            );
        };
        evaluate(4, !incremental, fn);
        insert(4, fn);

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with skewness > "
            << setw(3) << maxIntSkew
            << " (internal) or " << setw(3) << maxBounSkew
            << " (boundary) : " << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minWeight >= 0 && minWeight < 1)
    {
        auto fn = [&]
        (
            const labelList& faces,
            const List<labelPair>& bfs,
            labelHashSet& set
        )
        {
            checkFaceWeights
            (
                minWeight,
                mesh,
                cellCentres,
                faceCentres,
                faceAreas,
                faces,
                bfs,
                set
            );
        };
        evaluate(5, !incremental, fn);
        insert(5, fn);

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with interpolation weights (0..1)  < "
            << setw(5) << minWeight
            << "       : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVolRatio >= 0)
    {
        auto fn = [&]
        (
            const labelList& faces,
            const List<labelPair>& bfs,
            labelHashSet& set
        )
        {
            checkVolRatio(minVolRatio, mesh, cellVolumes, faces, bfs, set);
        };
        evaluate(6, !incremental, fn);
        insert(6, fn);

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with volume ratio of neighbour cells < "
            << setw(5) << minVolRatio
            << "     : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTwist > -1)
    {
        auto fn = [&]
        (
            const labelList& faces,
            const List<labelPair>&,
            labelHashSet& set
        )
        {
            checkFaceTwist
            (
                minTwist,
                mesh,
                cellCentres,
                faceAreas,
                faceCentres,
                p,
                faces,
                set
            );
        };
        evaluate(7, !incremental, fn);
        insert(7, fn);

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face twist < "
            << setw(5) << minTwist
            << "                          : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minFaceFlatness > -small)
    {
        auto fn = [&]
        (
            const labelList& faces,
            const List<labelPair>&,
            labelHashSet& set
        )
        {
            checkFaceFlatness
            (
                minFaceFlatness,
                mesh,
                faceAreas,
                faceCentres,
                p,
                faces,
                set
            );
        };
        evaluate(8, !incremental, fn);
        insert(8, fn);

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with flatness < "
            << setw(5) << minFaceFlatness
            << "                      : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minDet > -1)
    {
        // checkFaces = all faces: the affected cells are all cells, in
        // ascending order (see checkCellDeterminant)
        const labelList allCells(identityMap(mesh.nCells()));
        const labelList& evalCells =
            (incremental ? recheckCells_ : allCells);

        List<char> bad(evalCells.size());
        fsThreads::forRange
        (
            evalCells.size(),
            [&](const label i0, const label i1)
            {
                for (label i = i0; i < i1; i++)
                {
                    bad[i] = cellDetBad(cells[evalCells[i]], faceAreas, minDet);
                }
            }
        );
        forAll(evalCells, i)
        {
            cellBad_.set(evalCells[i], bad[i]);
        }

        forAll(cells, celli)
        {
            if (cellBad_.get(celli))
            {
                const cell& cFaces = cells[celli];
                forAll(cFaces, cFacei)
                {
                    wrongFaces.insert(cFaces[cFacei]);
                }
            }
        }

        const label nNewWrongFaces =
            returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces on cells with determinant < "
            << setw(5) << minDet << "                : "
            << nNewWrongFaces - nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    settings_ = settings;
    refVol_ = refVol;
    baffles_ = baffles;
    haveResults_ = true;

    if (incrMode() == "verify")
    {
        // Same initial table size as the caller's set, so the same
        // insertions give the same iteration order
        labelHashSet full(capacity0);
        fsMeshCheck::checkMesh(false, mesh, dict, checkFaces, baffles, full);
        const labelList a(wrongFaces.toc());
        const labelList b(full.toc());
        Info<< "fsIncrementalCheck check (" << (incremental ? "incr" : "full")
            << ", " << recheckFaces_.size() << " faces): order "
            << (a == b ? "same" : "DIFFERENT") << ", contents "
            << (wrongFaces.sortedToc() == full.sortedToc() ? "same" : "DIFFERENT")
            << " (" << a.size() << " vs " << b.size() << ")" << endl;
    }

    return nWrongFaces > 0;
}


// ************************************************************************* //
