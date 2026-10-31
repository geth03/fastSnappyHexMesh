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

  fastSnappyHexMesh - see fsGeometry.H
\*---------------------------------------------------------------------------*/

#include "fsGeometry.H"
#include "fsThreads.H"
#include "face.H"
#include "UIndirectList.H"
#include "Tuple2.H"

#include <thread>
#include <vector>
#include <algorithm>

// * * * * * * * * * * Access to primitiveMesh geometry pointers * * * * * * //
//
// Standard C++ technique: access checking does not apply to the template
// arguments of an explicit instantiation, so a pointer-to-member of a private
// member can be captured there and returned through a friend function.
// OpenFOAM itself is not modified.

namespace
{

template<class Tag, typename Tag::type Member>
struct robber
{
    friend typename Tag::type get(Tag)
    {
        return Member;
    }
};

struct faceCentresTag
{
    typedef Foam::vectorField* Foam::primitiveMesh::*type;
    friend type get(faceCentresTag);
};
struct faceAreasTag
{
    typedef Foam::vectorField* Foam::primitiveMesh::*type;
    friend type get(faceAreasTag);
};
struct magFaceAreasTag
{
    typedef Foam::scalarField* Foam::primitiveMesh::*type;
    friend type get(magFaceAreasTag);
};
struct cellCentresTag
{
    typedef Foam::vectorField* Foam::primitiveMesh::*type;
    friend type get(cellCentresTag);
};
struct cellVolumesTag
{
    typedef Foam::scalarField* Foam::primitiveMesh::*type;
    friend type get(cellVolumesTag);
};

} // End anonymous namespace

template struct robber<faceCentresTag, &Foam::primitiveMesh::faceCentresPtr_>;
template struct robber<faceAreasTag, &Foam::primitiveMesh::faceAreasPtr_>;
template struct robber
<
    magFaceAreasTag,
    &Foam::primitiveMesh::magFaceAreasPtr_
>;
template struct robber<cellCentresTag, &Foam::primitiveMesh::cellCentresPtr_>;
template struct robber<cellVolumesTag, &Foam::primitiveMesh::cellVolumesPtr_>;

namespace
{
struct edgesTag
{
    typedef Foam::edgeList* Foam::primitiveMesh::*type;
    friend type get(edgesTag);
};
struct pointEdgesTag
{
    typedef Foam::labelListList* Foam::primitiveMesh::*type;
    friend type get(pointEdgesTag);
};
struct nInternal0EdgesTag
{
    typedef Foam::label Foam::primitiveMesh::*type;
    friend type get(nInternal0EdgesTag);
};
struct nInternalPointsTag
{
    typedef Foam::label Foam::primitiveMesh::*type;
    friend type get(nInternalPointsTag);
};
}

namespace
{
#define FS_LLTAG(NAME, MEMBER)                                                 \
struct NAME                                                                    \
{                                                                              \
    typedef Foam::labelListList* Foam::primitiveMesh::*type;                   \
    friend type get(NAME);                                                     \
};
FS_LLTAG(ccTag, ccPtr_)
FS_LLTAG(ecTag, ecPtr_)
FS_LLTAG(pcTag, pcPtr_)
FS_LLTAG(efTag, efPtr_)
FS_LLTAG(pfTag, pfPtr_)
FS_LLTAG(ceTag, cePtr_)
FS_LLTAG(feTag, fePtr_)
FS_LLTAG(ppTag, ppPtr_)
FS_LLTAG(cpTag, cpPtr_)
#undef FS_LLTAG
}

namespace
{
struct cfTag
{
    typedef Foam::cellList* Foam::primitiveMesh::*type;
    friend type get(cfTag);
};
}
template struct robber<cfTag, &Foam::primitiveMesh::cfPtr_>;
template struct robber<ccTag, &Foam::primitiveMesh::ccPtr_>;
template struct robber<ecTag, &Foam::primitiveMesh::ecPtr_>;
template struct robber<pcTag, &Foam::primitiveMesh::pcPtr_>;
template struct robber<efTag, &Foam::primitiveMesh::efPtr_>;
template struct robber<pfTag, &Foam::primitiveMesh::pfPtr_>;
template struct robber<ceTag, &Foam::primitiveMesh::cePtr_>;
template struct robber<feTag, &Foam::primitiveMesh::fePtr_>;
template struct robber<ppTag, &Foam::primitiveMesh::ppPtr_>;
template struct robber<cpTag, &Foam::primitiveMesh::cpPtr_>;

template struct robber<edgesTag, &Foam::primitiveMesh::edgesPtr_>;
template struct robber<pointEdgesTag, &Foam::primitiveMesh::pePtr_>;
template struct robber
<
    nInternal0EdgesTag,
    &Foam::primitiveMesh::nInternal0Edges_
>;
template struct robber
<
    nInternalPointsTag,
    &Foam::primitiveMesh::nInternalPoints_
>;


namespace
{

// Run f(i0, i1) on nt threads over [0, n)
template<class F>
void parallelRange(const Foam::label n, const Foam::label nt, const F& f)
{
    const Foam::label chunk = (n + nt - 1)/nt;
    std::vector<std::thread> workers;
    for (Foam::label t = 0; t < nt; t++)
    {
        const Foam::label i0 = t*chunk;
        const Foam::label i1 = Foam::min(n, i0 + chunk);
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

} // End anonymous namespace


// * * * * * * * * * * * * * * * Functions  * * * * * * * * * * * * * * * * //

namespace
{
    // Approximate heap bytes of a labelListList (data + list headers)
    double llMB(const Foam::labelListList* p)
    {
        if (!p)
        {
            return 0;
        }
        double b = double(p->size())*(sizeof(Foam::labelList) + 16);
        forAll(*p, i)
        {
            b += double((*p)[i].size())*sizeof(Foam::label);
        }
        return b/(1024.0*1024.0);
    }
}


void Foam::fsGeometry::reportCaches
(
    const primitiveMesh& cmesh,
    const char* where
)
{
    primitiveMesh& mesh = const_cast<primitiveMesh&>(cmesh);

    const double mb = 1.0/(1024.0*1024.0);
    const double vMB = sizeof(vector)*mb;
    const double sMB = sizeof(scalar)*mb;

    double faces = double(mesh.nFaces())*(sizeof(face) + 16);
    forAll(mesh.faces(), facei)
    {
        faces += double(mesh.faces()[facei].size())*sizeof(label);
    }

    Info<< "[fsMem] " << where
        << " nCells " << mesh.nCells()
        << " | faces " << faces*mb
        << " points " << mesh.nPoints()*vMB
        << " own+nei " << (mesh.nFaces() + mesh.nInternalFaces())*sizeof(label)*mb
        << " | geom fC " << (mesh.*get(faceCentresTag()) ? mesh.nFaces()*vMB : 0)
        << " fA " << (mesh.*get(faceAreasTag()) ? mesh.nFaces()*vMB : 0)
        << " |fA| " << (mesh.*get(magFaceAreasTag()) ? mesh.nFaces()*sMB : 0)
        << " cC " << (mesh.*get(cellCentresTag()) ? mesh.nCells()*vMB : 0)
        << " cV " << (mesh.*get(cellVolumesTag()) ? mesh.nCells()*sMB : 0)
        << " | edges " << (mesh.*get(edgesTag()) ? mesh.nEdges()*2*sizeof(label)*mb : 0)
        << " pE " << llMB(mesh.*get(pointEdgesTag()))
        << " cells(cf) "
        << llMB(reinterpret_cast<const labelListList*>(mesh.*get(cfTag())))
        << " cc " << llMB(mesh.*get(ccTag()))
        << " ec " << llMB(mesh.*get(ecTag()))
        << " pc " << llMB(mesh.*get(pcTag()))
        << " ef " << llMB(mesh.*get(efTag()))
        << " pf " << llMB(mesh.*get(pfTag()))
        << " ce " << llMB(mesh.*get(ceTag()))
        << " fe " << llMB(mesh.*get(feTag()))
        << " pp " << llMB(mesh.*get(ppTag()))
        << " cp " << llMB(mesh.*get(cpTag()))
        << " MB" << endl;
}


void Foam::fsGeometry::clearUnusedAddressing(const primitiveMesh& cmesh)
{
    primitiveMesh& mesh = const_cast<primitiveMesh&>(cmesh);

    labelListList** ptrs[] =
    {
        &(mesh.*get(ccTag())),
        &(mesh.*get(ecTag())),
        &(mesh.*get(pcTag())),
        &(mesh.*get(efTag())),
        &(mesh.*get(pfTag())),
        &(mesh.*get(ceTag())),
        &(mesh.*get(feTag())),
        &(mesh.*get(ppTag())),
        &(mesh.*get(cpTag()))
    };

    for (labelListList** pp : ptrs)
    {
        if (*pp)
        {
            delete *pp;
            *pp = nullptr;
        }
    }
}


void Foam::fsGeometry::ensureEdges(const primitiveMesh& cmesh)
{
    primitiveMesh& mesh = const_cast<primitiveMesh&>(cmesh);

    edgeList*& edgesPtr = mesh.*get(edgesTag());
    labelListList*& pePtr = mesh.*get(pointEdgesTag());

    // Only the unordered-points case (nInternalPoints == -1, what
    // polyTopoChange without point ordering produces); otherwise OpenFOAM
    // builds the edges on demand as before.
    if (edgesPtr || pePtr || mesh.*get(nInternalPointsTag()) != -1)
    {
        return;
    }

    const faceList& fcs = mesh.faces();
    const label nPoints = mesh.nPoints();

    // 1. All face edges as (lower, higher) pairs, bucketed by lower point
    labelList nPer(nPoints + 1, 0);
    forAll(fcs, facei)
    {
        const face& f = fcs[facei];
        forAll(f, fp)
        {
            const label a = f[fp];
            const label b = f[f.fcIndex(fp)];
            nPer[min(a, b) + 1]++;
        }
    }
    for (label pointi = 0; pointi < nPoints; pointi++)
    {
        nPer[pointi + 1] += nPer[pointi];
    }
    labelList higher(nPer[nPoints]);
    {
        labelList fill(SubList<label>(nPer, nPoints));
        forAll(fcs, facei)
        {
            const face& f = fcs[facei];
            forAll(f, fp)
            {
                const label a = f[fp];
                const label b = f[f.fcIndex(fp)];
                higher[fill[min(a, b)]++] = max(a, b);
            }
        }
    }

    // 2. Per lower point: sort + unique -> edges in lexicographic order,
    //    which is OpenFOAM's final edge order for nInternalPoints == -1
    //    (per point, edges to higher neighbours in ascending order).
    label nEdges = 0;
    labelList edgeStart(nPoints + 1, 0);
    for (label pointi = 0; pointi < nPoints; pointi++)
    {
        label* b = higher.begin() + nPer[pointi];
        label* e = higher.begin() + nPer[pointi + 1];
        std::sort(b, e);
        e = std::unique(b, e);
        edgeStart[pointi] = nEdges;
        nEdges += label(e - b);
        edgeStart[pointi + 1] = nEdges;
    }

    edgeList* es = new edgeList(nEdges);
    {
        label edgei = 0;
        for (label pointi = 0; pointi < nPoints; pointi++)
        {
            const label n = edgeStart[pointi + 1] - edgeStart[pointi];
            const label* b = higher.begin() + nPer[pointi];
            for (label i = 0; i < n; i++)
            {
                (*es)[edgei++] = edge(pointi, b[i]);
            }
        }
    }

    // 3. pointEdges: ascending edge indices per point (edges are visited in
    //    ascending order, so appending keeps every list sorted). A
    //    degenerate edge (a == a) is listed once, as in OpenFOAM.
    labelList nPE(nPoints, 0);
    forAll(*es, edgei)
    {
        const edge& e = (*es)[edgei];
        nPE[e[0]]++;
        if (e[1] != e[0])
        {
            nPE[e[1]]++;
        }
    }
    labelListList* pe = new labelListList(nPoints);
    forAll(*pe, pointi)
    {
        (*pe)[pointi].setSize(nPE[pointi]);
    }
    nPE = 0;
    forAll(*es, edgei)
    {
        const edge& e = (*es)[edgei];
        (*pe)[e[0]][nPE[e[0]]++] = edgei;
        if (e[1] != e[0])
        {
            (*pe)[e[1]][nPE[e[1]]++] = edgei;
        }
    }

    edgesPtr = es;
    pePtr = pe;
    mesh.*get(nInternal0EdgesTag()) = nEdges;
}


namespace
{

// Face centre, area and |area| of facei, as makeFaceCentresAndAreas
inline void faceGeometry
(
    const Foam::faceList& fs,
    const Foam::pointField& p,
    const Foam::label facei,
    Foam::vectorField& fCtrs,
    Foam::vectorField& fAreas,
    Foam::scalarField& magfAreas
)
{
    using namespace Foam;

    const Tuple2<vector, point> areaAndCentre =
        face::areaAndCentre(UIndirectList<point>(p, fs[facei]));

    fCtrs[facei] = areaAndCentre.second();
    fAreas[facei] = areaAndCentre.first();
    magfAreas[facei] = max(mag(fAreas[facei]), rootVSmall);
}


// Cell centre and volume of celli, as makeCellCentresAndVols; the faces of
// the cell in OpenFOAM's accumulation order (cellFaceOrder)
inline void cellGeometry
(
    const Foam::labelList& own,
    const Foam::labelList& offsets,
    const Foam::labelList& cellFaces,
    const Foam::vectorField& fCtrs,
    const Foam::vectorField& fAreas,
    const Foam::label celli,
    Foam::vectorField& cCtrs,
    Foam::scalarField& cVols
)
{
    using namespace Foam;

    const label s = offsets[celli];
    const label e = offsets[celli + 1];

    // Estimated centre: average of face centres
    vector cEst(Zero);
    for (label k = s; k < e; k++)
    {
        cEst += fCtrs[cellFaces[k]];
    }
    const label nCellFaces = e - s;
    cEst /= nCellFaces;

    vector cc(Zero);
    scalar vol = 0.0;
    for (label k = s; k < e; k++)
    {
        const label facei = cellFaces[k];

        // Owner faces come first, so the owner test is the same as
        // OpenFOAM's separate own/nei loops
        const scalar pyr3Vol =
        (
            own[facei] == celli
          ? fAreas[facei] & (fCtrs[facei] - cEst)
          : fAreas[facei] & (cEst - fCtrs[facei])
        );

        const vector pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst;

        cc += pyr3Vol*pc;
        vol += pyr3Vol;
    }

    if (mag(vol) > vSmall)
    {
        cc /= vol;
    }
    else
    {
        cc = cEst;
    }

    cCtrs[celli] = cc;
    cVols[celli] = vol*(1.0/3.0);    // as cellVols *= 1/3
}

} // End anonymous namespace


void Foam::fsGeometry::cellFaceOrder
(
    const primitiveMesh& mesh,
    labelList& offsets,
    labelList& cellFaces
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const label nCells = mesh.nCells();

    // Faces per cell in OpenFOAM's accumulation order: owner faces
    // ascending, then neighbour faces ascending (counting sort).
    offsets.setSize(nCells + 1);
    offsets = 0;
    forAll(own, facei)
    {
        offsets[own[facei] + 1]++;
    }
    forAll(nei, facei)
    {
        offsets[nei[facei] + 1]++;
    }
    for (label celli = 0; celli < nCells; celli++)
    {
        offsets[celli + 1] += offsets[celli];
    }
    cellFaces.setSize(offsets[nCells]);
    {
        labelList fill(SubList<label>(offsets, nCells));
        forAll(own, facei)
        {
            cellFaces[fill[own[facei]]++] = facei;
        }
        forAll(nei, facei)
        {
            cellFaces[fill[nei[facei]]++] = facei;
        }
    }
}


void Foam::fsGeometry::ensure(const primitiveMesh& cmesh)
{
    const label nt = fsThreads::nThreads();
    if (nt <= 1)
    {
        return;
    }

    primitiveMesh& mesh = const_cast<primitiveMesh&>(cmesh);

    vectorField*& fCtrsPtr = mesh.*get(faceCentresTag());
    vectorField*& fAreasPtr = mesh.*get(faceAreasTag());
    scalarField*& magfAreasPtr = mesh.*get(magFaceAreasTag());
    vectorField*& cCtrsPtr = mesh.*get(cellCentresTag());
    scalarField*& cVolsPtr = mesh.*get(cellVolumesTag());

    // Face geometry (all three or none, as calcFaceCentresAndAreas)
    if (!fCtrsPtr && !fAreasPtr && !magfAreasPtr)
    {
        const faceList& fs = mesh.faces();
        const pointField& p = mesh.points();

        vectorField* fCtrs = new vectorField(mesh.nFaces());
        vectorField* fAreas = new vectorField(mesh.nFaces());
        scalarField* magfAreas = new scalarField(mesh.nFaces());

        parallelRange
        (
            fs.size(),
            nt,
            [&](const label i0, const label i1)
            {
                for (label facei = i0; facei < i1; facei++)
                {
                    faceGeometry(fs, p, facei, *fCtrs, *fAreas, *magfAreas);
                }
            }
        );

        fCtrsPtr = fCtrs;
        fAreasPtr = fAreas;
        magfAreasPtr = magfAreas;
    }

    // Cell geometry (both or none, as calcCellCentresAndVols)
    if (!cCtrsPtr && !cVolsPtr && fCtrsPtr && fAreasPtr)
    {
        const vectorField& fCtrs = *fCtrsPtr;
        const vectorField& fAreas = *fAreasPtr;
        const labelList& own = mesh.faceOwner();
        const label nCells = mesh.nCells();

        labelList offsets;
        labelList cellFaces;
        cellFaceOrder(mesh, offsets, cellFaces);

        vectorField* cCtrs = new vectorField(nCells);
        scalarField* cVols = new scalarField(nCells);

        parallelRange
        (
            nCells,
            nt,
            [&](const label i0, const label i1)
            {
                for (label celli = i0; celli < i1; celli++)
                {
                    cellGeometry
                    (
                        own,
                        offsets,
                        cellFaces,
                        fCtrs,
                        fAreas,
                        celli,
                        *cCtrs,
                        *cVols
                    );
                }
            }
        );

        cCtrsPtr = cCtrs;
        cVolsPtr = cVols;
    }
}


bool Foam::fsGeometry::take
(
    const primitiveMesh& cmesh,
    autoPtr<vectorField>& fCtrs,
    autoPtr<vectorField>& fAreas,
    autoPtr<scalarField>& magfAreas,
    autoPtr<vectorField>& cCtrs,
    autoPtr<scalarField>& cVols
)
{
    primitiveMesh& mesh = const_cast<primitiveMesh&>(cmesh);

    vectorField*& fCtrsPtr = mesh.*get(faceCentresTag());
    vectorField*& fAreasPtr = mesh.*get(faceAreasTag());
    scalarField*& magfAreasPtr = mesh.*get(magFaceAreasTag());
    vectorField*& cCtrsPtr = mesh.*get(cellCentresTag());
    scalarField*& cVolsPtr = mesh.*get(cellVolumesTag());

    if (!fCtrsPtr || !fAreasPtr || !magfAreasPtr || !cCtrsPtr || !cVolsPtr)
    {
        return false;
    }

    fCtrs.reset(fCtrsPtr);
    fAreas.reset(fAreasPtr);
    magfAreas.reset(magfAreasPtr);
    cCtrs.reset(cCtrsPtr);
    cVols.reset(cVolsPtr);

    fCtrsPtr = nullptr;
    fAreasPtr = nullptr;
    magfAreasPtr = nullptr;
    cCtrsPtr = nullptr;
    cVolsPtr = nullptr;

    return true;
}


void Foam::fsGeometry::give
(
    const primitiveMesh& cmesh,
    autoPtr<vectorField>& fCtrs,
    autoPtr<vectorField>& fAreas,
    autoPtr<scalarField>& magfAreas,
    autoPtr<vectorField>& cCtrs,
    autoPtr<scalarField>& cVols
)
{
    primitiveMesh& mesh = const_cast<primitiveMesh&>(cmesh);

    vectorField*& fCtrsPtr = mesh.*get(faceCentresTag());
    vectorField*& fAreasPtr = mesh.*get(faceAreasTag());
    scalarField*& magfAreasPtr = mesh.*get(magFaceAreasTag());
    vectorField*& cCtrsPtr = mesh.*get(cellCentresTag());
    scalarField*& cVolsPtr = mesh.*get(cellVolumesTag());

    if (fCtrsPtr || fAreasPtr || magfAreasPtr || cCtrsPtr || cVolsPtr)
    {
        FatalErrorInFunction
            << "Mesh already has geometry" << abort(FatalError);
    }

    fCtrsPtr = fCtrs.ptr();
    fAreasPtr = fAreas.ptr();
    magfAreasPtr = magfAreas.ptr();
    cCtrsPtr = cCtrs.ptr();
    cVolsPtr = cVols.ptr();
}


void Foam::fsGeometry::updateFaces
(
    const primitiveMesh& mesh,
    const labelUList& faces,
    vectorField& fCtrs,
    vectorField& fAreas,
    scalarField& magfAreas
)
{
    const faceList& fs = mesh.faces();
    const pointField& p = mesh.points();

    fsThreads::forRange
    (
        faces.size(),
        [&](const label i0, const label i1)
        {
            for (label i = i0; i < i1; i++)
            {
                faceGeometry(fs, p, faces[i], fCtrs, fAreas, magfAreas);
            }
        }
    );
}


void Foam::fsGeometry::updateCells
(
    const primitiveMesh& mesh,
    const labelList& offsets,
    const labelList& cellFaces,
    const labelUList& cells,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cCtrs,
    scalarField& cVols
)
{
    const labelList& own = mesh.faceOwner();

    fsThreads::forRange
    (
        cells.size(),
        [&](const label i0, const label i1)
        {
            for (label i = i0; i < i1; i++)
            {
                cellGeometry
                (
                    own,
                    offsets,
                    cellFaces,
                    fCtrs,
                    fAreas,
                    cells[i],
                    cCtrs,
                    cVols
                );
            }
        }
    );
}


// ************************************************************************* //
