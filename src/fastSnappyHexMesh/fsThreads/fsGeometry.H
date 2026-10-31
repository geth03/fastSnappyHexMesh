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
    Foam::fsGeometry

Description
    Threaded calculation of the primitiveMesh geometry (face centres, face
    areas, |face areas|, cell centres, cell volumes), bit-identical to
    OpenFOAM-13's primitiveMesh::makeFaceCentresAndAreas and
    makeCellCentresAndVols:
      - faces: face::areaAndCentre per face (independent)
      - cells: per cell, contributions summed in OpenFOAM's order (owner
        faces ascending, then neighbour faces ascending - the order of the
        face loops in makeCellCentresAndVols).
    The results are stored in the mesh's own (private) geometry pointers
    so every later faceCentres()/cellCentres()/... call uses them.

    ensure() does nothing if the geometry is already there, or when only one
    thread is available (fsThreads::nThreads() == 1) - then OpenFOAM computes
    it on demand as before.

SourceFiles
    fsGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef fsGeometry_H
#define fsGeometry_H

#include "primitiveMesh.H"

namespace Foam
{
namespace fsGeometry
{

//- Compute missing face and cell geometry of mesh with threads
void ensure(const primitiveMesh& mesh);

//- Faces per cell in OpenFOAM's accumulation order of the cell geometry
//  (owner faces ascending, then neighbour faces ascending), CSR
void cellFaceOrder
(
    const primitiveMesh& mesh,
    labelList& offsets,
    labelList& cellFaces
);

//- Move the mesh's face and cell geometry out of the mesh (false and
//  nothing moved if it is not complete). For an incremental update.
bool take
(
    const primitiveMesh& mesh,
    autoPtr<vectorField>& fCtrs,
    autoPtr<vectorField>& fAreas,
    autoPtr<scalarField>& magfAreas,
    autoPtr<vectorField>& cCtrs,
    autoPtr<scalarField>& cVols
);

//- Hand face and cell geometry to a mesh without geometry
void give
(
    const primitiveMesh& mesh,
    autoPtr<vectorField>& fCtrs,
    autoPtr<vectorField>& fAreas,
    autoPtr<scalarField>& magfAreas,
    autoPtr<vectorField>& cCtrs,
    autoPtr<scalarField>& cVols
);

//- Recompute the geometry of the given faces in place (same formula as
//  ensure, so the values are the ones a full calculation gives)
void updateFaces
(
    const primitiveMesh& mesh,
    const labelUList& faces,
    vectorField& fCtrs,
    vectorField& fAreas,
    scalarField& magfAreas
);

//- Recompute the geometry of the given cells in place (offsets/cellFaces
//  from cellFaceOrder)
void updateCells
(
    const primitiveMesh& mesh,
    const labelList& offsets,
    const labelList& cellFaces,
    const labelUList& cells,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cCtrs,
    scalarField& cVols
);

//- Build edges + pointEdges (identical to primitiveMesh::calcEdges for
//  meshes without point ordering) with flat arrays instead of per-point
//  lists and linear searches. Single-threaded; no-op if already there.
void ensureEdges(const primitiveMesh& mesh);

//- Print the size of the mesh and of every cached derived quantity [MB]
void reportCaches(const primitiveMesh& mesh, const char* where);

//- Delete cached derived addressing that is rebuilt on demand
//  (cellCells, edgeCells, pointCells, edgeFaces, pointFaces, cellEdges,
//  faceEdges, pointPoints, cellPoints). Keeps cells, edges, pointEdges and
//  geometry.
void clearUnusedAddressing(const primitiveMesh& mesh);

} // End namespace fsGeometry
} // End namespace Foam

#endif

// ************************************************************************* //
