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

Namespace
    Foam::fsMeshCheck

Description
    Threaded drop-in for the motion check meshCheck::checkMesh (the one with
    checkFaces and baffles, used by scaleMesh, merging and the final check).

    The per-face loops of the checks run on fsThreads::nThreads() threads.
    Every thread collects the faces it would insert, in its own checkFaces
    order; the lists are inserted chunk by chunk afterwards, so wrongFaces
    gets exactly the insertions, in exactly the order, of the original.
    Boundary swaps, baffle loops and all parallel reductions stay on the
    calling thread.

    All ten checks of the motion check are threaded (pyramids, tets and
    flatness too: layer settings that check validity use minVol and
    minTetQuality).

    Falls back to meshCheck::checkMesh when report is true (keeps the
    original messages) or with FS_FASTCHECK=0.

SourceFiles
    fsMeshCheck.C

\*---------------------------------------------------------------------------*/

#ifndef fsMeshCheck_H
#define fsMeshCheck_H

#include "polyMesh.H"
#include "HashSet.H"
#include "labelPair.H"
#include "PackedBoolList.H"
#include "autoPtr.H"

namespace Foam
{
namespace fsMeshCheck
{

/*---------------------------------------------------------------------------*\
                      Class incrementalCheck Declaration
\*---------------------------------------------------------------------------*/

//- Repeated motion checks of one mesh whose points move (motionSmootherAlgo
//  ::scaleMesh): geometry and check results of the faces and cells that no
//  point movement can have changed are taken over from the last check.
//
//  A face's check results depend only on its points and on the geometry of
//  its owner and neighbour cells; a cell's geometry only on its faces. So
//  after a move only the faces with a moved point, and all faces of cells
//  with such a face, are evaluated again; everything else is bit for bit
//  what a full evaluation gives. The wrong faces are inserted in the order
//  of the full check (check by check, faces ascending, then the baffles).
//
//  Used for serial runs with checkFaces = all faces; otherwise (and after
//  any topology change: clear()) the full check runs.
//  FS_INCRCHECK=0: off. FS_INCRCHECK=verify: also run the full geometry
//  and check and report differences.
class incrementalCheck
{
    // Topology of the cached state
    label nPoints_;
    label nFaces_;
    label nCells_;
    const void* facesData_;

    //- Points of the last check (empty: none)
    pointField points_;

    //- Geometry of the last check while the mesh moves; the arrays given
    //  to the mesh (to recognise them)
    autoPtr<vectorField> fCtrs_;
    autoPtr<vectorField> fAreas_;
    autoPtr<scalarField> magfAreas_;
    autoPtr<vectorField> cCtrs_;
    autoPtr<scalarField> cVols_;
    const vectorField* givenFCtrs_;
    bool haveGeometry_;

    //- Faces per cell in the geometry's accumulation order
    labelList offsets_;
    labelList cellFaces_;

    //- Faces and cells to evaluate again after the last move; -1: all
    labelList recheckFaces_;
    labelList recheckCells_;
    bool all_;

    //- Settings and results of the last check
    scalarList settings_;
    scalar refVol_;
    List<labelPair> baffles_;
    List<PackedBoolList> faceBad_;
    PackedBoolList cellBad_;
    bool haveResults_;

    void clearGeometry();


public:

    incrementalCheck();

    //- Forget everything (topology change)
    void clear();

    //- Before the points of mesh are set: keep its geometry
    void beforeMove(const polyMesh& mesh);

    //- After the points are set: geometry of the mesh (updated or full)
    void afterMove(const polyMesh& mesh);

    //- As checkMesh(false, mesh, dict, checkFaces, baffles, wrongFaces)
    bool check
    (
        const polyMesh& mesh,
        const dictionary& dict,
        const labelList& checkFaces,
        const List<labelPair>& baffles,
        labelHashSet& wrongFaces
    );
};


//- See meshCheck::checkMesh
bool checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& wrongFaces
);

//- See meshCheck::checkMesh
bool checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    labelHashSet& wrongFaces
);

//- See meshCheck::checkMesh
bool checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    labelHashSet& wrongFaces
);

} // End namespace fsMeshCheck
} // End namespace Foam

#endif

// ************************************************************************* //
