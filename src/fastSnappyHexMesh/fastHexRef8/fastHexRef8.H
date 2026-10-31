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
    Foam::fastHexRef8

Description
    Parallel-consistent copy of hexRef8::consistentSlowRefinement2
    (OpenFOAM-13), used by snappy's shell refinement (nCellsBetweenLevels > 2).

    The original treats every non-internal face as a wall:
      - a face in facesToCheck that is not internal gets a fake neighbour
        one level finer (ownLevel+1), also for processor/cyclic faces that
        have a real neighbour on the other side -> over-refinement next to
        processor boundaries;
      - level-jump seed faces are only searched among internal faces, so a
        jump lying on a processor face is not seeded.
    Since snappy rebalances between iterations the processor boundaries move,
    new faces get the wrong treatment every iteration and shell refinement
    "creeps" (~100 cells per iteration up to the iteration limit, more cells
    than serial).

    Here coupled faces are handled exactly like internal faces, using the
    neighbour level and cell centre from across the coupled patch. In serial
    without cyclics the result is identical to the original.

SourceFiles
    fastHexRef8.C

\*---------------------------------------------------------------------------*/

#ifndef fastHexRef8_H
#define fastHexRef8_H

#include "fsHexRef8.H"

namespace Foam
{
namespace fastHexRef8
{

//- See hexRef8::consistentSlowRefinement2
labelList consistentSlowRefinement2
(
    const fsHexRef8& cutter,
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck
);

//- Copy of the private hexRef8::faceConsistentRefinement
label faceConsistentRefinement
(
    const fsHexRef8& cutter,
    const bool maxSet,
    PackedBoolList& refineCell
);

//- false if env FS_CREEPFIX=0 (use the original OpenFOAM function)
bool creepFixEnabled();

} // End namespace fastHexRef8
} // End namespace Foam

#endif

// ************************************************************************* //
