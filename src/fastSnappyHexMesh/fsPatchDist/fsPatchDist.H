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
    Foam::fsPatchDist

Description
    snappySnapDriver::edgePatchDist (distance of every edge to the nearest
    patch point, as found by PointEdgeWave<pointEdgePoint>) on the threaded,
    exact wave of fsPointEdgeWave.H.

    Returns false (nothing done) when fsWave::usable() is false.

SourceFiles
    fsPatchDist.C

\*---------------------------------------------------------------------------*/

#ifndef fsPatchDist_H
#define fsPatchDist_H

#include "polyMesh.H"
#include "indirectPrimitivePatch.H"

namespace Foam
{
namespace fsPatchDist
{

//- Edge distances of snappySnapDriver::edgePatchDist; false: not done
bool edgePatchDist
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    scalarField& edgeDist
);

} // End namespace fsPatchDist
} // End namespace Foam

#endif

// ************************************************************************* //
