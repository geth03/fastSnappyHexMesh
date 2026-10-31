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

Class
    Foam::fastRegionSplit

Description
    Drop-in replacement for regionSplit (global regions) as used by
    snappyHexMesh: cells connected through faces that are not blocked get
    the same region. Regions are numbered compactly per originating
    processor, exactly like regionSplit (globalNumbering(), nRegions()).

    regionSplit in OpenFOAM-13 seeds every face with its own number and
    runs a FaceCellWave<minData> until the minimum has spread through each
    region; the work grows with nCells x region diameter. Here:
      1. local flood fill (stack based, O(nCells))
      2. min-label propagation on the much smaller graph of local regions
         across coupled faces (processor + cyclic), only if coupled
         faces exist
      3. compaction: the owner of a region's minimum label numbers it,
         the number is spread with the same propagation.

SourceFiles
    fastRegionSplit.C

\*---------------------------------------------------------------------------*/

#ifndef fastRegionSplit_H
#define fastRegionSplit_H

#include "polyMesh.H"
#include "globalIndex.H"
#include "labelPair.H"
#include "autoPtr.H"

namespace Foam
{

class fastRegionSplit
:
    public labelList
{
    // Private Data

        const polyMesh& mesh_;

        autoPtr<globalIndex> globalNumberingPtr_;


    // Private Member Functions

        //- Spread the minimum of regionLabel over local regions that are
        //  connected through unblocked coupled faces
        void propagateMin
        (
            const boolList& blockedFace,
            labelList& regionLabel
        ) const;


public:

    // Constructors

        fastRegionSplit
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const List<labelPair>& explicitConnections = List<labelPair>()
        );

        fastRegionSplit(const fastRegionSplit&) = delete;


    // Member Functions

        const polyMesh& mesh() const
        {
            return mesh_;
        }

        const globalIndex& globalNumbering() const
        {
            return globalNumberingPtr_();
        }

        label nLocalRegions() const
        {
            return globalNumbering().localSize(Pstream::myProcNo());
        }

        label nRegions() const
        {
            return globalNumbering().size();
        }


    // Member Operators

        void operator=(const fastRegionSplit&) = delete;
};


} // End namespace Foam

#endif

// ************************************************************************* //
