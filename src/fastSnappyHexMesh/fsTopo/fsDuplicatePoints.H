// fastSnappyHexMesh: copy of OpenFOAM-13 duplicatePoints.H, class renamed to fsDuplicatePoints and
// working on fsPolyTopoChange. Changes are marked "fastSnappy".
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2026 Seyit Temiz
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Modifications
    Modified in 2026 by Seyit Temiz for fastSnappyHexMesh, a faster version of
    snappyHexMesh. This file is not approved or endorsed by the OpenFOAM
    Foundation or OpenCFD Limited.

Class
    Foam::fsDuplicatePoints

Description
    Duplicate points

SourceFiles
    fsDuplicatePoints.C

\*---------------------------------------------------------------------------*/

#ifndef fsDuplicatePoints_H
#define fsDuplicatePoints_H

#include "labelList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class fsPolyTopoChange;
class edge;
class face;
class polyTopoChangeMap;
class localPointRegion;

/*---------------------------------------------------------------------------*\
                       Class fsDuplicatePoints Declaration
\*---------------------------------------------------------------------------*/

class fsDuplicatePoints
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Per regionSide.pointRegion() the points created.
        //  element 0 is the original point.
        //  Can have empty entries if a point only has one region.
        labelListList duplicates_;


public:

    //- Runtime type information
    ClassName("fsDuplicatePoints");


    // Constructors

        //- Construct from mesh
        fsDuplicatePoints(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        fsDuplicatePoints(const fsDuplicatePoints&) = delete;


    // Member Functions

        // Access

            //- Per point in regionSide.pointRegions() the points created.
            const labelListList& duplicates() const
            {
                return duplicates_;
            }


        // Topology changes

            //- Play commands into fsPolyTopoChange to duplicate points. Gets
            //  localPointRegion structure which is per non-manifold point
            //  the regions per point.
            void setRefinement
            (
                const localPointRegion& regionSide,
                fsPolyTopoChange&
            );

            //- Force recalculation of locally stored data on topological change
            void topoChange(const polyTopoChangeMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fsDuplicatePoints&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
