// fastSnappyHexMesh: copy of OpenFOAM-13 removeCells.H, class renamed to fsRemoveCells and
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
    Foam::fsRemoveCells

Description
    Given list of cells to remove insert all the topology changes.

    Works in two passes:
    - get faces that will become boundary faces
    - given these faces and the patches they have to go into make the
      changes.

SourceFiles
    fsRemoveCells.C

\*---------------------------------------------------------------------------*/

#ifndef fsRemoveCells_H
#define fsRemoveCells_H

#include "labelList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class fsPolyTopoChange;
class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                         Class fsRemoveCells Declaration
\*---------------------------------------------------------------------------*/

class fsRemoveCells
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Whether or not to synchronise parallel case.
        const bool syncPar_;


    // Private Member Functions

        //- Decrease count of elements of f
        static void uncount
        (
            const labelList& f,
            labelList& nUsage
        );


public:

    //- Runtime type information
    ClassName("fsRemoveCells");


    // Constructors

        //- Construct from mesh. syncPar: do parallel synchronisation.
        fsRemoveCells(const polyMesh& mesh, const bool syncPar = true);


    // Member Functions

        // Topology changes

            //- Get labels of exposed faces.
            //  These are
            //  - internal faces that become boundary faces
            //  - coupled faces that become uncoupled (since one of the sides
            //    gets deleted)
            labelList getExposedFaces(const labelList& cellsToRemove) const;

            //- Play commands into fsPolyTopoChange to remove cells.
            //  patchIDs is for every element in facesToExpose (see above) the
            //  patch it has to go into. This cannot be a coupled patch!
            void setRefinement
            (
                const labelList& cellsToRemove,
                const labelList& facesToExpose,
                const labelList& patchIDs,
                fsPolyTopoChange&
            ) const;

            //- Force recalculation of locally stored data on topological change
            void topoChange(const polyTopoChangeMap&)
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
