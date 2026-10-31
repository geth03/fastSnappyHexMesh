// fastSnappyHexMesh: copy of OpenFOAM-13 polyTopoChangeI.H, class renamed to fsPolyTopoChange.
// Changes against the original are marked "fastSnappy".
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

\*---------------------------------------------------------------------------*/

#include "face.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::label Foam::fsPolyTopoChange::nPointsAll() const
{
    return nLazyPoints_ + points_.size();
}


inline Foam::point Foam::fsPolyTopoChange::pointAt(const label pointi) const
{
    if (pointi < nLazyPoints_)
    {
        Map<point>::const_iterator iter = lazyPointValue_.find(pointi);
        if (iter != lazyPointValue_.end())
        {
            return iter();
        }
        return (*lazyMeshPoints_)[pointi];
    }
    return points_[pointi - nLazyPoints_];
}


inline Foam::label Foam::fsPolyTopoChange::pointMapAt(const label pointi) const
{
    if (pointi < nLazyPoints_)
    {
        Map<label>::const_iterator iter = lazyPointMap_.find(pointi);
        return (iter != lazyPointMap_.end() ? iter() : pointi);
    }
    return pointMap_[pointi - nLazyPoints_];
}


inline bool Foam::fsPolyTopoChange::pointRemoved(const label pointi) const
{
    const point pt = pointAt(pointi);

    return
        pt.x() > 0.5*vector::max.x()
     && pt.y() > 0.5*vector::max.y()
     && pt.z() > 0.5*vector::max.z();
}


inline bool Foam::fsPolyTopoChange::faceRemoved(const label facei) const
{
    return
        faces_[facei].empty()
     && !(facei < lazyFace_.size() && lazyFace_.get(facei));
}


inline bool Foam::fsPolyTopoChange::cellRemoved(const label celli) const
{
    return cellMap_[celli] == -2;
}


inline void Foam::fsPolyTopoChange::setNumPatches(const label nPatches)
{
    nPatches_ = nPatches;
}


// ************************************************************************* //
