// fastSnappyHexMesh: copy of OpenFOAM-13 polyTopoChangeTemplates.C, class renamed to fsPolyTopoChange.
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

#include "fsPolyTopoChange.H"

#include <type_traits>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::fsPolyTopoChange::reorder
(
    const labelList& oldToNew,
    DynamicList<T>& lst
)
{
    // fastSnappy: identity map -> nothing to do (refinement and baffling
    // remove no faces/points, so the compaction maps are identities)
    {
        bool identity = (oldToNew.size() <= lst.size());
        forAll(oldToNew, elemI)
        {
            if (oldToNew[elemI] != elemI)
            {
                identity = false;
                break;
            }
        }
        if (identity)
        {
            return;
        }
    }

    // Plain data (labels, points, bools): copy + scatter as the original,
    // which is just memory copies
    if (std::is_trivially_copyable<T>::value)
    {
        const List<T> oldLst(lst);

        forAll(oldToNew, elemI)
        {
            label newElemI = oldToNew[elemI];

            if (newElemI != -1)
            {
                lst[newElemI] = oldLst[elemI];
            }
        }
        return;
    }

    // fastSnappy: move instead of copy. The original copied the whole list
    // and then copy-assigned every element - for faces (each a heap list)
    // two allocations + copies per face, twice per changeMesh. Elements not
    // targeted by the map keep their old value as before.
    DynamicList<T> oldLst;
    oldLst.transfer(lst);
    lst.setSize(oldLst.size());

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        if (newElemI != -1)
        {
            lst[newElemI] = std::move(oldLst[elemI]);
        }
    }

    // Untargeted slots: restore the old element of the same index (the
    // original left the old value there). Moved-from slots are empty.
    labelList isTarget(lst.size(), 0);
    forAll(oldToNew, elemI)
    {
        if (oldToNew[elemI] != -1)
        {
            isTarget[oldToNew[elemI]] = 1;
        }
    }
    forAll(lst, i)
    {
        if (!isTarget[i])
        {
            lst[i] = std::move(oldLst[i]);
        }
    }
}


template<class T>
void Foam::fsPolyTopoChange::reorder
(
    const labelList& oldToNew,
    List<DynamicList<T>>& lst
)
{
    // Create copy
    List<DynamicList<T>> oldLst(lst);

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        if (newElemI != -1)
        {
            lst[newElemI].transfer(oldLst[elemI]);
        }
    }
}


template<class T>
void Foam::fsPolyTopoChange::renumberKey
(
    const labelList& oldToNew,
    Map<T>& elems
)
{
    Map<T> newElems(elems.size());

    forAllConstIter(typename Map<T>, elems, iter)
    {
        label newElem = oldToNew[iter.key()];

        if (newElem >= 0)
        {
            newElems.insert(newElem, iter());
        }
    }

    elems.transfer(newElems);
}


// ************************************************************************* //
