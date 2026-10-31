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
    Foam::fsWave

Description
    Threaded, exact replica of
        PointEdgeWave<Type, TrackingData> wave(mesh, seedPoints, seedInfo,
            allPointInfo, allEdgeInfo, 0, td);
        wave.iterate(maxIter);
    for serial meshes without coupled patches. Type's own equal(),
    updateEdge() and updatePoint() are used, so the update rules are the
    original ones.

    The wave's result depends on the order of the updates, so the order is
    kept: a pass over the changed points updates every edge with at most
    two events (one per changed end point) - the earlier one in a first
    parallel phase, the later one in a second (each chunk queues its second
    events in order). A pass over the changed edges updates every point from
    its changed edges in list order, on the thread of the first of them
    (claimed with an atomic minimum of the list positions). The next change
    lists are rebuilt chunk by chunk in the serial order. Every point and
    edge sees the same sequence of updates as in PointEdgeWave.

    usable() is false in parallel, with coupled patches (the processor,
    cyclic and collocated handling is not replicated), with one thread, or
    with FS_FASTWAVE=0; then use PointEdgeWave.

\*---------------------------------------------------------------------------*/

#ifndef fsPointEdgeWave_H
#define fsPointEdgeWave_H

#include "polyMesh.H"
#include "PointEdgeWave.H"
#include "fsThreads.H"

#include <atomic>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

namespace Foam
{
namespace fsWave
{

//- See description
bool usable(const polyMesh& mesh);

//- Number of chunks for n items
inline label nChunks(const label n)
{
    return fsThreads::nChunks(n, 4096);
}

using fsThreads::forChunks;


//- See description. Returns the number of iterations, as iterate().
template<class Type, class TrackingData>
label pointEdgeWave
(
    const polyMesh& mesh,
    const labelUList& seedPoints,
    const UList<Type>& seedInfo,
    UList<Type>& allPointInfo,
    UList<Type>& allEdgeInfo,
    const label maxIter,
    TrackingData& td
)
{
    const edgeList& edges = mesh.edges();
    const labelListList& pointEdges = mesh.pointEdges();
    mesh.points();
    const label nPoints = mesh.nPoints();
    const label nEdges = mesh.nEdges();

    const scalar tol = PointEdgeWave<Type, TrackingData>::propagationTol();

    // Changed points and edges in wave order, and their list positions
    DynamicList<label> chP(nPoints);
    labelList posP(nPoints, -1);
    DynamicList<label> chE(nEdges);
    labelList posE(nEdges, -1);

    // setPointInfo
    forAll(seedPoints, i)
    {
        const label pointi = seedPoints[i];

        allPointInfo[pointi] = seedInfo[i];

        if (posP[pointi] < 0)
        {
            posP[pointi] = chP.size();
            chP.append(pointi);
        }
    }

    // Pass in which the first/second event of an edge propagated; pass and
    // key (2*list position + end index) of the first propagating event of
    // a point
    labelList prop1(nEdges, -1);
    labelList prop2(nEdges, -1);
    labelList firstPass(nPoints, -1);
    labelList firstKey(nPoints, -1);
    label pass = 0;

    List<DynamicList<label>> parts;
    List<DynamicList<labelPair>> secs;

    // List position of the first changed edge of a point, labelMax: none
    std::vector<std::atomic<label>> owner(nPoints);
    for (std::atomic<label>& o : owner)
    {
        o.store(labelMax, std::memory_order_relaxed);
    }

    label iter = 0;

    while (iter < maxIter)
    {
        // pointToEdge
        // ~~~~~~~~~~~
        pass++;
        {
            const label m = chP.size();
            const label nt = nChunks(m);

            // Events of the changed points on their edges, in list order.
            // The event is the second of its edge when the other end point
            // is earlier in the list. Per chunk: the edges in event order
            // (second events as -(edge+1)) and the second events
            // (edge, point).
            parts.setSize(nt);
            secs.setSize(nt);
            forChunks
            (
                m,
                nt,
                [&](const label t, const label i0, const label i1)
                {
                    DynamicList<label>& cand = parts[t];
                    DynamicList<labelPair>& sec = secs[t];
                    cand.clear();
                    sec.clear();

                    for (label i = i0; i < i1; i++)
                    {
                        const label pointi = chP[i];
                        const Type& pInfo = allPointInfo[pointi];
                        const labelList& pEdges = pointEdges[pointi];

                        forAll(pEdges, k)
                        {
                            const label edgei = pEdges[k];
                            const label q = edges[edgei].otherVertex(pointi);

                            if (posP[q] >= 0 && posP[q] < i)
                            {
                                sec.append(labelPair(edgei, pointi));
                                cand.append(-(edgei + 1));
                            }
                            else
                            {
                                Type& eInfo = allEdgeInfo[edgei];

                                if
                                (
                                    !eInfo.equal(pInfo, td)
                                 && eInfo.updateEdge
                                    (
                                        mesh,
                                        edgei,
                                        pointi,
                                        pInfo,
                                        tol,
                                        td
                                    )
                                )
                                {
                                    prop1[edgei] = pass;
                                    cand.append(edgei);
                                }
                            }
                        }
                    }
                }
            );

            // Second events, then the changed edges of the chunk in the
            // order of their first propagating event: a first event that
            // propagated, or a second one after a first one that did not
            forChunks
            (
                m,
                nt,
                [&](const label t, const label, const label)
                {
                    const DynamicList<labelPair>& sec = secs[t];

                    forAll(sec, s)
                    {
                        const label edgei = sec[s].first();
                        const label pointi = sec[s].second();
                        const Type& pInfo = allPointInfo[pointi];
                        Type& eInfo = allEdgeInfo[edgei];

                        if
                        (
                            !eInfo.equal(pInfo, td)
                         && eInfo.updateEdge
                            (
                                mesh,
                                edgei,
                                pointi,
                                pInfo,
                                tol,
                                td
                            )
                        )
                        {
                            prop2[edgei] = pass;
                        }
                    }

                    DynamicList<label>& cand = parts[t];
                    label n = 0;
                    forAll(cand, c)
                    {
                        const label x = cand[c];

                        if (x >= 0)
                        {
                            cand[n++] = x;
                        }
                        else
                        {
                            const label edgei = -x - 1;

                            if (prop2[edgei] == pass && prop1[edgei] != pass)
                            {
                                cand[n++] = edgei;
                            }
                        }
                    }
                    cand.setSize(n);
                }
            );

            forAll(chP, i)
            {
                posP[chP[i]] = -1;
            }
            chP.clear();

            chE.clear();
            for (label t = 0; t < nt; t++)
            {
                chE.append(parts[t]);
            }
            forAll(chE, j)
            {
                posE[chE[j]] = j;
            }
        }

        if (chE.empty())
        {
            break;
        }


        // edgeToPoint
        // ~~~~~~~~~~~
        pass++;
        {
            const label n = chE.size();
            const label nt = nChunks(n);

            // Every point takes its events from its changed edges in list
            // order, on the thread of the first of them
            forChunks
            (
                n,
                nt,
                [&](const label, const label j0, const label j1)
                {
                    for (label j = j0; j < j1; j++)
                    {
                        const edge& e = edges[chE[j]];

                        forAll(e, eI)
                        {
                            std::atomic<label>& o = owner[e[eI]];
                            label cur = o.load(std::memory_order_relaxed);
                            while
                            (
                                j < cur
                             && !o.compare_exchange_weak
                                (
                                    cur,
                                    j,
                                    std::memory_order_relaxed
                                )
                            )
                            {}
                        }
                    }
                }
            );

            forChunks
            (
                n,
                nt,
                [&](const label, const label j0, const label j1)
                {
                    DynamicList<labelPair> events(16);

                    for (label j = j0; j < j1; j++)
                    {
                        const edge& e = edges[chE[j]];

                        forAll(e, eI)
                        {
                            const label pointi = e[eI];

                            if
                            (
                                owner[pointi].load(std::memory_order_relaxed)
                             != j
                            )
                            {
                                continue;
                            }

                            const labelList& pEdges = pointEdges[pointi];

                            events.clear();
                            forAll(pEdges, k)
                            {
                                const label pos = posE[pEdges[k]];
                                if (pos >= 0)
                                {
                                    events.append(labelPair(pos, pEdges[k]));
                                }
                            }

                            // Insertion sort on the list position
                            for (label a = 1; a < events.size(); a++)
                            {
                                const labelPair ev = events[a];
                                label b = a - 1;
                                while
                                (
                                    b >= 0
                                 && events[b].first() > ev.first()
                                )
                                {
                                    events[b+1] = events[b];
                                    b--;
                                }
                                events[b+1] = ev;
                            }

                            Type& pInfo = allPointInfo[pointi];
                            label key = -1;

                            forAll(events, evi)
                            {
                                const label edgei = events[evi].second();
                                const Type& eInfo = allEdgeInfo[edgei];

                                if (!pInfo.equal(eInfo, td))
                                {
                                    const bool propagate = pInfo.updatePoint
                                    (
                                        mesh,
                                        pointi,
                                        edgei,
                                        eInfo,
                                        tol,
                                        td
                                    );

                                    if (propagate && key < 0)
                                    {
                                        key =
                                            2*events[evi].first()
                                          + (edges[edgei][0] == pointi ? 0 : 1);
                                    }
                                }
                            }

                            firstPass[pointi] = pass;
                            firstKey[pointi] = key;
                        }
                    }
                }
            );

            // Changed points in the order of their first propagating event;
            // release the points
            parts.setSize(nt);
            forChunks
            (
                n,
                nt,
                [&](const label t, const label j0, const label j1)
                {
                    DynamicList<label>& part = parts[t];
                    part.clear();

                    for (label j = j0; j < j1; j++)
                    {
                        const edge& e = edges[chE[j]];

                        forAll(e, eI)
                        {
                            const label pointi = e[eI];

                            owner[pointi].store
                            (
                                labelMax,
                                std::memory_order_relaxed
                            );

                            if
                            (
                                firstPass[pointi] == pass
                             && firstKey[pointi] == 2*j + eI
                            )
                            {
                                part.append(pointi);
                            }
                        }
                    }
                }
            );

            forAll(chE, j)
            {
                posE[chE[j]] = -1;
            }
            chE.clear();

            chP.clear();
            for (label t = 0; t < nt; t++)
            {
                chP.append(parts[t]);
            }
            forAll(chP, i)
            {
                posP[chP[i]] = i;
            }
        }

        if (chP.empty())
        {
            break;
        }

        iter++;
    }

    return iter;
}


//- Number of points the wave did not reach (PointEdgeWave::getUnsetPoints
//  for work arrays that start unset)
template<class Type, class TrackingData>
label nUnset(const UList<Type>& allPointInfo, TrackingData& td)
{
    label n = 0;
    forAll(allPointInfo, pointi)
    {
        if (!allPointInfo[pointi].valid(td))
        {
            n++;
        }
    }
    return n;
}


//- PointEdgeWave(mesh, seedPoints, seedInfo, allPointInfo, allEdgeInfo, 0,
//  td).iterate(maxIter), threaded when usable(). FS_WAVECHECK=1: run both
//  and report differing point/edge data.
template<class Type, class TrackingData>
void run
(
    const polyMesh& mesh,
    const labelList& seedPoints,
    const List<Type>& seedInfo,
    List<Type>& allPointInfo,
    List<Type>& allEdgeInfo,
    const label maxIter,
    TrackingData& td
)
{
    if (!usable(mesh))
    {
        PointEdgeWave<Type, TrackingData> wave
        (
            mesh,
            seedPoints,
            seedInfo,
            allPointInfo,
            allEdgeInfo,
            0,
            td
        );
        wave.iterate(maxIter);
        return;
    }

    const char* e = std::getenv("FS_WAVECHECK");
    if (e && std::string(e) == "1")
    {
        List<Type> pointRef(allPointInfo);
        List<Type> edgeRef(allEdgeInfo);
        {
            PointEdgeWave<Type, TrackingData> wave
            (
                mesh,
                seedPoints,
                seedInfo,
                pointRef,
                edgeRef,
                0,
                td
            );
            wave.iterate(maxIter);
        }
        pointEdgeWave
        (
            mesh,
            seedPoints,
            seedInfo,
            allPointInfo,
            allEdgeInfo,
            maxIter,
            td
        );

        label nDiff = 0;
        forAll(pointRef, i)
        {
            if (!(pointRef[i] == allPointInfo[i]))
            {
                nDiff++;
            }
        }
        forAll(edgeRef, i)
        {
            if (!(edgeRef[i] == allEdgeInfo[i]))
            {
                nDiff++;
            }
        }
        Info<< "fsWave check: " << nDiff << " of "
            << pointRef.size() + edgeRef.size()
            << " point/edge values differ" << endl;
        return;
    }

    pointEdgeWave
    (
        mesh,
        seedPoints,
        seedInfo,
        allPointInfo,
        allEdgeInfo,
        maxIter,
        td
    );
}

} // End namespace fsWave
} // End namespace Foam

#endif

// ************************************************************************* //
