/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "motionSmootherAlgo.H"
#include "twoDPointCorrector.H"
#include "faceSet.H"
#include "pointSet.H"
#include "fixedValuePointPatchFields.H"
#include "pointConstraints.H"
#include "syncTools.H"
#include "meshTools.H"
#include "meshCheck.H"
#include "fsMeshCheck.H"
#include "OFstream.H"
#include "fsThreads.H"
#include "fsGeometry.H"
#include "fsProfiler.H"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSmootherAlgo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::motionSmootherAlgo::testSyncPositions
(
    const pointField& fld,
    const scalar maxMag
) const
{
    pointField syncedFld(fld);

    syncTools::syncPointPositions
    (
        mesh_,
        syncedFld,
        minEqOp<point>(),           // combine op
        point(great,great,great)    // null
    );

    forAll(syncedFld, i)
    {
        if (mag(syncedFld[i] - fld[i]) > maxMag)
        {
            FatalErrorInFunction
                << "On point " << i << " point:" << fld[i]
                << " synchronised point:" << syncedFld[i]
                << abort(FatalError);
        }
    }
}


void Foam::motionSmootherAlgo::checkFld(const pointScalarField& fld)
{
    forAll(fld, pointi)
    {
        const scalar val = fld[pointi];

        if ((val > -great) && (val < great))
        {}
        else
        {
            FatalErrorInFunction
                << "Problem : point:" << pointi << " value:" << val
                << abort(FatalError);
        }
    }
}


Foam::labelHashSet Foam::motionSmootherAlgo::getPoints
(
    const labelHashSet& faceLabels
) const
{
    labelHashSet usedPoints(mesh_.nPoints()/100);

    forAllConstIter(labelHashSet, faceLabels, iter)
    {
        const face& f = mesh_.faces()[iter.key()];

        forAll(f, fp)
        {
            usedPoints.insert(f[fp]);
        }
    }

    return usedPoints;
}


Foam::tmp<Foam::scalarField> Foam::motionSmootherAlgo::calcEdgeWeights
(
    const pointField& points
) const
{
    FS_PROF("calcEdgeWeights");
    const edgeList& edges = mesh_.edges();

    tmp<scalarField> twght(new scalarField(edges.size()));
    scalarField& wght = twght.ref();

    forAll(edges, edgeI)
    {
        wght[edgeI] = 1.0/(edges[edgeI].mag(points)+small);
    }
    return twght;
}


namespace
{

// fastSnappyHexMesh: minSmooth only reads avg() at the affected points. In
// serial without constraint (incl. coupled) patches the syncs and
// constrain() in avg() leave a scalar field unchanged, so avg() can be
// evaluated at those points only: the master edges of a point from
// pointEdges (ascending, the order in which avg() adds them).
bool pointAvgUsable(const Foam::polyMesh& mesh)
{
    using namespace Foam;

    if (Pstream::parRun())
    {
        return false;
    }

    const char* e = std::getenv("FS_FASTMINSMOOTH");
    if (e && std::string(e) == "0")
    {
        return false;
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchi)
    {
        if (polyPatch::constraintType(patches[patchi].type()))
        {
            return false;
        }
    }

    return true;
}


// avg(fld, edgeWeights) at pointi, as motionSmootherAlgo::avg
Foam::scalar pointAvg
(
    const Foam::polyMesh& mesh,
    const Foam::PackedBoolList& isMasterEdge,
    const Foam::scalarField& edgeWeights,
    const Foam::pointScalarField& fld,
    const Foam::label pointi
)
{
    using namespace Foam;

    const edgeList& edges = mesh.edges();
    const labelList& pEdges = mesh.pointEdges()[pointi];

    scalar s = 0;
    scalar sumWeight = 0;

    forAll(pEdges, i)
    {
        const label edgeI = pEdges[i];

        if (isMasterEdge.get(edgeI) == 1)
        {
            const scalar w = edgeWeights[edgeI];

            s += w*fld[edges[edgeI].otherVertex(pointi)];
            sumWeight += w;
        }
    }

    if (mag(sumWeight) < vSmall)
    {
        // Unconnected point. Take over original value
        return fld[pointi];
    }
    else
    {
        s /= sumWeight;
        return s;
    }
}

} // End anonymous namespace


// Smooth on selected points (usually patch points)
void Foam::motionSmootherAlgo::minSmooth
(
    const scalarField& edgeWeights,
    const PackedBoolList& isAffectedPoint,
    const labelList& meshPoints,
    const pointScalarField& fld,
    pointScalarField& newFld
) const
{
    FS_PROF("minSmooth");

    if (pointAvgUsable(mesh_))
    {
        forAll(meshPoints, i)
        {
            label pointi = meshPoints[i];
            if (isAffectedPoint.get(pointi) == 1)
            {
                newFld[pointi] = min
                (
                    fld[pointi],
                    0.5*fld[pointi]
                  + 0.5*pointAvg(mesh_, isMasterEdge_, edgeWeights, fld, pointi)
                );
            }
        }

        // Single and multi-patch constraints
        pointConstraints::New(pMesh()).constrain(newFld, false);
        return;
    }

    tmp<pointScalarField> tavgFld = avg
    (
        fld,
        edgeWeights // scalarField(mesh_.nEdges(), 1.0)    // uniform weighting
    );
    const pointScalarField& avgFld = tavgFld();

    forAll(meshPoints, i)
    {
        label pointi = meshPoints[i];
        if (isAffectedPoint.get(pointi) == 1)
        {
            newFld[pointi] = min
            (
                fld[pointi],
                0.5*fld[pointi] + 0.5*avgFld[pointi]
            );
        }
    }

    // Single and multi-patch constraints
    pointConstraints::New(pMesh()).constrain(newFld, false);
}


// Smooth on all internal points
void Foam::motionSmootherAlgo::minSmooth
(
    const scalarField& edgeWeights,
    const PackedBoolList& isAffectedPoint,
    const pointScalarField& fld,
    pointScalarField& newFld
) const
{
    FS_PROF("minSmooth");

    if (pointAvgUsable(mesh_))
    {
        forAll(fld, pointi)
        {
            if (isAffectedPoint.get(pointi) == 1 && isInternalPoint(pointi))
            {
                newFld[pointi] = min
                (
                    fld[pointi],
                    0.5*fld[pointi]
                  + 0.5*pointAvg(mesh_, isMasterEdge_, edgeWeights, fld, pointi)
                );
            }
        }

        // Single and multi-patch constraints
        pointConstraints::New(pMesh()).constrain(newFld, false);
        return;
    }

    tmp<pointScalarField> tavgFld = avg
    (
        fld,
        edgeWeights // scalarField(mesh_.nEdges(), 1.0)    // uniform weighting
    );
    const pointScalarField& avgFld = tavgFld();

    forAll(fld, pointi)
    {
        if (isAffectedPoint.get(pointi) == 1 && isInternalPoint(pointi))
        {
            newFld[pointi] = min
            (
                fld[pointi],
                0.5*fld[pointi] + 0.5*avgFld[pointi]
            );
        }
    }

   // Single and multi-patch constraints
    pointConstraints::New(pMesh()).constrain(newFld, false);

}


// Scale on all internal points
void Foam::motionSmootherAlgo::scaleField
(
    const labelHashSet& pointLabels,
    const scalar scale,
    pointScalarField& fld
) const
{
    forAllConstIter(labelHashSet, pointLabels, iter)
    {
        if (isInternalPoint(iter.key()))
        {
            fld[iter.key()] *= scale;
        }
    }

    // Single and multi-patch constraints
    pointConstraints::New(pMesh()).constrain(fld, false);
}


// Scale on selected points (usually patch points)
void Foam::motionSmootherAlgo::scaleField
(
    const labelList& meshPoints,
    const labelHashSet& pointLabels,
    const scalar scale,
    pointScalarField& fld
) const
{
    forAll(meshPoints, i)
    {
        label pointi = meshPoints[i];

        if (pointLabels.found(pointi))
        {
            fld[pointi] *= scale;
        }
    }
}


// Lower on internal points
void Foam::motionSmootherAlgo::subtractField
(
    const labelHashSet& pointLabels,
    const scalar f,
    pointScalarField& fld
) const
{
    forAllConstIter(labelHashSet, pointLabels, iter)
    {
        if (isInternalPoint(iter.key()))
        {
            fld[iter.key()] = max(0.0, fld[iter.key()]-f);
        }
    }

    // Single and multi-patch constraints
    pointConstraints::New(pMesh()).constrain(fld);
}


// Scale on selected points (usually patch points)
void Foam::motionSmootherAlgo::subtractField
(
    const labelList& meshPoints,
    const labelHashSet& pointLabels,
    const scalar f,
    pointScalarField& fld
) const
{
    forAll(meshPoints, i)
    {
        label pointi = meshPoints[i];

        if (pointLabels.found(pointi))
        {
            fld[pointi] = max(0.0, fld[pointi]-f);
        }
    }
}


bool Foam::motionSmootherAlgo::isInternalPoint(const label pointi) const
{
    return isInternalPoint_.get(pointi) == 1;
}


void Foam::motionSmootherAlgo::getAffectedFacesAndPoints
(
    const label nPointIter,
    const faceSet& wrongFaces,

    labelList& affectedFaces,
    PackedBoolList& isAffectedPoint
) const
{
    FS_PROF("getAffectedFacesAndPoints");
    isAffectedPoint.setSize(mesh_.nPoints());
    isAffectedPoint = 0;

    faceSet nbrFaces(mesh_, "checkFaces", wrongFaces);

    // Find possible points influenced by nPointIter iterations of
    // scaling and smoothing by doing pointCellpoint walk.
    // Also update faces-to-be-checked to extend one layer beyond the points
    // that will get updated.

    for (label i = 0; i < nPointIter; i++)
    {
        pointSet nbrPoints(mesh_, "grownPoints", getPoints(nbrFaces.toc()));

        forAllConstIter(pointSet, nbrPoints, iter)
        {
            const labelList& pCells = mesh_.pointCells(iter.key());

            forAll(pCells, pCelli)
            {
                const cell& cFaces = mesh_.cells()[pCells[pCelli]];

                forAll(cFaces, cFacei)
                {
                    nbrFaces.insert(cFaces[cFacei]);
                }
            }
        }
        nbrFaces.sync(mesh_);

        if (i == nPointIter - 2)
        {
            forAllConstIter(faceSet, nbrFaces, iter)
            {
                const face& f = mesh_.faces()[iter.key()];
                forAll(f, fp)
                {
                    isAffectedPoint.set(f[fp], 1);
                }
            }
        }
    }

    affectedFaces = nbrFaces.toc();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSmootherAlgo::motionSmootherAlgo
(
    polyMesh& mesh,
    pointMesh& pMesh,
    indirectPrimitivePatch& pp,
    pointVectorField& displacement,
    pointScalarField& scale,
    pointField& oldPoints,
    const labelList& adaptPatchIDs,
    const dictionary& paramDict
)
:
    mesh_(mesh),
    pMesh_(pMesh),
    pp_(pp),
    displacement_(displacement),
    scale_(scale),
    oldPoints_(oldPoints),
    adaptPatchIndices_(adaptPatchIDs),
    paramDict_(paramDict),
    isInternalPoint_(mesh_.nPoints(), 1)
{
    topoChange();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSmootherAlgo::~motionSmootherAlgo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::motionSmootherAlgo::mesh() const
{
    return mesh_;
}


const Foam::pointMesh& Foam::motionSmootherAlgo::pMesh() const
{
    return pMesh_;
}


const Foam::indirectPrimitivePatch& Foam::motionSmootherAlgo::patch() const
{
    return pp_;
}


const Foam::labelList& Foam::motionSmootherAlgo::adaptPatchIndices() const
{
    return adaptPatchIndices_;
}


const Foam::dictionary& Foam::motionSmootherAlgo::paramDict() const
{
    return paramDict_;
}


void Foam::motionSmootherAlgo::correct()
{
    oldPoints_ = mesh_.points();

    scale_ = 1.0;

    // No need to update twoDmotion corrector since only holds edge labels
    // which will remain the same as before. So unless the mesh was distorted
    // severely outside of motionSmootherAlgo there will be no need.
}


void Foam::motionSmootherAlgo::setDisplacementPatchFields
(
    const labelList& patchIDs,
    pointVectorField& displacement
)
{
    pointVectorField::Boundary& displacementBf =
        displacement.boundaryFieldRef();

    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches)
    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];

        displacementBf[patchi] ==
            displacementBf[patchi].patchInternalField();
    }

    // Make consistent with non-adapted bc's by evaluating those now and
    // resetting the displacement from the values.
    // Note that we're just doing a correctBoundaryConditions with
    // fixedValue bc's first.
    labelHashSet adaptPatchSet(patchIDs);

    const lduSchedule& patchSchedule = displacement.mesh().globalData().
        patchSchedule();

    forAll(patchSchedule, patchEvalI)
    {
        label patchi = patchSchedule[patchEvalI].patch;

        if (!adaptPatchSet.found(patchi))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacementBf[patchi]
                    .initEvaluate(Pstream::commsTypes::scheduled);
            }
            else
            {
                displacementBf[patchi]
                    .evaluate(Pstream::commsTypes::scheduled);
            }
        }
    }

    // Multi-patch constraints
    pointConstraints::New(displacement.mesh()).constrainCorners(displacement);

    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches) to take the changes caused
    // by multi-corner constraints into account.
    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];

        displacementBf[patchi] ==
            displacementBf[patchi].patchInternalField();
    }
}


void Foam::motionSmootherAlgo::setDisplacementPatchFields()
{
    setDisplacementPatchFields(adaptPatchIndices_, displacement_);
}


void Foam::motionSmootherAlgo::setDisplacement
(
    const labelList& patchIDs,
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    pointVectorField& displacement
)
{
    const polyMesh& mesh = displacement.mesh()();

    // See comment in .H file about shared points.
    // We want to disallow effect of loose coupled points - we only
    // want to see effect of proper fixedValue boundary conditions

    const labelList& cppMeshPoints =
        mesh.globalData().coupledPatch().meshPoints();

    forAll(cppMeshPoints, i)
    {
        displacement[cppMeshPoints[i]] = vector::zero;
    }

    const labelList& ppMeshPoints = pp.meshPoints();

    // Set internal point data from displacement on combined patch points.
    forAll(ppMeshPoints, patchPointi)
    {
        displacement[ppMeshPoints[patchPointi]] = patchDisp[patchPointi];
    }


    // Combine any coupled points
    syncTools::syncPointList
    (
        mesh,
        displacement,
        maxMagEqOp(),           // combine op
        vector::zero            // null value
    );


    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches)
    setDisplacementPatchFields(patchIDs, displacement);


    if (debug)
    {
        OFstream str(mesh.db().path()/"changedPoints.obj");
        label nVerts = 0;
        forAll(ppMeshPoints, patchPointi)
        {
            const vector& newDisp = displacement[ppMeshPoints[patchPointi]];

            if (mag(newDisp-patchDisp[patchPointi]) > small)
            {
                const point& pt = mesh.points()[ppMeshPoints[patchPointi]];

                meshTools::writeOBJ(str, pt);
                nVerts++;
                // Pout<< "Point:" << pt
                //    << " oldDisp:" << patchDisp[patchPointi]
                //    << " newDisp:" << newDisp << endl;
            }
        }
        Pout<< "Written " << nVerts << " points that are changed to file "
            << str.name() << endl;
    }

    // Now reset input displacement
    forAll(ppMeshPoints, patchPointi)
    {
        patchDisp[patchPointi] = displacement[ppMeshPoints[patchPointi]];
    }
}


void Foam::motionSmootherAlgo::setDisplacement(pointField& patchDisp)
{
    setDisplacement(adaptPatchIndices_, pp_, patchDisp, displacement_);
}


// correctBoundaryConditions with fixedValue bc's first.
void Foam::motionSmootherAlgo::correctBoundaryConditions
(
    pointVectorField& displacement
) const
{
    labelHashSet adaptPatchSet(adaptPatchIndices_);

    const lduSchedule& patchSchedule = mesh_.globalData().patchSchedule();

    pointVectorField::Boundary& displacementBf =
        displacement.boundaryFieldRef();

    // 1. evaluate on adaptPatches
    forAll(patchSchedule, patchEvalI)
    {
        label patchi = patchSchedule[patchEvalI].patch;

        if (adaptPatchSet.found(patchi))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacementBf[patchi]
                    .initEvaluate(Pstream::commsTypes::scheduled);
            }
            else
            {
                displacementBf[patchi]
                    .evaluate(Pstream::commsTypes::scheduled);
            }
        }
    }


    // 2. evaluate on non-AdaptPatches
    forAll(patchSchedule, patchEvalI)
    {
        label patchi = patchSchedule[patchEvalI].patch;

        if (!adaptPatchSet.found(patchi))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacementBf[patchi]
                    .initEvaluate(Pstream::commsTypes::scheduled);
            }
            else
            {
                displacementBf[patchi]
                    .evaluate(Pstream::commsTypes::scheduled);
            }
        }
    }

    // Multi-patch constraints
    pointConstraints::New(displacement.mesh()).constrainCorners(displacement);

    // Correct for problems introduced by corner constraints
    syncTools::syncPointList
    (
        mesh_,
        displacement,
        maxMagEqOp(),           // combine op
        vector::zero            // null value
    );
}



void Foam::motionSmootherAlgo::modifyMotionPoints(pointField& newPoints) const
{
    // Correct for 2-D motion
    const twoDPointCorrector& twoDCorrector = twoDPointCorrector::New(mesh_);

    if (twoDCorrector.required())
    {
        Info<< "Correcting 2-D mesh motion";

        if (mesh_.globalData().parallel())
        {
            WarningInFunction
                << "2D mesh-motion probably not correct in parallel" << endl;
        }

        // We do not want to move 3D planes so project all points onto those
        const pointField& oldPoints = mesh_.points();
        const edgeList& edges = mesh_.edges();

        const labelList& neIndices = twoDCorrector.normalEdgeIndices();
        const vector& pn = twoDCorrector.planeNormal();

        forAll(neIndices, i)
        {
            const edge& e = edges[neIndices[i]];

            point& pStart = newPoints[e.start()];

            pStart += pn*(pn & (oldPoints[e.start()] - pStart));

            point& pEnd = newPoints[e.end()];

            pEnd += pn*(pn & (oldPoints[e.end()] - pEnd));
        }

        // Correct tangentially
        twoDCorrector.correctPoints(newPoints);
        Info<< " ...done" << endl;
    }

    if (debug)
    {
        Pout<< "motionSmootherAlgo::modifyMotionPoints :"
            << " testing sync of newPoints."
            << endl;
        testSyncPositions(newPoints, 1e-6*mesh_.bounds().mag());
    }
}


void Foam::motionSmootherAlgo::movePoints()
{
    // Make sure to clear out tetPtIs since used in checks (sometimes, should
    // really check)
    mesh_.clearTetBasePtIs();
    pp_.clearGeom();
}


Foam::scalar Foam::motionSmootherAlgo::setErrorReduction
(
    const scalar errorReduction
)
{
    scalar oldErrorReduction = paramDict_.lookup<scalar>("errorReduction");

    paramDict_.remove("errorReduction");
    paramDict_.add("errorReduction", errorReduction);

    return oldErrorReduction;
}


bool Foam::motionSmootherAlgo::scaleMesh
(
    const labelList& checkFaces,
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    List<labelPair> emptyBaffles;
    return scaleMesh
    (
        checkFaces,
        emptyBaffles,
        smoothMesh,
        nAllowableErrors
    );
}


bool Foam::motionSmootherAlgo::scaleMesh
(
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    return scaleMesh
    (
        checkFaces,
        baffles,
        paramDict_,
        paramDict_,
        smoothMesh,
        nAllowableErrors
    );
}


Foam::tmp<Foam::pointField> Foam::motionSmootherAlgo::curPoints() const
{
    // Set newPoints as old + scale*displacement

    // Create overall displacement with same b.c.s as displacement_
    wordList actualPatchTypes;
    {
        const pointBoundaryMesh& pbm = displacement_.mesh().boundary();
        actualPatchTypes.setSize(pbm.size());
        forAll(pbm, patchi)
        {
            actualPatchTypes[patchi] = pbm[patchi].type();
        }
    }

    wordList actualPatchFieldTypes;
    {
        const pointVectorField::Boundary& pfld =
            displacement_.boundaryField();
        actualPatchFieldTypes.setSize(pfld.size());
        forAll(pfld, patchi)
        {
            if (isA<fixedValuePointPatchField<vector>>(pfld[patchi]))
            {
                // Get rid of funny
                actualPatchFieldTypes[patchi] =
                    fixedValuePointPatchField<vector>::typeName;
            }
            else
            {
                actualPatchFieldTypes[patchi] = pfld[patchi].type();
            }
        }
    }

    pointVectorField totalDisplacement
    (
        IOobject
        (
            "totalDisplacement",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        scale_*displacement_,
        actualPatchFieldTypes,
        actualPatchTypes
    );
    correctBoundaryConditions(totalDisplacement);

    if (debug)
    {
        Pout<< "scaleMesh : testing sync of totalDisplacement" << endl;
        testSyncField
        (
            totalDisplacement,
            maxMagEqOp(),
            vector::zero,   // null value
            1e-6*mesh_.bounds().mag()
        );
    }

    tmp<pointField> tnewPoints(oldPoints_ + totalDisplacement.primitiveField());

    // Correct for 2-D motion
    modifyMotionPoints(tnewPoints.ref());

    return tnewPoints;
}


bool Foam::motionSmootherAlgo::scaleMesh
(
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    const dictionary& paramDict,
    const dictionary& meshQualityDict,
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    if (!smoothMesh && adaptPatchIndices_.empty())
    {
        FatalErrorInFunction
            << "You specified both no movement on the internal mesh points"
            << " (smoothMesh = false)" << nl
            << "and no movement on the patch (adaptPatchIDs is empty)" << nl
            << "Hence nothing to adapt."
            << exit(FatalError);
    }

    if (debug)
    {
        // Had a problem with patches moved non-synced. Check transformations.
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        Pout<< "Entering scaleMesh : coupled patches:" << endl;
        forAll(patches, patchi)
        {
            if (patches[patchi].coupled())
            {
                const coupledPolyPatch& pp =
                    refCast<const coupledPolyPatch>(patches[patchi]);

                Pout<< '\t' << patchi << '\t' << pp.name()
                    << " transform:" << pp.transform()
                    << endl;
            }
        }
    }

    const scalar errorReduction =
        paramDict.lookup<scalar>("errorReduction");

    const label nSmoothScale =
        paramDict.lookup<label>("nSmoothScale");


    // Note: displacement_ should already be synced already from setDisplacement
    // but just to make sure.
    syncTools::syncPointList
    (
        mesh_,
        displacement_,
        maxMagEqOp(),
        vector::zero    // null value
    );

    Info<< "Moving mesh using displacement scaling :"
        << " min:" << gMin(scale_.primitiveField())
        << "  max:" << gMax(scale_.primitiveField())
        << endl;

    {
        FS_PROF("setPoints");

        // Get points using current displacement and scale. Optionally 2D
        // corrected.
        pointField newPoints(curPoints());

        // fastSnappy: keep the geometry of the last check for the update
        incrCheck_.beforeMove(mesh_);

        // Move. No need to do 2D correction since curPoints already done
        // that.
        mesh_.setPoints(newPoints);
        movePoints();
    }

    {
        FS_PROF("geometry");

        // Threaded, bit-identical geometry for the checks below; only the
        // faces and cells near moved points when possible
        incrCheck_.afterMove(mesh_);
    }


    // Check. Returns parallel number of incorrect faces.
    faceSet wrongFaces(mesh_, "wrongFaces", mesh_.nFaces()/100+100);
    {
        FS_PROF("checkMesh");
        incrCheck_.check
        (
            mesh_,
            meshQualityDict,
            checkFaces,
            baffles,
            wrongFaces
        );
    }

    if (returnReduce(wrongFaces.size(), sumOp<label>()) <= nAllowableErrors)
    {
        return true;
    }
    else
    {
        // Sync across coupled faces by extending the set.
        wrongFaces.sync(mesh_);

        // Special case:
        // if errorReduction is set to zero, extend wrongFaces
        // to face-Cell-faces to ensure quick return to previously valid mesh

        if (mag(errorReduction) < small)
        {
            labelHashSet newWrongFaces(wrongFaces);
            forAllConstIter(labelHashSet, wrongFaces, iter)
            {
                label own = mesh_.faceOwner()[iter.key()];
                const cell& ownFaces = mesh_.cells()[own];

                forAll(ownFaces, cfI)
                {
                    newWrongFaces.insert(ownFaces[cfI]);
                }

                if (iter.key() < mesh_.nInternalFaces())
                {
                    label nei = mesh_.faceNeighbour()[iter.key()];
                    const cell& neiFaces = mesh_.cells()[nei];

                    forAll(neiFaces, cfI)
                    {
                        newWrongFaces.insert(neiFaces[cfI]);
                    }
                }
            }
            wrongFaces.transfer(newWrongFaces);
            wrongFaces.sync(mesh_);
        }


        // Find out points used by wrong faces and scale displacement.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        pointSet usedPoints(mesh_, "usedPoints", getPoints(wrongFaces));
        usedPoints.sync(mesh_);



        // Grow a few layers to determine
        // - points to be smoothed
        // - faces to be checked in next iteration
        PackedBoolList isAffectedPoint(mesh_.nPoints());
        labelList affectedFaces(checkFaces);
        getAffectedFacesAndPoints
        (
            nSmoothScale,       // smoothing iterations
            wrongFaces,         // error faces
            affectedFaces,
            isAffectedPoint
        );

        if (debug)
        {
            Pout<< "Faces in error:" << wrongFaces.size()
                << "  with points:" << usedPoints.size()
                << endl;
        }

        if (adaptPatchIndices_.size())
        {
            // Scale conflicting patch points
            scaleField(pp_.meshPoints(), usedPoints, errorReduction, scale_);
            // subtractField(pp_.meshPoints(), usedPoints, 0.1, scale_);
        }
        if (smoothMesh)
        {
            // Scale conflicting internal points
            scaleField(usedPoints, errorReduction, scale_);
            // subtractField(usedPoints, 0.1, scale_);
        }

        scalarField eWeights(calcEdgeWeights(oldPoints_));

        for (label i = 0; i < nSmoothScale; i++)
        {
            if (adaptPatchIndices_.size())
            {
                // Smooth patch values
                pointScalarField oldScale(scale_);
                minSmooth
                (
                    eWeights,
                    isAffectedPoint,
                    pp_.meshPoints(),
                    oldScale,
                    scale_
                );
                checkFld(scale_);
            }
            if (smoothMesh)
            {
                // Smooth internal values
                pointScalarField oldScale(scale_);
                minSmooth(eWeights, isAffectedPoint, oldScale, scale_);
                checkFld(scale_);
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            scale_,
            maxEqOp<scalar>(),
            -great              // null value
        );


        if (debug)
        {
            Pout<< "scale_ after smoothing :"
                << " min:" << Foam::gMin(scale_)
                << " max:" << Foam::gMax(scale_)
                << endl;
        }

        return false;
    }
}


void Foam::motionSmootherAlgo::topoChange()
{
    incrCheck_.clear();

    const pointBoundaryMesh& patches = pMesh_.boundary();

    // Check whether displacement has fixed value b.c. on adaptPatchID
    forAll(adaptPatchIndices_, i)
    {
        label patchi = adaptPatchIndices_[i];

        if
        (
           !isA<fixedValuePointPatchVectorField>
            (
                displacement_.boundaryField()[patchi]
            )
        )
        {
            FatalErrorInFunction
                << "Patch " << patches[patchi].name()
                << " has wrong boundary condition "
                << displacement_.boundaryField()[patchi].type()
                << " on field " << displacement_.name() << nl
                << "Only type allowed is "
                << fixedValuePointPatchVectorField::typeName
                << exit(FatalError);
        }
    }


    // Determine internal points. Note that for twoD there are no internal
    // points so we use the points of adaptPatchIDs instead

    const labelList& meshPoints = pp_.meshPoints();

    forAll(meshPoints, i)
    {
        isInternalPoint_.unset(meshPoints[i]);
    }

    // Calculate master edge addressing
    isMasterEdge_ = syncTools::getMasterEdges(mesh_);
}


void Foam::motionSmootherAlgo::smoothIterations
(
    const scalarField& edgeWeight,
    const label nIter,
    pointVectorField& fld
) const
{
    const polyMesh& mesh = fld.mesh()();
    const label nPoints = mesh.nPoints();
    const edgeList& edges = mesh.edges();

    // Point-edge list over the master edges. Filled in edge order, so the
    // entries of every point are in ascending edge order, which is the
    // order in which avg() adds them.
    labelList offset(nPoints + 1, 0);
    forAll(edges, edgeI)
    {
        if (isMasterEdge_.get(edgeI) == 1)
        {
            offset[edges[edgeI][0] + 1]++;
            offset[edges[edgeI][1] + 1]++;
        }
    }
    for (label pointi = 0; pointi < nPoints; pointi++)
    {
        offset[pointi + 1] += offset[pointi];
    }

    labelList nbr(offset[nPoints]);
    scalarField wt(offset[nPoints]);
    {
        labelList fill(SubList<label>(offset, nPoints));
        forAll(edges, edgeI)
        {
            if (isMasterEdge_.get(edgeI) == 1)
            {
                const edge& e = edges[edgeI];
                const scalar w = edgeWeight[edgeI];

                label k = fill[e[0]]++;
                nbr[k] = e[1];
                wt[k] = w;

                k = fill[e[1]]++;
                nbr[k] = e[0];
                wt[k] = w;
            }
        }
    }

    // Interior points (on no point patch) renumbered in Morton (z-order)
    // order into the interior slots, so the neighbour reads of the gather
    // hit the cache. Patch points keep their index: constrain() and the
    // parallel syncs only touch patch points, so they work on the
    // renumbered values as they are. Every point keeps its own entries in
    // the same order, so the sums are the same. FS_SMOOTHORDER=0: off.
    labelList newIdx;
    {
        const char* e = std::getenv("FS_SMOOTHORDER");
        if (!(e && std::string(e) == "0"))
        {
            boolList onPatch(nPoints, false);
            const pointBoundaryMesh& pbm = fld.mesh().boundary();
            forAll(pbm, patchi)
            {
                UIndirectList<bool>(onPatch, pbm[patchi].meshPoints()) =
                    true;
            }
            const labelList interior(findIndices(onPatch, false));
            onPatch.clear();

            // Morton keys from 21 bits per direction
            auto spread = [](uint64_t x)
            {
                x &= 0x1fffff;
                x = (x | x << 32) & 0x1f00000000ffffULL;
                x = (x | x << 16) & 0x1f0000ff0000ffULL;
                x = (x | x << 8) & 0x100f00f00f00f00fULL;
                x = (x | x << 4) & 0x10c30c30c30c30c3ULL;
                x = (x | x << 2) & 0x1249249249249249ULL;
                return x;
            };

            const pointField& pts = mesh.points();
            const boundBox bb(pts, false);
            const vector span(bb.span());
            const scalar maxQ = 2097151.0;

            std::vector<std::pair<uint64_t, label>> keys(interior.size());
            fsThreads::forRange
            (
                interior.size(),
                [&](const label i0, const label i1)
                {
                    for (label i = i0; i < i1; i++)
                    {
                        const vector d(pts[interior[i]] - bb.min());
                        uint64_t q[3];
                        for (direction cmpt = 0; cmpt < 3; cmpt++)
                        {
                            const scalar s =
                                maxQ*d[cmpt]/(span[cmpt] + vSmall);
                            q[cmpt] = uint64_t(min(max(s, 0.0), maxQ));
                        }
                        keys[i] = std::make_pair
                        (
                            spread(q[0])
                          | spread(q[1]) << 1
                          | spread(q[2]) << 2,
                            interior[i]
                        );
                    }
                }
            );
            std::sort(keys.begin(), keys.end());

            newIdx = identityMap(nPoints);
            forAll(interior, i)
            {
                newIdx[keys[i].second] = interior[i];
            }
        }
    }

    vectorField& nf = fld.primitiveFieldRef();
    List<char> internal(nPoints);

    if (newIdx.size())
    {
        labelList oldIdx(nPoints);
        forAll(newIdx, pointi)
        {
            oldIdx[newIdx[pointi]] = pointi;
        }

        labelList offsetN(nPoints + 1);
        offsetN[0] = 0;
        forAll(oldIdx, i)
        {
            const label o = oldIdx[i];
            offsetN[i+1] = offsetN[i] + offset[o+1] - offset[o];
        }
        labelList nbrN(nbr.size());
        scalarField wtN(wt.size());
        fsThreads::forRange
        (
            nPoints,
            [&](const label i0, const label i1)
            {
                for (label i = i0; i < i1; i++)
                {
                    const label o = oldIdx[i];
                    label k = offsetN[i];
                    for (label j = offset[o]; j < offset[o+1]; j++)
                    {
                        nbrN[k] = newIdx[nbr[j]];
                        wtN[k] = wt[j];
                        k++;
                    }
                    internal[i] = isInternalPoint(o);
                }
            }
        );
        offset.transfer(offsetN);
        nbr.transfer(nbrN);
        wt.transfer(wtN);

        const vectorField old(nf);
        forAll(old, pointi)
        {
            nf[newIdx[pointi]] = old[pointi];
        }
    }
    else
    {
        forAll(internal, pointi)
        {
            internal[pointi] = isInternalPoint(pointi);
        }
    }

    const bool parRun = Pstream::parRun();
    const pointConstraints& pcs = pointConstraints::New(fld.mesh());

    // Reused for all iterations. A fresh avg() field differs only in its
    // values, which are all rewritten below.
    tmp<pointVectorField> tres
    (
        pointVectorField::New
        (
            "avg("+fld.name()+')',
            fld.mesh(),
            dimensioned<vector>("zero", fld.dimensions(), Zero)
        )
    );
    pointVectorField& res = tres.ref();
    vectorField& r = res.primitiveFieldRef();

    scalarField sumWeight(nPoints);

    for (label iter = 0; iter < nIter; iter++)
    {
        if ((iter % 10) == 0)
        {
            Info<< "Iteration " << iter << endl;
        }

        // fld still holds the old values until the update below, so it
        // serves as the old field of smooth()
        const vectorField& f = fld.primitiveField();

        {
            FS_PROF("gather");

            fsThreads::forRange
            (
                nPoints,
                [&](const label p0, const label p1)
                {
                    for (label pointi = p0; pointi < p1; pointi++)
                    {
                        vector s(Zero);
                        scalar sw = 0.0;
                        for
                        (
                            label k = offset[pointi];
                            k < offset[pointi+1];
                            k++
                        )
                        {
                            s += wt[k]*f[nbr[k]];
                            sw += wt[k];
                        }

                        if (parRun)
                        {
                            r[pointi] = s;
                            sumWeight[pointi] = sw;
                        }
                        else if (mag(sw) < vSmall)
                        {
                            r[pointi] = f[pointi];
                        }
                        else
                        {
                            s /= sw;
                            r[pointi] = s;
                        }
                    }
                },
                label(20000)
            );
        }

        if (parRun)
        {
            syncTools::syncPointList
            (
                mesh,
                r,
                plusEqOp<vector>(),
                vector::zero
            );
            syncTools::syncPointList
            (
                mesh,
                sumWeight,
                plusEqOp<scalar>(),
                scalar(0)
            );

            fsThreads::forRange
            (
                nPoints,
                [&](const label p0, const label p1)
                {
                    for (label pointi = p0; pointi < p1; pointi++)
                    {
                        if (mag(sumWeight[pointi]) < vSmall)
                        {
                            r[pointi] = f[pointi];
                        }
                        else
                        {
                            r[pointi] /= sumWeight[pointi];
                        }
                    }
                },
                label(20000)
            );
        }

        pcs.constrain(res, false);

        // Each point reads only its own old value, so the update can be in
        // place
        {
            FS_PROF("update");

            fsThreads::forRange
            (
                nPoints,
                [&](const label p0, const label p1)
                {
                    for (label pointi = p0; pointi < p1; pointi++)
                    {
                        if (internal[pointi])
                        {
                            nf[pointi] = 0.5*nf[pointi] + 0.5*r[pointi];
                        }
                    }
                },
                label(20000)
            );
        }

        {
            FS_PROF("constrain");
            pcs.constrain(fld, false);
        }
    }

    if (newIdx.size())
    {
        const vectorField cur(nf);
        forAll(newIdx, pointi)
        {
            nf[pointi] = cur[newIdx[pointi]];
        }
    }
}


// ************************************************************************* //
