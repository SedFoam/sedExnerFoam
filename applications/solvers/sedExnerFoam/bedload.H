/*---------------------------------------------------------------------------*\
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of ScourFOAM.

    ScourFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ScourFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ScourFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// access to saturated bedload and bedload fields
areaVectorField& qsat = qsatPtr.ref();
areaVectorField& qb = qbPtr.ref();

dimensionedScalar smallVal(dimless, SMALL);
// unit vector gravity oriented
vector eg = (g/mag(g)).value();

//Info<<"Critical shear stress calculation"<<endl;
// access criticalShieldsModel
const criticalShieldsModels::criticalShieldsModel& critShieldsMod =
    bed.critShieldsModel();

// access critical Shields number
areaScalarField& critShields = critShieldsPtr.ref();

if (not readCritShields)
{
    // shear stress direction
    vectorField stressDir = shields / (Foam::mag(shields) + SMALL);
    // slope direction
    const vectorField& normalFaces =
        bed.aMesh().faceAreaNormals().internalField();
    // bed faces angle
    scalarField slopeAngle = Foam::acos(normalFaces & eg);
    // slope direction
    vectorField slopeDir = eg - normalFaces * (eg & normalFaces);
    slopeDir /= (Foam::mag(slopeDir) + SMALL);
    //slopeDir /= (Foam::mag(slopeDir) + SMALL);
    // Apply slope correction on critical Shields number
    critShieldsMod.slopeCorrection(
        critShields, stressDir, slopeAngle, slopeDir, betaRep);
}

//Info<<"qsat calculation"<<endl;
// compute bedload from Shields number
qsat = bedloadMod.qb(
    shields, critShields, gradZb, rhoS, rhoF, mag(g), dS, nuF);

qsat.correctBoundaryConditions();
// Refresh coupled values before using qb as the state in the saturation
// transport equation.
qb.correctBoundaryConditions();

// test 3D
/*
qsat.correctBoundaryConditions();
// access to projected faMesh
const projectedFaMesh& meshProj = bed.aProjMesh();
forAll(meshProj.S(), facei)
{
    qsat[facei].y() = 0;
    qsat[facei].z() = 0;
}
*/
// apply saturation time to compute bedload
if (saturationType == "implicit1D")
{
    Info << "Implicit saturation" << endl;

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "saturation type 'implicit1D' solves a dense system using "
            << "only the local finite-area face ordering. It is not "
            << "distributed across processor boundaries. Use 'explicit' "
            << "or 'semiImplicit' for parallel runs."
            << exit(FatalError);
    }

    // access to projected faMesh
    const projectedFaMesh& meshProj = bed.aProjMesh();
    // projected face centers coordinates
    const vectorField& xFacesProj = meshProj.areaCentres();
    // face owner list
    const labelList& owner = meshProj.edgeOwner();
    // face neighbour list
    const labelList& neighbour = meshProj.edgeNeighbour();
    label Nf = meshProj.nFaces();  // number of faces
    Info << "n faces = " << Nf << endl;

    if (Nf < 2)
    {
        FatalErrorInFunction
            << "saturation type 'implicit1D' requires at least two "
            << "finite-area faces."
            << exit(FatalError);
    }

    scalarSquareMatrix Amat(meshProj.nFaces());
    Amat = Zero;

    scalar dl = 0;
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        scalar newdl = Foam::mag(
            xFacesProj[neighbour[edgei]] - xFacesProj[owner[edgei]]);
        if (newdl>dl)
        {
            dl = newdl;
        }
    }

    if (dl <= VSMALL)
    {
        FatalErrorInFunction
            << "saturation type 'implicit1D' requires a non-zero "
            << "internal-edge distance."
            << exit(FatalError);
    }

    scalar Mcoef = Lsat->value() / (2*dl);
    scalarField qsx(Nf);

    for (label i=1; i<Nf-1; i++)
    {
        Amat[i][i] = 1;
        Amat[i][i+1] = Mcoef;
        Amat[i][i-1] = -Mcoef;
        qsx[i] = qsat[i].x();
    }
    Amat[0][0] = 1 - Mcoef;
    Amat[0][1] = Mcoef;
    Amat[Nf-1][Nf-1] = 1 + Mcoef;
    Amat[Nf-1][Nf-2] = -Mcoef;
    qsx[0] = qsat[0].x();
    qsx[Nf-1] = qsat[Nf-1].x();

    LUscalarMatrix LU(Amat);
    scalarField qbx(LU.solve(qsx));

    forAll(meshProj.S(), facei)
    {
        qb[facei].x() = qbx[facei];
    }
}
else if (saturationType == "explicit")
{
    Info << "Explicit saturation" << endl;
    // time step
    const scalar& deltaT = runTime.time().deltaT().value();
    scalar TsatLim = 2 * deltaT;
    scalar TsatCorr = max(Tsat->value(), TsatLim);
    if (Tsat->value() < TsatLim)
    {
        Info << "WARNING: Tsat < TsatLim, new Tsat = " << TsatCorr << endl;
    }
    // access to projected faMesh
    const projectedFaMesh& meshProj = bed.aProjMesh();
    // projected face centers coordinates
    const vectorField& xFacesProj = meshProj.areaCentres();
    // projected edge centers coordinates
    const vectorField& xEdgesProj = meshProj.edgeCentres();
    // projected edge vectors, oriented toward owner face
    const vectorField& LeProj = meshProj.Le();
    //Info << "Le: " << LeProj << endl;
    // area of projected mesh faces
    const scalarField& SfProj = meshProj.S();
    // face owner list
    const labelList& owner = meshProj.edgeOwner();
    // face neighbour list
    const labelList& neighbour = meshProj.edgeNeighbour();

    // project qb on horizontal plane
    const vectorField qbProjOld = meshProj.project(qb);

    // project qsat on horizontal plane
    const vectorField qsatProj = meshProj.project(qsat);

    // new projected bedload values
    vectorField qbProj = meshProj.project(qb);

    // saturation equation, internal edges contribution
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        // distance between edgei center and owner face center
        scalar distOwner =
            Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
        // distance between edgei center and neighbour center
        scalar distNeighbour =
            Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);
        // linear interpolate qb on edge i (centered scheme)
        vector qbEdgei =
            (
                distNeighbour * qbProjOld[owner[edgei]]
                + distOwner * qbProjOld[neighbour[edgei]]
            )
            * (1 / (distOwner + distNeighbour));
        // upwind interpolate qb on edge i
        scalar faceFlux = qbEdgei & LeProj[edgei];
        if (faceFlux >= 0)
        {
            qbEdgei = qbProjOld[owner[edgei]];
        }
        else
        {
            qbEdgei = qbProjOld[neighbour[edgei]];
        }
        vector divLsQb = (deltaT * Lsat->value() / TsatCorr)
            * (qbEdgei & LeProj[edgei])
            * (qbEdgei / (Foam::mag(qbEdgei) + SMALL));

        qbProj[owner[edgei]] -= divLsQb / SfProj[owner[edgei]];

        qbProj[neighbour[edgei]] += divLsQb / SfProj[neighbour[edgei]];
    }
    qbProj += (deltaT / TsatCorr) * (qsatProj - qbProjOld);

    // Saturation transport across coupled edges. The same upwind edge value
    // is used on both sides of a processor interface, giving equal and
    // opposite contributions to the two local faces.
    forAll(bed.aMesh().boundary(), patchi)
    {
        const faPatch& patch = bed.aMesh().boundary()[patchi];

        if (not patch.coupled())
        {
            continue;
        }

        const faPatchField<Vector<scalar>>& QbBound =
            qb.boundaryField()[patchi];

        label edgei = patch.start();
        forAll(patch, patchEdgei)
        {
            const label own = owner[edgei];
            const vector qbOwner = qbProjOld[own];
            const vector qbNeighbour = meshProj.project(QbBound[patchEdgei]);

            vector qbEdgei = 0.5*(qbOwner + qbNeighbour);
            if ((qbEdgei & LeProj[edgei]) >= scalar(0))
            {
                qbEdgei = qbOwner;
            }
            else
            {
                qbEdgei = qbNeighbour;
            }

            const vector divLsQb = (deltaT * Lsat->value() / TsatCorr)
                * (qbEdgei & LeProj[edgei])
                * (qbEdgei / (Foam::mag(qbEdgei) + SMALL));

            qbProj[own] -= divLsQb / SfProj[own];
            edgei++;
        }
    }

    // normal to faMesh faces
    const vectorField& nFaces =
        bed.aMesh().faceAreaNormals().internalField();

    // reconstruct qb from qbProj on bed surface
    forAll(meshProj.S(), facei)
    {
        /*
        const vector& nf = nFaces[facei];
        vector qbReproj = qbProj[facei] - nf * (qbProj[facei] & nf);
        qbReproj *= Foam::mag(qbProj[facei]) / (Foam::mag(qbReproj) + SMALL);
            qb[facei] = qbReproj;
        */

        qb[facei] = qbProj[facei];
    }

}
else if (saturationType == "semiImplicit")
{
    const scalar& deltaT = runTime.time().deltaT().value();

    // access to projected faMesh
    const projectedFaMesh& meshProj = bed.aProjMesh();
    // projected face centers coordinates
    const vectorField& xFacesProj = meshProj.areaCentres();
    // projected edge centers coordinates
    const vectorField& xEdgesProj = meshProj.edgeCentres();
    // projected edge vectors, oriented toward owner face
    const vectorField& LeProj = meshProj.Le();
    // area of projected mesh faces
    const scalarField& SfProj = meshProj.S();
    // face owner list
    const labelList& owner = meshProj.edgeOwner();
    // face neighbour list
    const labelList& neighbour = meshProj.edgeNeighbour();

    // project qb on horizontal plane
    const vectorField qbProjOld = meshProj.project(qb);

    // project qsat on horizontal plane
    const vectorField qsatProj = meshProj.project(qsat);

    // new projected bedload values
    vectorField qbProj = meshProj.project(qb);

    scalar localMaxCFL = -1;

    // saturation equation, internal edges contribution
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {

        scalar localCFL = Lsat->value()/Tsat->value() * deltaT
            * Foam::mag(LeProj[edgei]) / SfProj[owner[edgei]];

        localMaxCFL = max(localMaxCFL, localCFL);

        // distance between edgei center and owner face center
        scalar distOwner =
            Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
        // distance between edgei center and neighbour center
        scalar distNeighbour =
            Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);
        // linear interpolate qb on edge i
        vector qbEdgei =
            (
                distNeighbour * qbProjOld[owner[edgei]]
                + distOwner * qbProjOld[neighbour[edgei]]
            )
            * (1 / (distOwner + distNeighbour));
        scalar faceFlux = qbEdgei & LeProj[edgei];
        if (faceFlux >= 0)
        {
            qbEdgei = qbProjOld[owner[edgei]];
        }
        else
        {
            qbEdgei = qbProjOld[neighbour[edgei]];
        }
        vector divLsQb = (Lsat->value()*deltaT / (deltaT + Tsat->value()))
            * (qbEdgei & LeProj[edgei])
            * (qbEdgei / (Foam::mag(qbEdgei) + SMALL));

        qbProj[owner[edgei]] -= divLsQb / SfProj[owner[edgei]];

        qbProj[neighbour[edgei]] += divLsQb / SfProj[neighbour[edgei]];
    }
    qbProj += ((1/(1 + deltaT / Tsat->value())) - 1)
        * qbProjOld + (deltaT / (deltaT + Tsat->value())) * qsatProj;

    // Saturation equation on boundary edges. Coupled patches use the same
    // upwind flux on both sides; physical patches keep their prescribed
    // boundary value.
    forAll(bed.aMesh().boundary(), patchi)
    {
        const faPatch& patch = bed.aMesh().boundary()[patchi];
        const faPatchField<Vector<scalar>>& QbBound =
            qb.boundaryField()[patchi];

        if (patch.type() != "empty")
        {
            label edgei = patch.start();
            forAll(patch, patchEdgei)
            {
                const label own = owner[edgei];

                const scalar localCFL = Lsat->value()/Tsat->value()*deltaT
                    * Foam::mag(LeProj[edgei])/SfProj[own];
                localMaxCFL = max(localMaxCFL, localCFL);

                vector qbEdgei = meshProj.project(QbBound[patchEdgei]);

                if (patch.coupled())
                {
                    const vector qbOwner = qbProjOld[own];
                    const vector qbNeighbour = qbEdgei;
                    qbEdgei = 0.5*(qbOwner + qbNeighbour);

                    if ((qbEdgei & LeProj[edgei]) >= scalar(0))
                    {
                        qbEdgei = qbOwner;
                    }
                    else
                    {
                        qbEdgei = qbNeighbour;
                    }
                }

                const vector divLsQb = (Lsat->value()*deltaT
                    / (deltaT + Tsat->value()))
                    * (qbEdgei & LeProj[edgei])
                    * (qbEdgei / (Foam::mag(qbEdgei) + SMALL));

                qbProj[own] -= divLsQb / SfProj[own];
                edgei++;
            }
        }
    }

    // normal to faMesh faces
    const vectorField& nFaces =
        bed.aMesh().faceAreaNormals().internalField();

    // reconstruct qb from qbProj on bed surface
    forAll(meshProj.S(), facei)
    {
        /*const vector& nf = nFaces[facei];
        vector qbReproj = qbProj[facei] - nf * (qbProj[facei] & nf);
        qbReproj *= Foam::mag(qbProj[facei]) / (Foam::mag(qbReproj) + SMALL);
            qb[facei] = qbReproj;*/


        qb[facei] = qbProj[facei];
    }

    reduce(localMaxCFL, maxOp<scalar>());

    Info << "max saturation CFL = " << localMaxCFL << endl;

}
else if (saturationType == "none")
{
    Info << "No saturation" << endl;
    qb = qsat;
}
else
{
    FatalError << "in bedloadProperties/saturation,"
        << "wrong keyword for entry type "
        << saturationType
        << "\nvalid options are: implicit1D explicit semiImplicit none" << endl;
    Info << abort(FatalError) << endl;
}

qb.correctBoundaryConditions();

Info << "Avalanche calculation" << endl;
// additional bedload due to avalanche
if (bedloadMod.avalanche())
{
    #include "sedimAvalanche.H"
}
