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

// reset phiqbProj to 0
scalarField& phiqbProj = phiqbProjPtr.ref();
phiqbProj = Zero;
// reset dVsed to 0
areaScalarField& dVsed = dVsedPtr.ref();
// store old dVsed field
dVsed.storeOldTimes();
forAll(bed.aMesh().areaCentres(), facei)
{
    dVsed[facei] = 0;
}

// compute gradient of bed level, slope
areaVectorField& gradZb = gradZbPtr.ref();

// slope direction
const vectorField& nFaces =
    bed.aMesh().faceAreaNormals().internalField();
// bed faces angle
scalarField slopeAngle = Foam::acos(nFaces & eg);
// slope direction
vectorField slopeDir = eg - nFaces * (eg & nFaces);
slopeDir /= (Foam::mag(slopeDir) + SMALL);

gradZb.primitiveFieldRef() = -slopeAngle * slopeDir;

// compute bed angle beta
areaScalarField& beta = betaPtr.ref();
beta = Foam::mag(gradZb);

areaVectorField& shields = shieldsPtr.ref();

// access bedloadModel
const bedloadModels::bedloadModel& bedloadMod = bed.bedloadModel();

// compute wallShearStress on each patch belonging to sediment bed
for (const label patchID : bed.bedPatchesID())
{
    // compute wall shear stress on patch
    vectorField& ssp = shieldsVf.boundaryFieldRef()[patchID];
    const vectorField& Sfp = mesh.Sf().boundaryField()[patchID];
    const scalarField& magSfp =
        mesh.magSf().boundaryField()[patchID];
    volSymmTensorField Reff =
        rhoF * turbulence->nuEff() * dev(twoSymm(fvc::grad(U)));
    const symmTensorField& Reffp = Reff.boundaryField()[patchID];
    // compute shields number
    ssp = bedloadMod.coefShields()
        * ((-Sfp / magSfp) & Reffp)
        / (mag(g).value() * dS.value()
        * (rhoS.value() - rhoF.value()));

    // remove bed normal component, tangential component only
    ssp += (ssp & (-Sfp / magSfp)) * (Sfp / magSfp);

    if (filterShields)
    {

        // access to projected faMesh
        const projectedFaMesh& meshProj = bed.aProjMesh();

        // projected face centers coordinates
        const vectorField& xFacesProj = meshProj.areaCentres();
        // face owner list
        const labelList& owner = meshProj.edgeOwner();
        // face neighbour list
        const labelList& neighbour = meshProj.edgeNeighbour();

        const scalarField& SfProj = meshProj.S();
        const vectorField& LeProj = meshProj.Le();

        if (meshProj.nInternalEdges() > 0)
        {

            scalar dl = Foam::VGREAT;

            for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
            {
                scalar newdl = Foam::mag(
                xFacesProj[neighbour[edgei]] - xFacesProj[owner[edgei]]);

                if (newdl < dl)
                {
                    dl = newdl;
                }

            }

            scalar Kex = Foam::sqr(dl) * alphaShields;


                ssp = filterField(
                    NShields,
                    Kex,
                    ssp,
                    meshProj,
                    SfProj,
                    xFacesProj,
                    LeProj,
                    owner,
                    neighbour);

        }

        // remove bed normal component, tangential component only
        ssp += (ssp & (-Sfp / magSfp)) * (Sfp / magSfp);
    }

    // print min/max of Shields number
    Info << "Shields number: max:" << gMax(mag(ssp));
    Info << " min:" << gMin(mag(ssp)) << endl;
}

// map volumic shields number to area shields number
shields.primitiveFieldRef() =
    bed.vsm.ref().mapToSurface<vector>
    (
        shieldsVf.boundaryFieldRef()
    );
