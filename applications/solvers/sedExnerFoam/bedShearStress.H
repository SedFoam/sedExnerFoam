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

}

// map volumic shields number to area shields number
shields.primitiveFieldRef() =
    bed.vsm.ref().mapToSurface<vector>
    (
        shieldsVf.boundaryFieldRef()
    );

if (filterShields)
{
    const projectedFaMesh& meshProj = bed.aProjMesh();

    if (projectedGlobalEdgeCount(meshProj) > 0)
    {
        const scalar dl = projectedMinEdgeDistance(meshProj);
        const scalar Kex = Foam::sqr(dl) * alphaShields;

        shields.primitiveFieldRef() = filterField
        (
            NShields,
            Kex,
            shields.primitiveField(),
            meshProj,
            meshProj.S(),
            meshProj.areaCentres(),
            meshProj.Le(),
            meshProj.edgeOwner(),
            meshProj.edgeNeighbour(),
            shields
        );
    }

    // Remove any bed-normal component introduced by the interpolation.
    const vectorField& faceNormals =
        bed.aMesh().faceAreaNormals().internalField();
    forAll(shields, facei)
    {
        shields[facei] -= (shields[facei] & faceNormals[facei])
            * faceNormals[facei];
    }

    shields.correctBoundaryConditions();
    bed.vsm.ref().mapToVolume
    (
        shields,
        shieldsVf.boundaryFieldRef()
    );
}

for (const label patchID : bed.bedPatchesID())
{
    const vectorField& ssp = shieldsVf.boundaryField()[patchID];
    Info << "Shields number: max:" << gMax(mag(ssp));
    Info << " min:" << gMin(mag(ssp)) << endl;
}
