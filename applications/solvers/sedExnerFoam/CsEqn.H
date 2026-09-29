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

// Access variables relative to suspension
volScalarField& Cs = CsPtr.ref();  // suspended load volumic fraction
volVectorField& Ws = WsPtr.ref();  // settling velocity
volScalarField& diffSed = diffSedPtr.ref();  // turbulent diffusivity
surfaceScalarField& phip = phipPtr.ref();  // flux(U+Ws)
surfaceScalarField& bedPhi = bedPhiPtr.ref();  // mesh flux through faces
// reference to erosion terms
volScalarField& Erosion = erosionPtr.ref();

// turbulent sediment diffusivity
if (CsDiffusivity=="turbulent")
{
    diffSed = turbulence->nut()/Sc;
}
else if (CsDiffusivity=="effective")
{
    diffSed = turbulence->nuEff()/Sc;
}


if (bed.exist())
{
    // reference to faMesh
    const faMesh& aMesh = bed.aMesh();
    // Shields number
    areaVectorField& shields = shieldsPtr.ref();
    // critical Shields stress
    areaScalarField& critShields = critShieldsPtr.ref();

    if (isWallDiff)
    {
        // additional near bed patches turbulent diffusivity for suspension
        // access wallDist object for distance to bed patches
        wallDist& bedDist = bedDistPtr.ref();
        diffSed += epsWall0 * 0.5
           *(
                1 - Foam::tanh(betaWall * (bedDist.y() - ks) / ks)
            );
    }
    diffSed.correctBoundaryConditions();

    // handle erosion and deposition if bedload and suspension are switched on
    for (const label patchID : bed.bedPatchesID())
    {
        // erosion boundary field
        scalarField& Ep = Erosion.boundaryFieldRef()[patchID];
        // reference to bed fvPatch
        const fvPatch& bedPatch = mesh.boundary()[patchID];

        // if dynamicMesh, meshPhi is needed for erosion flux correction
        if (mesh.changing())
        {
            // get meshPhi on bed patch
            bedPhi.boundaryFieldRef()[patchID] =
                mesh.setPhi().ref().boundaryFieldRef()[patchID];
        }
        const scalarField& bedPhiPatch = bedPhi.boundaryFieldRef()[patchID];
        // reference level
        scalar zbRef = 2.5 * dS.value();

        forAll(bedPatch, facei)
        {
            // label of neighbour cell
            label idNCell = bedPatch.faceCells()[facei];
            // shear friction dimension less parameter
            scalar Tshear =
                Foam::mag(shields[facei])/(critShields[facei] + SMALL) - 1;
            // compute reference concentration, Van Rijn 1984
            scalar cbRef =
                Foam::min(
                    0.015 * (dS.value()/zbRef)
                    * Foam::pow(Tshear * Foam::pos(Tshear), 1.5)
                    * Foam::pow(Dstar.value(), -0.3),
                    CsbMax.value()
                );
            // coordinates of cell centers located above bed
            const vector xCell = mesh.C()[idNCell];
            // face coordinates
            const vector xFace = aMesh.areaCentres()[facei];
            const vector nFace = aMesh.faceAreaNormals()[facei];
            const scalar Sface = aMesh.S()[facei];
            const scalar Vcell = mesh.V()[idNCell];
            scalar zi = -(xCell - xFace) & nFace;
            //scalar ci = Cs[idNCell];
            vector wsi = Ws[idNCell];
            scalar bedPhiFace = bedPhiPatch[facei];
            scalar Roi = Foam::mag(wsi) / Foam::max(diffSed[idNCell], 1e-6);
            //scalar cb = ci + cbRef * (
            //    1 - Foam::exp(-Roi * (Foam::mag(zi-zbRef))));
            // implicit, zeroGradient bc for Ws and Cs on bed needed
            // extrapolate cbRef at level
            scalar cbRefZ1 = Foam::min
                (
                    cbRef * Foam::exp(-Roi * (zi-zbRef)),
                    CsbMax.value()
                );
            Ep[facei] = (cbRefZ1 / Vcell)
                * ((wsi & nFace) * Sface - bedPhiFace);
        }
        // version proposée par Matthias mais non testée
        if (bed.rigidBed())
        {
            areaScalarField& HsedBed = HsedBedPtr.ref();

            // current time step
            scalar deltaT = runTime.time().deltaT().value();

            forAll(bedPatch, facei)
            {
                // label de la cellule voisine
                label idNCell = bedPatch.faceCells()[facei];
                const scalar Sface = aMesh.S()[facei];
                const scalar Vcell = mesh.V()[idNCell];
                
                if (HsedBed[facei]<=0)
                {
                    Ep[facei] = 0;
                }
                else
                {
                    Ep[facei] = Foam::min
                        (
                            Ep[facei],
                            (HsedBed[facei] * Sface) / (deltaT * Vcell)
                        );
                }
            }
        }
        /* version code Matthias original
        if (bed.rigidBed())
        {
            areaScalarField& HsedBed = HsedBedPtr.ref();

            // current time step
            scalar deltaT = runTime.time().deltaT().value();

            forAll(bedPatch, facei)
            {
                Ep[facei] = Foam::min
                    (
                        Ep[facei],
                        HsedBed[facei] / deltaT
                    );
            }
        }
         */
        forAll(bedPatch, facei)
        {
            // label of neighbour cell
            label idNCell = bedPatch.faceCells()[facei];
            Erosion[idNCell] = Ep[facei];
        }
    }
}

// MULES FORMULATION FOR SUSPENTION TRANSPORT
if (switchMULES=="on")
{
    word scheme("div(phi,Cs)");
    word schemer("div(phir,Cs)");
    surfaceScalarField phiCs
    (
        fvc::flux
        (
            phi,
            Cs,
            scheme
        )
        + fvc::flux
        (
            fvc::flux(Ws),
            Cs,
            schemer
        )
    );

    MULES::explicitSolve
    (
        geometricOneField(),
        Cs,
        phi,
        phiCs,
        zeroField(),
        zeroField(),
        UniformField<scalar>(0.99*0.7),
        zeroField()
    );

    // Apply the diffusion term separately to allow implicit solution
    // and boundedness of the explicit advection

    fvScalarMatrix CsEqn
    (
        fvm::ddt(Cs) - fvc::ddt(Cs)
      - fvm::laplacian(diffSed, Cs)
        ==
        fvOptions(Cs)
    );

    fvOptions.constrain(CsEqn);

    CsEqn.solve(mesh.solver("CsDiffusion"));

    phiCs += CsEqn.flux();

    Cs.correctBoundaryConditions();

    fvOptions.correct(Cs);

    // Add Boussinesq term
    // moved to Seqn rhok = 1 + Cs * (rhoS - rhoF) / rhoF ;
}

// CLASSICAL FORMULATION FOR SUSPENTION TRANSPORT
else if (switchMULES=="off")
{
    fvScalarMatrix CsEqn
    (
        fvm::ddt(Cs)
      + fvm::div(phip, Cs)
      - fvm::laplacian(diffSed, Cs)
        ==
        fvOptions(Cs)
        + Erosion
    );

    fvOptions.constrain(CsEqn);

    CsEqn.solve();

    Cs.correctBoundaryConditions();

    fvOptions.correct(Cs);
}

Ws = settlingModelPtr->Ufall(
    Cs, CsMax, dS, rhoS, rhoF, nuF, mag(g))*g/mag(g);

Ws.correctBoundaryConditions();

phip = fvc::flux(U + Ws);

// Make the  flux relative to the mesh motion
fvc::makeRelative(phip, U);

Info << "sediment phase volume fraction = "
        << Cs.weightedAverage(mesh.V()).value()
        << "  Min(Cs) = " << gMin(Cs)
        << "  Max(Cs) = " << gMax(Cs)
        << endl;
