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

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading/calculating face flux field phi\n" << endl;

#include "createPhi.H"

#include "createRDeltaT.H"

#include "readGravitationalAcceleration.H"

vector eg = (g / Foam::mag(g)).value();
Info << "unit vector oriented with gradvity: eg = " << eg << endl;

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// fluid and sediment properties
// fluid kinematic viscosity
dimensionedScalar nuF("nu", transportProperties);
// water density (kg.m-3), default is 1000 (water)
dimensionedScalar* rhoFptr = new dimensionedScalar(dimDensity, 1000);
if (transportProperties.found("rhoF"))
{
    rhoFptr = new dimensionedScalar("rhoF", dimDensity, transportProperties);
}
dimensionedScalar& rhoF = *rhoFptr;
// sediment density (kg.m-3)
dimensionedScalar rhoS("rhoS", dimDensity, transportProperties);
// density ratio
dimensionedScalar sRatio = rhoS / rhoF;
// sediment diameter (m)
dimensionedScalar dS("dS", dimLength, transportProperties);
// dimension less sediment diameter
dimensionedScalar Dstar = dS
    * Foam::pow(
        (sRatio - 1) * mag(g) / Foam::pow(nuF, 2),
        1.0/3.0
    );
Info << "dimensionless diameter, Dstar = " << Dstar.value() << endl;

// maximum sediment volumic fraction (1 - porosity)
dimensionedScalar CsMax("CsMax", dimless, transportProperties);
Info << "maximum possible sediment volume fraction: "
    << CsMax.value() << endl;
// Schmidt number, default value 1
dimensionedScalar Sc("Sc", dimless, 1);
if (transportProperties.found("Sc"))
{
    Sc.value() = readScalar(transportProperties.lookup("Sc"));
}

Info<< "Reading suspensionProperties\n" << endl;
IOdictionary suspensionProperties
(
    IOobject
    (
        "suspensionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Check if suspension is on or off
const word switchSuspension(suspensionProperties.lookup("suspension"));
if (switchSuspension!="on" and switchSuspension!="off")
{
    FatalError << switchSuspension
        << " is not a valid entry for suspension"
        << ", valid options are: on off" << endl;
    Info << abort(FatalError) << endl;
}
Info << "suspension is " << switchSuspension << endl;

word switchMULES = "off";
if (suspensionProperties.found("MULESformulation"))
{
    switchMULES = word(suspensionProperties.lookup("MULESformulation"));
}
if (switchMULES!="on" and switchMULES!="off")
{
    FatalError << switchMULES
        << " is not a valid entry for MULESformulation"
        << ", valid options are: on off" << endl;
    Info << abort(FatalError) << endl;
}

// use nut/Sc or nu + nut/Sc for diffusivity, options are effective/turbulent
const word CsDiffusivity =
    suspensionProperties.lookupOrDefault<word>
    (
        "diffusivity",
        "turbulent"
    );
if (CsDiffusivity != "turbulent" and CsDiffusivity != "effective")
{
    FatalError << CsDiffusivity
        << " is not a valid entry for diffusivity"
        << ", valid options are: turbulent effective" << endl;
    Info << abort(FatalError) << endl;
}

// add wall diffusivity for cells close to bed patches
bool isWallDiff = false;
dimensionedScalar epsWall0(dimViscosity, 0);
dimensionedScalar betaWall(dimless, 0);
dimensionedScalar gammaWall(dimless, 2.5);
if (suspensionProperties.found("wallDiffusivity"))
{
    isWallDiff = true;
    const dictionary& wallDiffDict
        (
            suspensionProperties.subDict("wallDiffusivity")
        );
    epsWall0.value() = wallDiffDict.lookupOrDefault<scalar>("alpha", 0);
    betaWall.value() = wallDiffDict.lookupOrDefault<scalar>("beta", 5);
    gammaWall.value() = wallDiffDict.lookupOrDefault<scalar>("gamma", 2.5);
}
if (isWallDiff)
{
    Info << "wall diffusivity (saltation) enabled, "
        << "alpha = " << epsWall0.value()
        << ", beta = " << betaWall.value()
        << ", gamma = " << gammaWall.value()
        << endl;
}
else
{
    Info << "wall diffusivity off" << endl;
}
// Nikuradse equivalent roughness height, skin friction
dimensionedScalar ks =  gammaWall * dS;

// maximum reference concentration for erosion condition
// sediment vol fraction, default is 0.3 , user modifiable
dimensionedScalar CsbMax(dimless, 0.3);
if (suspensionProperties.found("CsbMax"))
{
    scalar csbmaxVal = readScalar(suspensionProperties.lookup("CsbMax"));
    CsbMax.value() = Foam::min(CsMax.value(), csbmaxVal);
}
Info << "maximum reference sediment volumic fraction: "
    << CsbMax.value() << endl;

// Instantiate pointers to fields associated with suspension
autoPtr<volScalarField> CsPtr;  // suspended load concentration
autoPtr<volVectorField> WsPtr;  // sediment settling velocity
autoPtr<volScalarField> diffSedPtr;  // turbulent diffusivity for sediment
autoPtr<volScalarField> erosionPtr;  // erosion flux
autoPtr<volScalarField> depositionPtr;  // deposition flux
autoPtr<surfaceScalarField> phipPtr;  // flux of suspended sediment
autoPtr<wallDist> bedDistPtr;  // wallDist for bed patches

autoPtr<settlingModel> settlingModelPtr;

if (switchSuspension=="on")
{
    Info<< "Reading field Cs" << endl;
    CsPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "Cs",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    Info<< "Reading field Ws\n" << endl;
    WsPtr.reset
    (
        new volVectorField
        (
            IOobject
            (
                "Ws",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    volScalarField& Cs = CsPtr.ref();  // suspended load
    volVectorField& Ws = WsPtr.ref();  // settling velocity

    // Initialize Settling Model
    Info << "Initialize settling Model\n" << endl;
    settlingModelPtr.reset(new settlingModel(suspensionProperties));

    Info<< "Initialize sediment settling velocity\n" << endl;
    Ws = settlingModelPtr->Ufall
        (
            Cs,
            CsMax,
            dS,
            rhoS,
            rhoF,
            nuF,
            mag(g)
        ) *g/mag(g);

    Ws.correctBoundaryConditions();

    erosionPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "Erosion",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimensionSet(0, 0, -1, 0, 0, 0, 0), Zero)
        )
    );

    depositionPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "Deposition",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimensionSet(0, 0, -1, 0, 0, 0, 0), Zero)
        )
    );

    Info << "Reading/calculating face flux field phi\n" << endl;
    phipPtr.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phip",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(U+Ws)
            //fvc::interpolate(U+Ws) & mesh.Sf()
        )
    );

    Info << "Initialize turbulent diffusivity for sediment\n" << endl;
    diffSedPtr.reset(
        new volScalarField
        (
            IOobject
            (
                "diffSed",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimViscosity, Zero)
        )
    );
}

//phip = fvc::flux(U+Ws);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

volVectorField shieldsVf
(
    IOobject
    (
        "shieldsVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimless, Zero)
);

volScalarField critShieldsVf
(
    IOobject
    (
        "critShieldsVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

volVectorField qsatVf
(
    IOobject
    (
        "qsatVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(0, 2, -1, 0, 0, 0, 0), Zero)
);

volVectorField qavVf
(
    IOobject
    (
        "qavVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(0, 2, -1, 0, 0, 0, 0), Zero)
);

volVectorField qbVf
(
    IOobject
    (
        "qbVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(0, 2, -1, 0, 0, 0, 0), Zero)
);

volScalarField betaVf
(
     IOobject
    (
        "betaVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

volScalarField betaTildeVf
(
     IOobject
    (
        "betaTildeVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

volVectorField VelBfVf
(
    IOobject
    (
        "VelBfVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(0, 2, -1, 0, 0, 0, 0), Zero)
);

#include "createMRF.H"
#include "createFvOptions.H"
