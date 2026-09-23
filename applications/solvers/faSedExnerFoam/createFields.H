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

autoPtr<settlingModel> settlingModelPtr;

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

#include "createMRF.H"
#include "createFvOptions.H"
