/*---------------------------------------------------------------------------*\
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of exnerFOAM.

    exnerFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    exnerFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with exnerFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

Info<< "Reading bedloadProperties\n" << endl;

IOdictionary bedloadProperties
(
    IOobject
    (
        "bedloadProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Instantiate sedimentBed object
Info << "Instantiate sedimentBed" << endl;
sedimentBed bed(bedloadProperties, mesh, g);

// Const refrence to finite area mesh
const faMesh& aMesh = bed.aMesh();

// Const reference to projected finite area mesh
const projectedFaMesh& meshProj = bed.aProjMesh();

// bed elevation increment on points
vectorField dispVals(aMesh.nPoints(), Zero);
scalarField dHpoints(aMesh.nPoints(), Zero);

const dictionary bedloadModelDict = bedloadProperties.subDict("bedloadModel");

// Read alpha and beta coefficients for bedload flux
//dimensionedScalar alphaQb("alpha", bedloadModelDict);
scalar alphaQb(bedloadModelDict.get<scalar>("alpha"));
scalar betaQb(bedloadModelDict.get<scalar>("beta"));
Info << "bedload transport law, qb = alpha*U^beta,"
    << "alpha = " << alphaQb << ", beta = " << betaQb << endl;

// divergence scheme for qb, linear, upwind or linearUpwind
word qbScheme(bedloadProperties.get<word>("qbScheme"));
Info << "div(qb) scheme: " << qbScheme << endl;

// Order of explicit Adams-Bashforth time scheme
// default is 1, euler explicit, options are 1 and 2
label ABorder(bedloadProperties.get<label>("ABorder"));

// points displacements values for mesh motion
autoPtr<vectorField> dispValsPtr;

scalar betaRepDeg =
(
    transportProperties.lookupOrDefault<scalar>("reposeAngle", 32)
);

Info << "granular material repose angle: " << betaRepDeg
    << " degrees" << endl;
if (betaRepDeg <= 0 or betaRepDeg >= 90)
{
    FatalError
        << "repose angle must be in the interval 0 - 90" << endl;
    Info << abort(FatalError) << endl;
}

scalarField zb(aMesh.nFaces(), Zero);

areaVectorField Q
(
    IOobject
    (
        "Q",
        runTime.timeName(),
        aMesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    aMesh
);

areaScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        aMesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    aMesh
);

areaVectorField gradZb
(
    IOobject
    (
        "gradZb",
        runTime.timeName(),
        aMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedVector(dimless, Zero)
);

areaScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        aMesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh,
    dimensionedScalar(dimless, Zero)
);

areaScalarField dVsed
(
    IOobject
    (
        "dVsed",
        runTime.timeName(),
        aMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedScalar(dimensionSet(0, 3, -1, 0, 0, 0, 0), Zero)
);

areaVectorField qb
(
    IOobject
    (
        "qb",
        runTime.timeName(),
        aMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh
);

scalarField phiqbProj(aMesh.nEdges(), 0);

tensorField gradQbProj(aMesh.nFaces(), Zero);

areaScalarField dH
(
    IOobject
    (
        "deltaH",
        runTime.timeName(),
        aMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedScalar(dimLength, Zero)
);
