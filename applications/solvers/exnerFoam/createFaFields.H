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

// bed elevation increment on points
vectorField dispVals(bed.aMesh().nPoints(), Zero);
scalarField dzbPoints(bed.aMesh().nPoints(), Zero);

const dictionary bedloadModelDict = bedloadProperties.subDict("bedloadModel");

// Read alpha and beta coefficients for bedload flux
//dimensionedScalar alphaQb("alpha", bedloadModelDict);
scalar alphaQb(bedloadModelDict.get<scalar>("alpha"));
scalar betaQb(bedloadModelDict.get<scalar>("beta"));

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
scalar betaRep = degToRad(betaRepDeg);

areaVectorField Q
(
    IOobject
    (
        "Q",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    bed.aMesh()
);

areaScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    bed.aMesh()
);

areaScalarField zb
(
    IOobject
    (
        "zb",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    bed.aMesh(),
    dimensionedScalar(dimLength, Zero)
);

areaScalarField dzb
(
    IOobject
    (
        "dzb",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    bed.aMesh(),
    dimensionedScalar(dimLength, Zero)
);

areaVectorField gradZb
(
    IOobject
    (
        "gradZb",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    bed.aMesh(),
    dimensionedVector(dimless, Zero)
);

areaScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    bed.aMesh(),
    dimensionedScalar(dimless, Zero)
);

areaVectorField qb
(
    IOobject
    (
        "qb",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    bed.aMesh()
);

edgeScalarField phiqb
(
    IOobject
    (
        "phiqb",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fac::interpolate(qb) & bed.aMesh().Le()
);

areaVectorField qav
(
    IOobject
    (
        "qav",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    bed.aMesh(),
    dimensionedVector(dimensionSet(0, 2, -1, 0, 0, 0, 0), Zero)
);

edgeScalarField phiqav
(
    IOobject
    (
        "phiqav",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fac::interpolate(qav) & bed.aMesh().Le()
);

areaScalarField dH
(
    IOobject
    (
        "deltaH",
        runTime.timeName(),
        bed.aMesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    bed.aMesh(),
    dimensionedScalar(dimLength, Zero)
);
