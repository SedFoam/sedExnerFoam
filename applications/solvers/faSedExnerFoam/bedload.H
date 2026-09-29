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

// access to bedload field
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
qb = bedloadMod.qb(
    shields, critShields, gradZb, rhoS, rhoF, mag(g), dS, nuF);

qb.correctBoundaryConditions();

Info << "Avalanche calculation" << endl;
// additional bedload due to avalanche
if (bedloadMod.avalanche())
{
    #include "sedimAvalanche.H"
}
