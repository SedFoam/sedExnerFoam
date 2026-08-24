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

pointVectorField& pointDisplacement = const_cast<pointVectorField&>
(
    mesh.objectRegistry::lookupObject<pointVectorField>
    (
        "pointDisplacement"
    )
);

vectorField& dispVals = dispValsPtr.ref();
for (const label patchID : bed.bedPatchesID())
{
    dispVals =
        pointDisplacement.boundaryFieldRef()[patchID].patchInternalField();
}

scalarField& dHpoints = dHpointsPtr.ref();
dHpoints = Zero;

//interpolate bed level increment from face center to vertices
bed.interpFaceToVertices(dH, dHpoints);

// zero displacements for faPatches
const labelList& fixedPointsLabels = fixedPatchPointLabelsPtr.ref();

for (const label i : fixedPointsLabels)
{
    dHpoints[i] = 0.;
}

dispVals -= eg * dHpoints;


// interpolate displacement on points, vertical orientation
//dispVals -= patchInterp.faceToPointInterpolate(dH) * eg;

for (const label patchID : bed.bedPatchesID())
{
    Info << "mesh vertices motion, max:" << gMax(dispVals & -eg);
    Info << " m , min:" << gMin(dispVals & -eg);
    Info << " m , mean:" << gAverage(-dispVals & eg) << " m" << endl;
    pointDisplacement.boundaryFieldRef()[patchID] == dispVals;
    pointDisplacement.correctBoundaryConditions();
}
