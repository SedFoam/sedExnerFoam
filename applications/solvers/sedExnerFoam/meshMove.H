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

pointScalarField& pointDisplacement = const_cast<pointScalarField&>
(
    mesh.objectRegistry::lookupObject<pointScalarField>
    (
        "pointDisplacementz"
    )
);

// The imposed displacement is rebuilt at every PIMPLE sub-iteration from the
// displacement at the beginning of the time step (dispVals0, saved in
// storeBedState.H) and the current dH. Starting from the current
// pointDisplacement would add dH once per sub-iteration.
scalarField& dispVals = dispValsPtr.ref();
scalarField& dHpoints = dHpointsPtr.ref();

// Initialiser dispVals avec les déplacements initiaux (corrigés)
dispVals = dispVals0;


// Forcer à zéro sur TOUS les points qui ne sont PAS des patches du lit
forAll(pointDisplacement.boundaryField(), patchi)
{
    // Si ce n'est PAS un patch du lit, forcer à zéro
    if (findIndex(bed.bedPatchesID(), patchi) == -1)
    {
        pointDisplacement.boundaryFieldRef()[patchi] = 0;
    }
}

// Maintenant calculer le nouveau déplacement
dHpoints = Zero;
//interpolate bed level increment from face center to vertices
bed.interpFaceToVertices(dH, dHpoints);

// zero displacements for faPatches
const labelList& fixedPointsLabels = fixedPatchPointLabelsPtr.ref();

for (const label i : fixedPointsLabels)
{
    dHpoints[i] = 0.;
}

dispVals += dHpoints;

// interpolate displacement on points, vertical orientation
//dispVals -= patchInterp.faceToPointInterpolate(dH) * eg;

// Appliquer uniquement aux patches du lit
for (const label patchID : bed.bedPatchesID())
{
    Info << "mesh vertices motion, max:" << gMax(-dispVals);
    Info << " m , min:" << gMin(-dispVals);
    Info << " m , mean:" << gAverage(-dispVals) << " m" << endl;
    pointDisplacement.boundaryFieldRef()[patchID] == dispVals;
}
pointDisplacement.correctBoundaryConditions();
