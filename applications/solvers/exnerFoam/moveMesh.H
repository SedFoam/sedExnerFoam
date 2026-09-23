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

pointVectorField& pointDisplacement = const_cast<pointVectorField&>
(
    mesh.objectRegistry::lookupObject<pointVectorField>
    (
        "pointDisplacement"
    )
);

for (const label patchID : bed.bedPatchesID())
{
    dispVals =
        pointDisplacement.boundaryFieldRef()[patchID].patchInternalField();
}

dzbPoints = Zero;

//interpolate bed level increment from face center to vertices
bed.interpFaceToVertices(dzb, dzbPoints);

dispVals -= eg * dzbPoints;

for (const label patchID : bed.bedPatchesID())
{
    Info << "mesh vertices motion, max:" << gMax(dispVals & -eg);
    Info << " m , min:" << gMin(dispVals & -eg);
    Info << " m , mean:" << gAverage(-dispVals & eg) << " m" << endl;
    pointDisplacement.boundaryFieldRef()[patchID] == dispVals;
}

pointDisplacement.correctBoundaryConditions();
