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

Description
    Save, at the first PIMPLE outer iteration, the bed state at the beginning
    of the time step. Every sub-iteration PIMPLE <-> Exner restarts from it.
    Must be included before the first mesh update of the time step.

\*---------------------------------------------------------------------------*/

// sub-iteration counter of the current time step
zbIter = 0;

if (switchMULES == "on" && switchSuspension == "on")
{
    CsStart = CsPtr.ref().primitiveField();

Info<< "CsStart min = " << gMin(CsStart) << endl;
Info<< "CsStart max = " << gMax(CsStart) << endl;
}

if (bed.exist())
{
    // bedload at t: explicit / semiImplicit saturation is an evolution equation
    // from qb(t), it must not be advanced once per sub-iteration
    qbStart = qbPtr.ref().primitiveField();

    // bed points displacement at t: the displacement imposed at each
    // sub-iteration is dispVals0 - eg*dH, not cumulative
//if (bed.bedMotion())
//{
//    const pointVectorField& pointDisplacement =
//  mesh.objectRegistry::lookupObject<pointVectorField>("pointDisplacement");
//
//    // Sauvegarder TOUT le champ, pas juste le dernier patch
//    dispVals0 = pointDisplacement.primitiveField();
//}

    if (bed.bedMotion())
    {
        const pointScalarField& pointDisplacement =
            mesh.objectRegistry::lookupObject<pointScalarField>
            (
                "pointDisplacementz"
            );

        for (const label patchID : bed.bedPatchesID())
        {
            dispVals0 = //pointDisplacement.primitiveField();
                pointDisplacement.boundaryField()[patchID].patchInternalField();
        }
    }
}
