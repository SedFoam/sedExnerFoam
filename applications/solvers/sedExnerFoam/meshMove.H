/*---------------------------------------------------------------------------*\
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of sedExnerFOAM.

    sedExnerFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    sedExnerFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ScourFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Move the bed mesh points at each PIMPLE <-> Exner sub-iteration.

    The imposed displacement of the bed patches is rebuilt from the
    displacement at the beginning of the time step (dispVals0, saved in
    storeBedState.H) plus the bed level increment dH interpolated from face
    centres to vertices, so that it is not cumulative across sub-iterations.

    Two mesh motion solvers are supported, detected from the displacement
    field registered on the mesh:
    - displacementComponentLaplacian : pointScalarField
      "pointDisplacementz" (vertical component only)
    - displacementLaplacian          : pointVectorField
      "pointDisplacement" (displacement applied along eg)

    The displacement is forced to zero on all non-bed patches, and set to
    zero on the points of the faPatches (fixedPointsLabels).
    Boundary conditions are corrected on the field used by the solver.

\*---------------------------------------------------------------------------*/

scalarField& dispVals = dispValsPtr.ref();
scalarField& dHpoints = dHpointsPtr.ref();

// Initiate dispVals with initialdisplacements
dispVals = dispVals0;

// -------------------------------------------------------------------
// Detect the mesh motion solver from the displacement field it provides
// - displacementComponentLaplacian (component z) -> pointScalarField
//   "pointDisplacementz"
// - displacementLaplacian                          -> pointVectorField
//   "pointDisplacement"
// -------------------------------------------------------------------
const bool componentLaplacian =
    mesh.foundObject<pointScalarField>("pointDisplacementz");

// Force to zero for all patches not cooresponding to sediment bed
if (componentLaplacian)
{
    pointScalarField& pointDispz =
        mesh.lookupObjectRef<pointScalarField>("pointDisplacementz");

    forAll(pointDispz.boundaryField(), patchi)
    {
        // force to 0 if patch is not a bedPatch
        if (not bed.bedPatchesID().found(patchi))
        {
            pointDispz.boundaryFieldRef()[patchi] == scalar(0);
        }
    }
}
else
{
    pointVectorField& pointDisp =
        mesh.lookupObjectRef<pointVectorField>("pointDisplacement");

    forAll(pointDisp.boundaryField(), patchi)
    {
        // force to 0 if patch is not a bedPatch
        if (not bed.bedPatchesID().found(patchi))
        {
            pointDisp.boundaryFieldRef()[patchi] == vector(0, 0, 0);
        }
    }
}

// Maintenant calculer le nouveau déplacement
dHpoints = Zero;
// interpolate bed level increment from face center to vertices
bed.interpFaceToVertices(dH, dHpoints);

// zero displacements for faPatches
const labelList& fixedPointsLabels = fixedPatchPointLabelsPtr.ref();

for (const label i : fixedPointsLabels)
{
    dHpoints[i] = 0.;
}

dispVals += dHpoints;

// Appliquer uniquement aux patches du lit
for (const label patchID : bed.bedPatchesID())
{
    Info << "mesh vertices motion, max:" << gMax(-dispVals);
    Info << " m , min:" << gMin(-dispVals);
    Info << " m , mean:" << gAverage(-dispVals) << " m" << endl;

    if (componentLaplacian)
    {
        pointScalarField& pointDispz =
            mesh.lookupObjectRef<pointScalarField>("pointDisplacementz");

        pointDispz.boundaryFieldRef()[patchID] == dispVals;
    }
    else
    {
        pointVectorField& pointDisp =
            mesh.lookupObjectRef<pointVectorField>("pointDisplacement");

        // vertical displacement only
        pointDisp.boundaryFieldRef()[patchID] == -dispVals * eg;
    }
}

// Correct boundary conditions on the field actually used by the solver
if (componentLaplacian)
{
    mesh.lookupObjectRef<pointScalarField>("pointDisplacementz")
        .correctBoundaryConditions();
}
else
{
    mesh.lookupObjectRef<pointVectorField>("pointDisplacement")
        .correctBoundaryConditions();
}
