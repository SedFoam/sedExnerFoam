/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    suspensionFoam

Group
    grpBasicSolvers

Description
    Scalar transport and incompressible turbulent flow solver.

\*---------------------------------------------------------------------------*/
//    \heading Solver details
//    The equation is given by:
//
//    \f[
//        \ddt{T} + \div \left(\vec{U} T\right) - \div \left(D_T \grad T \right)
//        = S_{T}
//    \f]

//    Where:
//    \vartable
//        Cs         | Volumic fraction
//        \vec{U}    | Velocity
//        \vec{R}    | Stress tensor
//        p          | Pressure
//        \vec{S}_U  | Momentum source
//    \endvartable

//    \heading Required fields
//    \plaintable
//        Cs      | Passive scalar
//        U       | Velocity [m/s]
//        p       | Kinematic pressure, p/rho [m2/s2]
//        \<turbulence fields\> | As required by user selection
//    \endplaintable


#include "fvCFD.H"
#include "faCFD.H"
#include "CMULES.H"
#include "dynamicFvMesh.H"
#include "meshTools.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "unitConversion.H"
#include "fvcSmooth.H"
#include "primitivePatchInterpolation.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "wallDist.H"

#include "MatrixTools.H"
#include "LUscalarMatrix.H"

#include "settlingModel.H"
#include "criticalShieldsModel.H"
#include "bedloadModel.H"
#include "sedimentBed.H"
#include "projectedFaMesh.H"
#include "filter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Scalar transport and incompressible turbulent flow solver."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFaFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();
    
    if (not LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                        if (switchSuspension=="on")
                        {
                            surfaceScalarField& phip = phipPtr.ref();
                            fvc::makeRelative(phip, U);
                        }
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                    if (bed.rigidBed())
                    {
                        areaVectorField& rigidBed = rigidBedPtr.ref();
                        areaScalarField& HsedBed = HsedBedPtr.ref();

                        HsedBed =
                            (
                                rigidBed
                                - bed.aMesh().areaCentres()
                            ) & (g/mag(g)).value();

                        Info << "width of sediment layer, max: "
                            << max(HsedBed).value() << " m, min: "
                            << min(HsedBed).value() << " m, mean: "
                            << average(HsedBed).value() << endl;
                    }
                }
            }
            
            #include "UEqn.H"
            
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        if (bed.exist())
        {
            #include "bedShearStress.H"

            #include "bedload.H"
        }

        if (switchSuspension=="on")
        {
            #include "CsEqn.H"
        }

        if (bed.exist())
        {
            #include "exnerEqn.H"
        }

        if (runTime.writeTime())
        {
            if (bed.exist())
            {
                areaVectorField& shields = shieldsPtr.ref();
                areaScalarField& critShields = critShieldsPtr.ref();
                areaVectorField& qsat = qsatPtr.ref();
                areaVectorField& qav = qavPtr.ref();
                areaVectorField& qb = qbPtr.ref();
                areaScalarField& beta = betaPtr.ref();
                
                // map areaFields to volFields for vizualisation
                bed.vsm.ref().mapToVolume
                    (
                        shields,
                        shieldsVf.boundaryFieldRef()
                    );
                bed.vsm.ref().mapToVolume
                    (
                        critShields,
                        critShieldsVf.boundaryFieldRef()
                    );
                bed.vsm.ref().mapToVolume
                    (
                        qsat,
                        qsatVf.boundaryFieldRef()
                    );
                bed.vsm.ref().mapToVolume
                    (
                        qav,
                        qavVf.boundaryFieldRef()
                    );
                bed.vsm.ref().mapToVolume
                    (
                        qb,
                        qbVf.boundaryFieldRef()
                    );
                bed.vsm.ref().mapToVolume
                    (
                        beta,
                        betaVf.boundaryFieldRef()
                    );
            }
            runTime.write();

            runTime.printExecutionTime(Info);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
