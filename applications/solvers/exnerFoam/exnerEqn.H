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

// area of projected mesh faces
const scalarField& SfProj = meshProj.S();

#include "balanceBedload.H"

scalar deltaT = runTime.time().deltaT().value();  // current time step

// solve exner equation explicitly with Adams-Bashforth scheme
if (ABorder==1)
{
    // Adams-Bashforth order 1, euler explicit scheme
    forAll(aMesh.areaCentres(), facei)
    {
        dH[facei] =
            (deltaT / SfProj[facei])
            * dVsed[facei];
    }
}
else if (ABorder==2)
{
    scalar deltaT0 = runTime.time().deltaT0().value();  // old time step

    scalar coefT0 = deltaT * (1 + 0.5 * (deltaT / deltaT0));
    scalar coefT00 = - 0.5 * deltaT * deltaT / deltaT0;

    // Adams-Bashforth order 1
    forAll(aMesh.areaCentres(), facei)
    {
        dH[facei] =
            (coefT0 * dVsed[facei] + coefT00 * dVsed.oldTime()[facei])
            * (1 / SfProj[facei]);
    }
}
else
{
    FatalError
        << "Invalid option " << ABorder << " for ABorder" << endl
        << "valid options are: 1 2" << endl;
    Info << abort(FatalError) << endl;
}

#include "moveMesh.H"
